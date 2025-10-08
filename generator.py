"""
Decoy molecule generator.

This module contains the main decoy generation logic, including the
DecoysGenerator class and related processing functions.
"""

import os
import pandas as pd
import numpy as np
import time
from multiprocessing import Pool
from numpy.lib.format import open_memmap
from typing import Dict, Tuple, Any

from models import Molecule
from database import DatabaseProperties
from utils import get_charge, predict_charge, sim_filter


def process_molecule_batch(args: Tuple) -> Tuple[str, Molecule]:
    """
    Process a single molecule batch without class dependencies to avoid pickling issues.
    
    This function calculates similarity scores between a target molecule and 
    database molecules, filters by threshold and charge, and samples decoys.
    
    Args:
        args: Tuple containing (name, molecule, df, arr, config, db_props, n_oversampled_decoys)
        
    Returns:
        Tuple[str, Molecule]: (molecule_name, updated_molecule)
    """
    name, molecule, df, arr, config, db_props, n_oversampled_decoys = args
    
    # Calculate similarity scores without copying the dataframe yet
    scores = (
        (np.absolute(arr - molecule.properties) / db_props.queried_enamine_std)
        .sum(axis=1)
    )
    
    # Filter by threshold first to reduce memory usage
    threshold_mask = scores < molecule.threshold
    
    # Only proceed if we have molecules that pass the threshold
    if not threshold_mask.any():
        return name, molecule  # No qualifying molecules, return unchanged
    
    # Now create a copy only of the filtered subset
    df_copy = df[threshold_mask].copy()
    df_copy['score'] = scores[threshold_mask]
    
    # Add molecular properties to the dataframe
    filtered_arr = arr[threshold_mask]
    active_property_names = db_props.get_active_property_names()
    
    # Add each property as a column
    for i, prop_name in enumerate(active_property_names):
        df_copy[prop_name] = filtered_arr[:, i]
    
    # Filter by charge if enabled
    if config.get('Charge', False):
        if config.get('Protonate', False):
            df_copy = df_copy[
                df_copy.smi.apply(lambda x: predict_charge(x) == molecule.charge)
            ]
        else:
            df_copy = df_copy[
                df_copy.smi.apply(lambda x: get_charge(x) == molecule.charge)
            ]
    
    # Combine with existing decoys
    if molecule.decoys.empty:
        # No existing decoys, use the new dataframe as-is
        combined_df = df_copy
    else:
        # Add missing property columns to existing decoys if they don't exist
        existing_decoys = molecule.decoys.copy()
        for prop_name in active_property_names:
            if prop_name not in existing_decoys.columns:
                # Calculate properties for existing decoys
                existing_decoys[prop_name] = existing_decoys['smi'].apply(
                    lambda smi: db_props.calculate_queried_properties(smi)[
                        active_property_names.index(prop_name)
                    ]
                )
        
        combined_df = pd.concat([existing_decoys, df_copy], ignore_index=True)
    
    df_copy = combined_df
    
    # Sample if we have too many decoys
    if df_copy.shape[0] > n_oversampled_decoys:
        # Take the best (lowest) scores instead of random sampling
        df_copy = df_copy.nsmallest(n_oversampled_decoys, 'score')
        # Update threshold to the maximum score in the selected set
        molecule.threshold = df_copy.score.max()
    
    # Update molecule with new decoys
    molecule.decoys = df_copy
    
    return name, molecule


class DecoysGenerator:
    """
    Main class for generating molecular decoys.
    
    This class orchestrates the decoy generation process by loading target molecules,
    processing database bundles, and generating decoys based on molecular similarity.
    """
    
    def __init__(self, enamine_std_path: str, prop_bundle_dir: str, smi_bundle_dir: str,
                 prop_prefix: str, smi_prefix: str, config_file: str, ligand_smi_file: str,
                 oversample: int = 10, n_proc: int = 8, batch_size: int = 10000, verbose: bool = False,
                 apply_similarity_filter: bool = True):
        """
        Initialize the decoy generator.
        
        Args:
            enamine_std_path: Path to Enamine standard deviations file
            prop_bundle_dir: Directory containing property bundles
            smi_bundle_dir: Directory containing SMILES bundles
            prop_prefix: Prefix for property bundle files
            smi_prefix: Prefix for SMILES bundle files
            config_file: Path to configuration file
            ligand_smi_file: Path to file containing target ligand SMILES
            oversample: Oversampling factor for decoy generation
            n_proc: Number of parallel processes
            batch_size: Batch size for processing database molecules
            verbose: Enable verbose output with timing information
            apply_similarity_filter: Whether to apply similarity filtering after each bundle
        """
        self.config_file = config_file
        self.ligand_smi_file = ligand_smi_file
        self.oversample = oversample
        self.n_proc = n_proc
        self.batch_size = batch_size
        self.verbose = verbose
        self.apply_similarity_filter = apply_similarity_filter
        
        # Initialize database handler
        self.db_props = DatabaseProperties(
            enamine_std_path, prop_bundle_dir, smi_bundle_dir,
            prop_prefix, smi_prefix, config_file
        )
        
        # Load configuration
        self.config = self.db_props.config
        
        # Calculate target number of oversampled decoys
        self.n_oversampled_decoys = self.config['Max_decoys_per_ligand'] * self.oversample
        
        # Load target molecules
        self.actives = self._load_ligands()
    
    def _load_ligands(self) -> Dict[str, Molecule]:
        """
        Load target ligand molecules from file.
        
        Returns:
            Dict[str, Molecule]: Dictionary mapping molecule names to Molecule objects
            
        Raises:
            FileNotFoundError: If the ligand file doesn't exist
            ValueError: If the file format is invalid
        """
        if not os.path.exists(self.ligand_smi_file):
            raise FileNotFoundError(f"Ligand file not found: {self.ligand_smi_file}")
        
        actives = {}
        
        try:
            with open(self.ligand_smi_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue  # Skip empty lines and comments
                    
                    parts = line.split()
                    if len(parts) != 2:
                        raise ValueError(
                            f"Invalid format on line {line_num}. Expected: 'SMILES NAME'"
                        )
                    
                    smiles, name = parts
                    
                    # Calculate charge and properties
                    if self.config.get('Protonate', False):
                        charge = predict_charge(smiles)  # pH-dependent protonation
                    else:
                        charge = get_charge(smiles)  # Use SMILES as-is
                    all_properties = self.db_props.calculate_all_properties(smiles)
                    queried_props = all_properties[self.db_props.column_mask]
                    
                    # Calculate threshold
                    threshold = self.config['Threshold'] * len(queried_props)
                    
                    # Create molecule object
                    molecule = Molecule(
                        smiles=smiles,
                        charge=charge,
                        properties=queried_props,
                        threshold=threshold,
                        decoys=pd.DataFrame()
                    )
                    
                    actives[name] = molecule
                    
        except Exception as e:
            raise ValueError(f"Error loading ligands from {self.ligand_smi_file}: {e}")
        
        if not actives:
            raise ValueError(f"No valid ligands found in {self.ligand_smi_file}")
        
        return actives
    
    def process_bundle(self, bundle_index: int) -> None:
        """
        Process a single database bundle with optimized file reading.
        
        Performance optimizations by database fraction:
        - fraction >= 1.0: Fast path, no sampling overhead
        - 0.1 < fraction < 1.0: Chunked random sampling (50K chunks)
        - 0.01 < fraction <= 0.1: Geometric skip patterns for efficiency
        - fraction <= 0.01: Ultra-fast seeking method (avoids reading unwanted lines)
        
        For 2B+ line files with small fractions, the seeking method can be 
        100-1000x faster than sequential reading.
        
        Args:
            bundle_index: Index of the bundle to process
            
        Raises:
            FileNotFoundError: If bundle files don't exist
        """
        # Get bundle filenames
        smi_file, prop_file = self.db_props.get_bundle_filenames(bundle_index)
        
        # Validate files exist
        if not self.db_props.validate_bundle_files(bundle_index):
            raise FileNotFoundError(f"Bundle files not found for index {bundle_index}")
        
        # Load property array as memory map
        prop_arr = self.db_props.load_property_bundle(bundle_index)
        
        # Process SMILES file in batches
        database_fraction = self.config.get('Database_fraction', 1.0)
        
        # For extremely large files (>1B lines), consider line-counting approach
        if database_fraction <= 0.01:  # Less than 1% sampling - use ultra-fast seeking
            return self._process_bundle_with_seeking(smi_file, prop_arr, bundle_index)
        
        with open(smi_file, 'r') as f:
            df_batch = []
            selected_indices = []
            total_processed = 0
            
            if database_fraction >= 1.0:
                # Fast path: use all lines
                for line_count, line in enumerate(f):
                    # Parse SMILES line
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        df_batch.append(parts[:2])
                        selected_indices.append(total_processed)
                    elif len(parts) == 1:
                        df_batch.append([parts[0], f"mol_{line_count}"])
                        selected_indices.append(total_processed)
                    else:
                        total_processed += 1
                        continue
                    
                    total_processed += 1
                    
                    # Process batch when it reaches batch_size
                    if len(df_batch) >= self.batch_size:
                        if df_batch:
                            self._process_batch(df_batch, prop_arr, selected_indices)
                        df_batch = []
                        selected_indices = []
            else:
                # Ultra-optimized for very large files (2B+ lines)
                # Use deterministic sampling for reproducibility
                np.random.seed(hash(smi_file) % (2**32))  # Seed based on filename
                
                if database_fraction <= 0.1:  # For small fractions, use line skipping
                    # Calculate average lines to skip between samples
                    skip_avg = int(1.0 / database_fraction) if database_fraction > 0 else 1000000
                    
                    # Pre-generate skip patterns for reproducibility
                    skip_pattern_size = 1000
                    skip_offsets = np.random.geometric(database_fraction, skip_pattern_size)
                    skip_idx = 0
                    
                    line_count = 0
                    next_sample_line = skip_offsets[0]
                    
                    for line in f:
                        line_count += 1
                        
                        # Skip until we reach the next sample line
                        if line_count < next_sample_line:
                            total_processed += 1
                            continue
                        
                        # Update next sample line
                        skip_idx = (skip_idx + 1) % skip_pattern_size
                        next_sample_line += skip_offsets[skip_idx]
                        
                        # Parse SMILES line
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            df_batch.append(parts[:2])
                            selected_indices.append(total_processed)
                        elif len(parts) == 1:
                            df_batch.append([parts[0], f"mol_{line_count-1}"])
                            selected_indices.append(total_processed)
                        else:
                            total_processed += 1
                            continue
                        
                        total_processed += 1
                        
                        # Process batch when it reaches batch_size
                        if len(df_batch) >= self.batch_size:
                            if df_batch:
                                self._process_batch(df_batch, prop_arr, selected_indices)
                            df_batch = []
                            selected_indices = []
                else:
                    # For larger fractions (>10%), use chunked random sampling
                    # Generate sampling decisions in larger chunks for better performance
                    chunk_size = 50000  # Larger chunks for better performance
                    sampling_chunk = np.random.random(chunk_size) <= database_fraction
                    chunk_idx = 0
                    
                    for line_count, line in enumerate(f):
                        # Refresh sampling chunk if needed
                        if chunk_idx >= chunk_size:
                            sampling_chunk = np.random.random(chunk_size) <= database_fraction
                            chunk_idx = 0
                        
                        # Skip line based on pre-generated sampling decision
                        if not sampling_chunk[chunk_idx]:
                            chunk_idx += 1
                            total_processed += 1
                            continue
                        
                        chunk_idx += 1
                        
                        # Parse SMILES line
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            df_batch.append(parts[:2])
                            selected_indices.append(total_processed)
                        elif len(parts) == 1:
                            df_batch.append([parts[0], f"mol_{line_count}"])
                            selected_indices.append(total_processed)
                        else:
                            total_processed += 1
                            continue
                        
                        total_processed += 1
                        
                        # Process batch when it reaches batch_size
                        if len(df_batch) >= self.batch_size:
                            if df_batch:
                                self._process_batch(df_batch, prop_arr, selected_indices)
                            df_batch = []
                            selected_indices = []
            
            # Process any remaining batch
            if df_batch:
                self._process_batch(df_batch, prop_arr, selected_indices)
    
    def _process_bundle_with_seeking(self, smi_file: str, prop_arr: np.ndarray, bundle_index: int) -> None:
        """
        Ultra-fast processing for very small database fractions using true file seeking.
        
        For extremely small fractions (≤1%), this method:
        1. Estimates file size and line positions using byte offsets
        2. Pre-calculates exact byte positions to seek to
        3. Seeks directly to those positions without reading intermediate lines
        4. Perfect for 0.01% sampling on 2B line files (reads only ~200K lines)
        
        Args:
            smi_file: Path to SMILES file
            prop_arr: Property array for the current bundle
            bundle_index: Bundle index for progress reporting
        """
        database_fraction = self.config.get('Database_fraction', 1.0)
        
        # Use deterministic sampling for reproducibility
        np.random.seed(hash(smi_file) % (2**32))
        
        if self.verbose:
            print(f"  Using ultra-fast seeking method for fraction {database_fraction:.4f}")
        
        # Step 1: Build a map of line positions by sampling the file
        line_positions = []
        
        with open(smi_file, 'rb') as f:
            # Sample every 1000th line to build position map efficiently
            pos = 0
            line_count = 0
            sample_interval = 1000
            
            while True:
                if line_count % sample_interval == 0:
                    line_positions.append(pos)
                
                line = f.readline()
                if not line:
                    break
                    
                pos = f.tell()
                line_count += 1
                
                # For very large files, we don't need to map every position
                # Stop mapping after we have enough reference points
                if len(line_positions) > 10000:  # 10M lines mapped
                    # Estimate remaining file
                    f.seek(0, 2)  # Seek to end
                    total_file_size = f.tell()
                    avg_line_size = pos / line_count if line_count > 0 else 50
                    estimated_total_lines = int(total_file_size / avg_line_size)
                    break
            else:
                # File was fully read
                estimated_total_lines = line_count
        
        if self.verbose:
            print(f"  Estimated {estimated_total_lines:,} total lines, mapped {len(line_positions):,} positions")
        
        # Step 2: Calculate which lines to sample
        target_samples = int(estimated_total_lines * database_fraction)
        target_samples = min(target_samples, self.batch_size * 100)  # Reasonable upper limit
        
        if target_samples == 0:
            if self.verbose:
                print(f"  No samples needed for fraction {database_fraction}")
            return
        
        if self.verbose:
            print(f"  Generating {target_samples:,} random sample positions...")
        
        # Generate random line numbers to sample efficiently
        # For sparse sampling on huge files, use geometric distribution for speed
        if estimated_total_lines > 1000000 and target_samples < estimated_total_lines // 100:
            if self.verbose:
                print(f"    Using geometric sampling for sparse coverage...")
            
            # Ultra-fast geometric sampling - no duplicates by design
            sample_line_numbers = []
            current_pos = 0
            avg_gap = estimated_total_lines / target_samples
            
            # Generate positions using exponential gaps
            for i in range(target_samples):
                # Add some randomness to the gap while maintaining coverage
                gap_variation = avg_gap * 0.5  # 50% variation around average
                gap = max(1, int(np.random.exponential(avg_gap - gap_variation) + gap_variation))
                current_pos += gap
                
                if current_pos >= estimated_total_lines:
                    break
                    
                sample_line_numbers.append(current_pos)
                
                # Quick progress for very large samples
                if self.verbose and i > 0 and i % 5000 == 0:
                    print(f"      Generated {i:,}/{target_samples:,} positions")
            
            sample_line_numbers = np.array(sample_line_numbers)
        else:
            # For smaller ranges or dense sampling, use numpy choice
            sample_line_numbers = np.sort(np.random.choice(
                estimated_total_lines, 
                size=target_samples, 
                replace=False
            ))
        
        if self.verbose:
            print(f"  Will sample {len(sample_line_numbers):,} lines")
        
        # Step 3: Use seeking to read only the selected lines
        df_batch = []
        selected_indices = []
        lines_processed = 0
        
        # Progress reporting
        total_samples = len(sample_line_numbers)
        next_progress_report = max(1, total_samples // 20)  # Report every 5%
        
        if self.verbose:
            seeking_start = time.time()
            print(f"  Starting to seek and process {total_samples:,} target lines...")
        else:
            print(f"  Processing {total_samples:,} sampled lines...")
        
        with open(smi_file, 'r') as f:
            for i, sample_line_num in enumerate(sample_line_numbers):
                # Progress reporting
                if i > 0 and (i % next_progress_report == 0 or i == total_samples - 1):
                    progress_pct = (i + 1) * 100.0 / total_samples
                    if self.verbose:
                        print(f"    Progress: {progress_pct:.1f}% ({i+1:,}/{total_samples:,} lines sampled)")
                    else:
                        print(f"    Seeking progress: {progress_pct:.0f}%")
                
                # Calculate approximate byte position for this line
                # Calculate approximate byte position for this line
                if sample_line_num < len(line_positions) * sample_interval:
                    # We have a mapped position
                    pos_index = sample_line_num // sample_interval
                    if pos_index < len(line_positions):
                        # Seek to the mapped position
                        f.seek(line_positions[pos_index])
                        
                        # Read forward to the exact line
                        lines_to_skip = sample_line_num % sample_interval
                        for _ in range(lines_to_skip):
                            line = f.readline()
                            if not line:
                                break
                        
                        # Read the target line
                        target_line = f.readline()
                        if target_line:
                            parts = target_line.strip().split()
                            if len(parts) >= 2:
                                df_batch.append(parts[:2])
                                selected_indices.append(min(sample_line_num, len(prop_arr)-1))
                            elif len(parts) == 1:
                                df_batch.append([parts[0], f"mol_{sample_line_num}"])
                                selected_indices.append(min(sample_line_num, len(prop_arr)-1))
                else:
                    # Line is beyond our mapped region, estimate position
                    if line_positions:
                        avg_line_size = line_positions[-1] / (len(line_positions) * sample_interval)
                        estimated_pos = int(sample_line_num * avg_line_size)
                        
                        # Seek to estimated position and find next line boundary
                        f.seek(max(0, estimated_pos - 100))  # Seek a bit before to find line start
                        f.readline()  # Skip partial line
                        
                        # Read the line at this position
                        target_line = f.readline()
                        if target_line:
                            parts = target_line.strip().split()
                            if len(parts) >= 2:
                                df_batch.append(parts[:2])
                                selected_indices.append(min(sample_line_num, len(prop_arr)-1))
                            elif len(parts) == 1:
                                df_batch.append([parts[0], f"mol_{sample_line_num}"])
                                selected_indices.append(min(sample_line_num, len(prop_arr)-1))
                
                # Process batch when it reaches batch_size
                if len(df_batch) >= self.batch_size:
                    if df_batch:
                        if self.verbose:
                            print(f"      Processing batch of {len(df_batch)} molecules...")
                        self._process_batch(df_batch, prop_arr, selected_indices)
                        lines_processed += len(df_batch)
                        if self.verbose:
                            print(f"      Batch processed. Total molecules processed: {lines_processed}")
                    df_batch = []
                    selected_indices = []
        
        # Process any remaining batch
        if df_batch:
            if self.verbose:
                print(f"    Processing final batch of {len(df_batch)} molecules...")
            self._process_batch(df_batch, prop_arr, selected_indices)
            lines_processed += len(df_batch)
        
        if self.verbose:
            print(f"  Successfully sampled {len(sample_line_numbers):,} lines using seeking")
            print(f"  Total valid molecules processed: {lines_processed}")
        else:
            print(f"  Seeking completed: {lines_processed} molecules processed")
    
    def _process_batch(self, df_batch: list, prop_arr: np.ndarray, 
                      selected_indices: list) -> None:
        """
        Process a batch of molecules.
        
        Args:
            df_batch: List of [SMILES, name] pairs
            prop_arr: Property array for the current bundle
            selected_indices: List of indices in prop_arr corresponding to df_batch entries
        """
        if not df_batch or not selected_indices:
            return  # Nothing to process
        
        if self.verbose:
            batch_start = time.time()
        
        # Create DataFrame from batch
        df_batch_df = pd.DataFrame(df_batch, columns=["smi", "name"])
        
        # Get corresponding property slice using selected indices
        arr_batch = prop_arr[selected_indices, :][:, self.db_props.column_mask]
        
        # Prepare arguments for multiprocessing
        args_list = [
            (name, molecule, df_batch_df, arr_batch, self.config,
             self.db_props, self.n_oversampled_decoys)
            for name, molecule in self.actives.items()
        ]
        
        if self.verbose:
            print(f"        Starting parallel processing of {len(df_batch)} molecules for {len(self.actives)} targets...")
        
        # Process in parallel
        with Pool(self.n_proc) as pool:
            results = pool.map(process_molecule_batch, args_list)
        
        # Update actives with results
        self.actives = {name: molecule for name, molecule in results}
        
        if self.verbose:
            batch_time = time.time() - batch_start
            print(f"        Batch processing completed in {batch_time:.2f}s")
            
            # Show current decoy counts
            for name, molecule in self.actives.items():
                print(f"        {name}: {molecule.num_decoys} decoys found")
    
    def _apply_similarity_filtering(self) -> None:
        """
        Apply similarity filtering to all molecules' decoys.
        """
        for name, molecule in self.actives.items():
            if molecule.num_decoys > 0 and 'smi' in molecule.decoys.columns:
                # Get current decoys
                current_decoys = molecule.decoys
                
                # Apply similarity filter
                smi_list = current_decoys['smi'].tolist()
                target_size = min(len(smi_list), self.config['Max_decoys_per_ligand'] * self.oversample)
                
                if len(smi_list) > 1:  # Only filter if we have multiple molecules
                    # Sort by score first to prioritize best molecules
                    if 'score' in current_decoys.columns:
                        sorted_decoys = current_decoys.sort_values('score')
                        sorted_smiles = sorted_decoys['smi'].tolist()
                    else:
                        sorted_smiles = smi_list
                    
                    # Apply similarity filter
                    keep_mask = sim_filter(sorted_smiles, target_size, self.config['Max_tc'])
                    
                    if 'score' in current_decoys.columns:
                        # Apply mask to sorted dataframe
                        filtered_decoys = sorted_decoys.iloc[keep_mask].copy()
                    else:
                        # Apply mask to original dataframe
                        filtered_indices = current_decoys.index[keep_mask]
                        filtered_decoys = current_decoys.loc[filtered_indices].copy()
                    
                    # Update molecule decoys
                    molecule.decoys = filtered_decoys
                    
                    if self.verbose and len(smi_list) != len(filtered_decoys):
                        print(f"    {name}: {len(smi_list)} -> {len(filtered_decoys)} decoys after similarity filtering")
    
    def generate_decoys(self) -> Dict[str, Molecule]:
        """
        Generate decoys by processing all available bundles.
        
        Returns:
            Dict[str, Molecule]: Updated dictionary of molecules with decoys
        """
        bundle_index = 0
        processed_bundles = 0
        total_start_time = time.time() if self.verbose else None
        
        while True:
            if not self.db_props.validate_bundle_files(bundle_index):
                if self.verbose:
                    print(f"Bundle {bundle_index} not found, stopping.")
                else:
                    print(f"Bundle {bundle_index} not found, stopping.")
                break
            
            if self.verbose:
                print(f"Processing bundle {bundle_index}...")
                bundle_start = time.time()
            else:
                print(f"Processing bundle {bundle_index}...")
            
            try:
                self.process_bundle(bundle_index)
                
                # Apply similarity filtering after each bundle if enabled
                if self.apply_similarity_filter:
                    if self.verbose:
                        filter_start = time.time()
                        print(f"  Applying similarity filtering (cutoff: {self.config['Max_tc']})...")
                    
                    self._apply_similarity_filtering()
                    
                    if self.verbose:
                        filter_time = time.time() - filter_start
                        print(f"  Similarity filtering completed ({filter_time:.2f}s)")
                
                processed_bundles += 1
                
                if self.verbose:
                    bundle_time = time.time() - bundle_start
                    elapsed_total = time.time() - total_start_time
                    avg_time_per_bundle = elapsed_total / (processed_bundles)
                    
                    print(f"Completed bundle {bundle_index} ({bundle_time:.2f}s, "
                          f"avg: {avg_time_per_bundle:.2f}s/bundle)")
                else:
                    print(f"Completed bundle {bundle_index}")
                    
            except Exception as e:
                print(f"Error processing bundle {bundle_index}: {e}")
                break
            
            bundle_index += 1
        
        if self.verbose and total_start_time:
            total_time = time.time() - total_start_time
            avg_time = total_time / processed_bundles if processed_bundles > 0 else 0
            print(f"Processed {processed_bundles} bundles (total: {total_time:.2f}s, avg: {avg_time:.2f}s/bundle)")
        else:
            print(f"Processed {processed_bundles} bundles")
        
        return self.actives
    
    def finalize_decoys(self) -> Dict[str, Molecule]:
        """
        Finalize decoy selection by applying similarity filtering to reach the target number.
        
        Returns:
            Dict[str, Molecule]: Final molecules with filtered decoys
        """
        target_decoys = self.config['Max_decoys_per_ligand']
        
        for name, molecule in self.actives.items():
            if molecule.num_decoys > 0 and 'smi' in molecule.decoys.columns:
                # Get current decoys
                current_decoys = molecule.decoys
                smi_list = current_decoys['smi'].tolist()
                
                if len(smi_list) > 1:
                    # Sort by score first to prioritize best molecules
                    if 'score' in current_decoys.columns:
                        sorted_decoys = current_decoys.sort_values('score')
                        sorted_smiles = sorted_decoys['smi'].tolist()
                    else:
                        sorted_smiles = smi_list
                    
                    # Apply similarity filter to reach target size
                    keep_mask = sim_filter(sorted_smiles, target_decoys, self.config['Max_tc'])
                    
                    if 'score' in current_decoys.columns:
                        # Apply mask to sorted dataframe
                        filtered_decoys = sorted_decoys.iloc[keep_mask].copy()
                    else:
                        # Apply mask to original dataframe
                        filtered_indices = current_decoys.index[keep_mask]
                        filtered_decoys = current_decoys.loc[filtered_indices].copy()
                    
                    # Update molecule decoys
                    molecule.decoys = filtered_decoys
                    
                    if self.verbose and len(smi_list) != len(filtered_decoys):
                        print(f"Final filtering {name}: {len(smi_list)} -> {len(filtered_decoys)} decoys")
        
        return self.actives
    
    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of the decoy generation results.
        
        Returns:
            Dict[str, Any]: Summary statistics
        """
        summaries = {name: mol.get_summary() for name, mol in self.actives.items()}
        
        total_decoys = sum(mol.num_decoys for mol in self.actives.values())
        
        return {
            'num_target_molecules': len(self.actives),
            'total_decoys_generated': total_decoys,
            'average_decoys_per_molecule': total_decoys / len(self.actives) if self.actives else 0,
            'config_summary': self.db_props.get_config_summary(),
            'molecule_summaries': summaries
        }