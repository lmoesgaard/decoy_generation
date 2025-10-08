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
        Process a single database bundle.
        
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
        with open(smi_file, 'r') as f:
            df_batch = []
            selected_indices = []  # Track which indices were selected
            total_processed = 0
            
            for line_count, line in enumerate(f):
                # Check if we should include this molecule based on database fraction
                if np.random.random() > self.config.get('Database_fraction', 1.0):
                    total_processed += 1
                    continue  # Skip this line based on sampling fraction
                
                # Parse SMILES line
                parts = line.strip().split()
                if len(parts) >= 2:
                    df_batch.append(parts[:2])  # Take only SMILES and name
                    selected_indices.append(total_processed)
                elif len(parts) == 1:
                    df_batch.append([parts[0], f"mol_{line_count}"])  # Generate name if missing
                    selected_indices.append(total_processed)
                else:
                    total_processed += 1
                    continue  # Skip invalid lines
                
                total_processed += 1
                
                # Process batch when it reaches batch_size or end of file
                if (len(df_batch) >= self.batch_size or 
                    total_processed >= prop_arr.shape[0]):
                    
                    if df_batch:  # Only process if we have data
                        self._process_batch(df_batch, prop_arr, selected_indices)
                    
                    df_batch = []
                    selected_indices = []
            
            # Process any remaining batch
            if df_batch:
                self._process_batch(df_batch, prop_arr, selected_indices)
    
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
        
        # Process in parallel
        with Pool(self.n_proc) as pool:
            results = pool.map(process_molecule_batch, args_list)
        
        # Update actives with results
        self.actives = {name: molecule for name, molecule in results}
    
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