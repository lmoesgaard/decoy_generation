"""
Parallel Bundle Generator - Clean Architecture

This module implements parallel bundle processing from the ground up,
designed specifically for processing multiple massive SMILES bundles simultaneously.

Key principles:
1. Shared memory for target molecules and configuration
2. Simple, pickleable worker functions  
3. Streaming results aggregation
4. Minimal object copying
5. Clear separation of concerns
"""

import os
import pandas as pd
import numpy as np
import time
from multiprocessing import Pool, Manager, Value, Array
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
import json

from models import Molecule
from database import DatabaseProperties
from utils import get_charge, predict_charge, sim_filter


class BundleResult:
    """Simple container for bundle processing results."""
    
    def __init__(self, bundle_id: int):
        self.bundle_id = bundle_id
        self.decoys = {}  # target_name -> list of (smiles, name, score, properties)
        self.stats = {
            'lines_processed': 0,
            'molecules_found': 0,
            'processing_time': 0.0
        }
    
    def add_decoy(self, target_name: str, smiles: str, name: str, score: float, properties: np.ndarray):
        """Add a decoy to the results."""
        if target_name not in self.decoys:
            self.decoys[target_name] = []
        self.decoys[target_name].append({
            'smiles': smiles,
            'name': name, 
            'score': score,
            'properties': properties
        })


def process_bundle_worker(args: Tuple) -> BundleResult:
    """
    Worker function to process a single bundle in parallel.
    
    This function is designed to be simple and pickleable.
    All complex objects are reconstructed within the worker.
    
    Args:
        args: Tuple containing:
            - bundle_id: Bundle index
            - bundle_config: Simple dict with paths and settings
            - target_molecules: List of simple target molecule dicts
            - sampling_config: Dict with sampling parameters
    
    Returns:
        BundleResult: Container with found decoys and statistics
    """
    bundle_id, bundle_config, target_molecules, sampling_config = args
    
    start_time = time.time()
    result = BundleResult(bundle_id)
    
    try:
        # Reconstruct necessary objects within worker
        db_props = DatabaseProperties(
            bundle_config['enamine_std_path'],
            bundle_config['prop_bundle_dir'],
            bundle_config['smi_bundle_dir'],
            bundle_config['prop_prefix'],
            bundle_config['smi_prefix'],
            bundle_config['config_file']
        )
        
        # Get bundle files
        smi_file, prop_file = db_props.get_bundle_filenames(bundle_id)
        if not os.path.exists(smi_file) or not os.path.exists(prop_file):
            return result
        
        # Load property array
        prop_arr = db_props.load_property_bundle(bundle_id)
        total_lines = prop_arr.shape[0]
        
        # Sample lines to process
        sample_indices = _generate_sample_indices(
            total_lines, 
            sampling_config['database_fraction'],
            smi_file  # for deterministic seeding
        )
        
        if not sample_indices:
            return result
        
        # Process sampled lines
        molecules_data = _read_sampled_molecules(smi_file, sample_indices)
        result.stats['lines_processed'] = len(molecules_data)
        
        if not molecules_data:
            return result
        
        # Extract properties for sampled molecules
        sampled_properties = prop_arr[sample_indices, :][:, db_props.column_mask]
        
        # Process each target molecule
        for target in target_molecules:
            target_decoys = _find_decoys_for_target(
                target, molecules_data, sampled_properties, 
                db_props, sampling_config
            )
            
            if target_decoys:
                result.decoys[target['name']] = target_decoys
                result.stats['molecules_found'] += len(target_decoys)
        
        result.stats['processing_time'] = time.time() - start_time
        
    except Exception as e:
        print(f"Worker error in bundle {bundle_id}: {e}")
        result.stats['processing_time'] = time.time() - start_time
    
    return result


def _generate_sample_indices(total_lines: int, database_fraction: float, smi_file: str) -> List[int]:
    """Generate indices of lines to sample from the SMILES file."""
    if database_fraction >= 1.0:
        return list(range(total_lines))
    
    # Use deterministic sampling
    np.random.seed(hash(smi_file) % (2**32))
    
    target_samples = int(total_lines * database_fraction)
    if target_samples == 0:
        return []
    
    if database_fraction <= 0.01:
        # Geometric sampling for very sparse sampling
        sample_indices = []
        current_pos = 0
        avg_gap = total_lines / target_samples
        
        for _ in range(target_samples):
            gap_variation = avg_gap * 0.5
            gap = max(1, int(np.random.exponential(avg_gap - gap_variation) + gap_variation))
            current_pos += gap
            
            if current_pos >= total_lines:
                break
                
            sample_indices.append(current_pos)
        
        return sample_indices
    else:
        # Random sampling for larger fractions
        return sorted(np.random.choice(total_lines, size=target_samples, replace=False))


def _read_sampled_molecules(smi_file: str, sample_indices: List[int]) -> List[Tuple[str, str]]:
    """Read molecules at specific line indices from SMILES file."""
    if not sample_indices:
        return []
    
    molecules = []
    sample_set = set(sample_indices)
    
    with open(smi_file, 'r') as f:
        for line_idx, line in enumerate(f):
            if line_idx in sample_set:
                parts = line.strip().split()
                if len(parts) >= 2:
                    molecules.append((parts[0], parts[1]))  # smiles, name
                elif len(parts) == 1:
                    molecules.append((parts[0], f"mol_{line_idx}"))
                
                # Early exit if we've found all samples
                if len(molecules) >= len(sample_indices):
                    break
    
    return molecules


def _find_decoys_for_target(target: dict, molecules_data: List[Tuple[str, str]], 
                           properties: np.ndarray, db_props, config: dict) -> List[dict]:
    """Find decoys for a specific target molecule."""
    target_properties = np.array(target['properties'])
    threshold = target['threshold']
    
    # Calculate similarity scores
    scores = (
        (np.absolute(properties - target_properties) / db_props.queried_enamine_std)
        .sum(axis=1)
    )
    
    # Filter by threshold
    valid_mask = scores < threshold
    if not valid_mask.any():
        return []
    
    # Extract valid molecules and their data
    valid_molecules = [molecules_data[i] for i in range(len(molecules_data)) if valid_mask[i]]
    valid_scores = scores[valid_mask]
    valid_properties = properties[valid_mask]
    
    # Apply charge filtering if enabled
    if config.get('charge_filter', False):
        charge_filtered = []
        for i, (smiles, name) in enumerate(valid_molecules):
            if config.get('protonate', False):
                mol_charge = predict_charge(smiles)
            else:
                mol_charge = get_charge(smiles)
            
            if mol_charge == target['charge']:
                charge_filtered.append({
                    'smiles': smiles,
                    'name': name,
                    'score': valid_scores[i],
                    'properties': valid_properties[i]
                })
        return charge_filtered
    else:
        # Return all valid molecules
        return [
            {
                'smiles': smiles,
                'name': name, 
                'score': valid_scores[i],
                'properties': valid_properties[i]
            }
            for i, (smiles, name) in enumerate(valid_molecules)
        ]


class ParallelBundleGenerator:
    """
    Clean, parallel bundle processor designed for massive SMILES databases.
    
    This class processes multiple bundles simultaneously using a simple,
    shared-nothing architecture that avoids complex object serialization.
    """
    
    def __init__(self, bundle_config: dict, target_molecules: Dict[str, Molecule], 
                 sampling_config: dict, n_processes: int = 8, verbose: bool = False):
        """
        Initialize the parallel generator.
        
        Args:
            bundle_config: Dict with database paths and settings
            target_molecules: Dict of target molecules
            sampling_config: Dict with sampling parameters
            n_processes: Number of parallel processes
            verbose: Enable verbose output
        """
        self.bundle_config = bundle_config
        self.target_molecules = target_molecules
        self.sampling_config = sampling_config
        self.n_processes = n_processes
        self.verbose = verbose
        
        # Convert target molecules to simple dicts for serialization
        self.target_data = []
        for name, molecule in target_molecules.items():
            self.target_data.append({
                'name': name,
                'smiles': molecule.smiles,
                'charge': molecule.charge,
                'properties': molecule.properties,
                'threshold': molecule.threshold
            })
    
    def find_available_bundles(self) -> List[int]:
        """Find all available bundle files."""
        available_bundles = []
        bundle_index = 0
        
        # Create temporary DatabaseProperties to check files
        db_props = DatabaseProperties(
            self.bundle_config['enamine_std_path'],
            self.bundle_config['prop_bundle_dir'],
            self.bundle_config['smi_bundle_dir'],
            self.bundle_config['prop_prefix'],
            self.bundle_config['smi_prefix'],
            self.bundle_config['config_file']
        )
        
        while True:
            if db_props.validate_bundle_files(bundle_index):
                available_bundles.append(bundle_index)
                bundle_index += 1
            else:
                break
        
        return available_bundles
    
    def process_bundles_parallel(self) -> Dict[str, Molecule]:
        """
        Process all available bundles in parallel.
        
        Returns:
            Dict[str, Molecule]: Updated target molecules with decoys
        """
        available_bundles = self.find_available_bundles()
        
        if not available_bundles:
            print("No bundles found!")
            return self.target_molecules
        
        print(f"Found {len(available_bundles)} bundles to process")
        
        if len(available_bundles) == 1:
            print("Single bundle detected, using sequential processing")
            return self._process_single_bundle(available_bundles[0])
        
        # Prepare arguments for parallel processing
        args_list = []
        for bundle_id in available_bundles:
            args_list.append((
                bundle_id,
                self.bundle_config,
                self.target_data,
                self.sampling_config
            ))
        
        # Process bundles in parallel
        if self.verbose:
            print(f"Processing {len(available_bundles)} bundles with {self.n_processes} processes...")
            start_time = time.time()
        
        with Pool(self.n_processes) as pool:
            results = pool.map(process_bundle_worker, args_list)
        
        # Aggregate results
        total_decoys_found = 0
        bundle_stats = []
        
        for result in results:
            bundle_stats.append({
                'bundle_id': result.bundle_id,
                'processing_time': result.stats['processing_time'],
                'lines_processed': result.stats['lines_processed'],
                'molecules_found': result.stats['molecules_found']
            })
            
            # Merge decoys into target molecules
            for target_name, decoys in result.decoys.items():
                if target_name in self.target_molecules:
                    self._merge_decoys(target_name, decoys)
                    total_decoys_found += len(decoys)
        
        # Print statistics
        if self.verbose:
            total_time = time.time() - start_time
            print(f"Parallel processing completed in {total_time:.2f}s")
            print(f"Total decoys found: {total_decoys_found}")
            
            for stats in bundle_stats:
                print(f"  Bundle {stats['bundle_id']}: {stats['molecules_found']} decoys "
                      f"({stats['processing_time']:.2f}s, {stats['lines_processed']} lines)")
        else:
            print(f"Parallel processing completed: {total_decoys_found} total decoys found")
        
        return self.target_molecules
    
    def _merge_decoys(self, target_name: str, new_decoys: List[dict]):
        """Merge new decoys into a target molecule."""
        molecule = self.target_molecules[target_name]
        
        if not new_decoys:
            return
        
        # Convert to DataFrame
        new_df = pd.DataFrame(new_decoys)
        
        # Rename properties columns to match expected format
        property_names = ['NumHeavyAtoms', 'logP', 'NumHAcceptors', 'NumHDonors', 'Charge']  # Adjust as needed
        if 'properties' in new_df.columns:
            # Split properties array into separate columns
            props_df = pd.DataFrame(
                new_df['properties'].tolist(), 
                columns=property_names[:len(new_df['properties'].iloc[0])]
            )
            new_df = pd.concat([new_df.drop('properties', axis=1), props_df], axis=1)
        
        # Merge with existing decoys
        if molecule.decoys.empty:
            molecule.decoys = new_df
        else:
            # Combine and keep best decoys
            combined = pd.concat([molecule.decoys, new_df], ignore_index=True)
            
            # Sort by score and keep best
            max_decoys = self.sampling_config.get('max_decoys_per_target', 500)
            if len(combined) > max_decoys:
                combined = combined.nsmallest(max_decoys, 'score')
                molecule.threshold = combined.score.max()
            
            molecule.decoys = combined
    
    def _process_single_bundle(self, bundle_id: int) -> Dict[str, Molecule]:
        """Process a single bundle (fallback for single bundle datasets)."""
        print(f"Processing single bundle {bundle_id}...")
        
        args = (bundle_id, self.bundle_config, self.target_data, self.sampling_config)
        result = process_bundle_worker(args)
        
        # Merge results
        total_decoys = 0
        for target_name, decoys in result.decoys.items():
            if target_name in self.target_molecules:
                self._merge_decoys(target_name, decoys)
                total_decoys += len(decoys)
        
        print(f"Single bundle processing completed: {total_decoys} decoys found")
        return self.target_molecules