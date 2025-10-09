#!/usr/bin/env python3
"""
Minimal Decoy Generation System

A clean, simple implementation focused on processing massive SMILES databases
efficiently with 0.01% sampling for molecular decoy generation.

Key principles:
- Single file with clear, linear flow
- Minimal dependencies (RDKit, pandas, numpy)
- Optimized for 2B+ line files with sparse sampling
- Simple data structures and functions
- No unnecessary abstractions
"""

import os
import sys
import json
import time
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Dict, NamedTuple, Optional
from collections import namedtuple

# Simple data structures
Target = namedtuple('Target', ['name', 'smiles', 'properties', 'threshold', 'charge'])
Decoy = namedtuple('Decoy', ['smiles', 'name', 'score', 'properties', 'target_name'])


def calculate_properties(smiles: str) -> np.ndarray:
    """Calculate molecular properties using RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        # Calculate core properties
        props = [
            Descriptors.HeavyAtomCount(mol),
            Descriptors.MolLogP(mol), 
            rdMolDescriptors.CalcNumHBA(mol),
            rdMolDescriptors.CalcNumHBD(mol),
            Chem.GetFormalCharge(mol)
        ]
        
        return np.array(props, dtype=np.float32)
        
    except ImportError:
        print("Error: RDKit not available. Please install: pip install rdkit")
        sys.exit(1)
    except Exception:
        return None


def get_charge_pattern(smiles: str) -> List[bool]:
    """Get charge pattern [has_positive, has_negative] for SMILES string."""
    has_positive = '+' in smiles and not any(nitro in smiles for nitro in ['[N+](=O)[O-]', '[n+]([o-])'])
    has_negative = '-' in smiles
    return [has_positive, has_negative]


def load_targets(smiles_file: str, property_weights: np.ndarray, threshold_factor: float = 0.2) -> List[Target]:
    """Load target molecules from SMILES file."""
    if not os.path.exists(smiles_file):
        raise FileNotFoundError(f"Target file not found: {smiles_file}")
    
    targets = []
    
    with open(smiles_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            parts = line.split()
            if len(parts) != 2:
                print(f"Warning: Invalid format on line {line_num}, skipping")
                continue
            
            smiles, name = parts
            
            # Calculate properties
            properties = calculate_properties(smiles)
            if properties is None:
                print(f"Warning: Could not calculate properties for {name}, skipping")
                continue
            
            # Calculate charge
            charge = get_charge_pattern(smiles)
            
            # Calculate threshold
            threshold = threshold_factor * len(properties)
            
            targets.append(Target(
                name=name,
                smiles=smiles, 
                properties=properties,
                threshold=threshold,
                charge=charge
            ))
    
    print(f"Loaded {len(targets)} target molecules")
    return targets


def generate_sample_indices(total_lines: int, fraction: float, seed_string: str) -> List[int]:
    """Generate indices of lines to sample using optimized geometric sampling."""
    if fraction >= 1.0:
        return list(range(total_lines))
    
    # Deterministic sampling
    np.random.seed(hash(seed_string) % (2**32))
    
    target_samples = int(total_lines * fraction)
    if target_samples == 0:
        return []
    
    if fraction <= 0.01:
        # Ultra-fast geometric sampling for sparse sampling
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


def read_sampled_molecules(smiles_file: str, sample_indices: List[int]) -> List[Tuple[str, str]]:
    """Read molecules at specific line indices."""
    if not sample_indices:
        return []
    
    molecules = []
    sample_set = set(sample_indices)
    
    with open(smiles_file, 'r') as f:
        for line_idx, line in enumerate(f):
            if line_idx in sample_set:
                parts = line.strip().split()
                if len(parts) >= 2:
                    molecules.append((parts[0], parts[1]))
                elif len(parts) == 1:
                    molecules.append((parts[0], f"mol_{line_idx}"))
                
                if len(molecules) >= len(sample_indices):
                    break
    
    return molecules


def find_decoys_for_targets(molecules: List[Tuple[str, str]], targets: List[Target], 
                           property_weights: np.ndarray, max_decoys: int = 500,
                           charge_filter: bool = False) -> Dict[str, List[Decoy]]:
    """Find decoys for all target molecules."""
    results = {target.name: [] for target in targets}
    
    # Calculate properties for all molecules
    valid_molecules = []
    all_properties = []
    
    print(f"Calculating properties for {len(molecules)} molecules...")
    
    for smiles, name in molecules:
        props = calculate_properties(smiles)
        if props is not None:
            valid_molecules.append((smiles, name))
            all_properties.append(props)
    
    if not valid_molecules:
        print("No valid molecules found after property calculation")
        return results
    
    all_properties = np.array(all_properties)
    print(f"Valid molecules with properties: {len(valid_molecules)}")
    
    # Process each target
    for target in targets:
        print(f"Finding decoys for {target.name}...")
        
        # Calculate similarity scores
        scores = (
            (np.abs(all_properties - target.properties) / property_weights)
            .sum(axis=1)
        )
        
        # Filter by threshold
        valid_mask = scores < target.threshold
        if not valid_mask.any():
            print(f"  No decoys found for {target.name} (threshold too strict)")
            continue
        
        # Apply charge filtering if enabled
        if charge_filter:
            charge_mask = []
            for i, (smiles, _) in enumerate(valid_molecules):
                if valid_mask[i]:
                    mol_charge = get_charge_pattern(smiles)
                    charge_mask.append(mol_charge == target.charge)
                else:
                    charge_mask.append(False)
            charge_mask = np.array(charge_mask)
            valid_mask = valid_mask & charge_mask
        
        if not valid_mask.any():
            print(f"  No decoys found for {target.name} after charge filtering")
            continue
        
        # Collect valid decoys
        target_decoys = []
        valid_indices = np.where(valid_mask)[0]
        
        for idx in valid_indices:
            smiles, name = valid_molecules[idx]
            score = scores[idx]
            properties = all_properties[idx]
            
            target_decoys.append(Decoy(
                smiles=smiles,
                name=name,
                score=score,
                properties=properties,
                target_name=target.name
            ))
        
        # Sort by score and limit
        target_decoys.sort(key=lambda x: x.score)
        target_decoys = target_decoys[:max_decoys]
        
        results[target.name] = target_decoys
        print(f"  Found {len(target_decoys)} decoys for {target.name}")
    
    return results


def process_bundle(bundle_dir: str, bundle_prefix: str, bundle_id: int, 
                  targets: List[Target], property_weights: np.ndarray,
                  sampling_fraction: float = 0.0001, max_decoys: int = 500,
                  charge_filter: bool = False) -> Dict[str, List[Decoy]]:
    """Process a single SMILES bundle."""
    
    # Construct file paths
    smiles_file = os.path.join(bundle_dir, f"{bundle_prefix}{bundle_id}.smi")
    props_file = os.path.join(bundle_dir.replace('smiles', 'props'), f"props_{bundle_id}.npy")
    
    if not os.path.exists(smiles_file):
        print(f"Bundle file not found: {smiles_file}")
        return {target.name: [] for target in targets}
    
    print(f"Processing bundle {bundle_id}: {smiles_file}")
    
    # Get exact line count from property file if available
    if os.path.exists(props_file):
        try:
            props_array = np.load(props_file, mmap_mode='r')
            total_lines = props_array.shape[0]
            print(f"  Exact line count from properties: {total_lines:,}")
        except Exception:
            print(f"  Could not load properties file, estimating line count...")
            total_lines = sum(1 for _ in open(smiles_file))
    else:
        print(f"  Properties file not found, counting lines...")
        total_lines = sum(1 for _ in open(smiles_file))
    
    # Generate sample indices
    sample_indices = generate_sample_indices(total_lines, sampling_fraction, smiles_file)
    print(f"  Sampling {len(sample_indices):,} lines ({sampling_fraction:.4%})")
    
    if not sample_indices:
        return {target.name: [] for target in targets}
    
    # Read sampled molecules
    molecules = read_sampled_molecules(smiles_file, sample_indices)
    print(f"  Read {len(molecules)} sampled molecules")
    
    if not molecules:
        return {target.name: [] for target in targets}
    
    # Find decoys
    return find_decoys_for_targets(molecules, targets, property_weights, max_decoys, charge_filter)


def save_results(results: Dict[str, List[Decoy]], output_dir: str = "results"):
    """Save decoy results to CSV files."""
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    total_decoys = 0
    
    for target_name, decoys in results.items():
        if not decoys:
            print(f"No decoys to save for {target_name}")
            continue
        
        # Convert to DataFrame
        data = []
        for decoy in decoys:
            row = {
                'smi': decoy.smiles,
                'name': decoy.name,
                'score': decoy.score,
                'target': decoy.target_name
            }
            # Add property columns
            prop_names = ['NumHeavyAtoms', 'logP', 'NumHAcceptors', 'NumHDonors', 'Charge']
            for i, prop_name in enumerate(prop_names):
                if i < len(decoy.properties):
                    row[prop_name] = decoy.properties[i]
            
            data.append(row)
        
        df = pd.DataFrame(data)
        output_file = output_path / f"{target_name}_decoys.csv"
        df.to_csv(output_file, index=False)
        
        total_decoys += len(decoys)
        print(f"Saved {len(decoys)} decoys for {target_name} to {output_file}")
    
    print(f"Total decoys saved: {total_decoys}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Minimal molecular decoy generation for massive SMILES databases",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("--targets", required=True, help="Target molecules file (SMILES NAME format)")
    parser.add_argument("--bundle-dir", required=True, help="Directory containing SMILES bundles")
    parser.add_argument("--bundle-prefix", default="smiles_", help="Prefix for bundle files")
    
    # Optional arguments  
    parser.add_argument("--output-dir", default="results", help="Output directory")
    parser.add_argument("--sampling-fraction", type=float, default=0.0001, help="Fraction of database to sample")
    parser.add_argument("--max-decoys", type=int, default=500, help="Maximum decoys per target")
    parser.add_argument("--threshold", type=float, default=0.2, help="Similarity threshold factor")
    parser.add_argument("--charge-filter", action="store_true", help="Enable charge filtering")
    parser.add_argument("--bundles", help="Specific bundle IDs to process (e.g., '0,1,2' or '0-5')")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.targets):
        print(f"Error: Target file not found: {args.targets}")
        sys.exit(1)
    
    if not os.path.exists(args.bundle_dir):
        print(f"Error: Bundle directory not found: {args.bundle_dir}")
        sys.exit(1)
    
    print("=" * 60)
    print("Minimal Decoy Generation System")
    print("=" * 60)
    
    # Standard property weights (from Enamine or empirical data)
    property_weights = np.array([10.0, 2.0, 3.0, 3.0, 1.0], dtype=np.float32)
    
    print(f"Configuration:")
    print(f"  Sampling fraction: {args.sampling_fraction:.4%}")
    print(f"  Max decoys per target: {args.max_decoys}")
    print(f"  Similarity threshold factor: {args.threshold}")
    print(f"  Charge filtering: {'enabled' if args.charge_filter else 'disabled'}")
    
    # Load targets
    print(f"\nLoading targets from {args.targets}...")
    targets = load_targets(args.targets, property_weights, args.threshold)
    
    if not targets:
        print("No valid targets found")
        sys.exit(1)
    
    # Determine bundles to process
    if args.bundles:
        if '-' in args.bundles:
            start, end = map(int, args.bundles.split('-'))
            bundle_ids = list(range(start, end + 1))
        else:
            bundle_ids = [int(x.strip()) for x in args.bundles.split(',')]
    else:
        # Auto-discover bundles
        bundle_ids = []
        bundle_id = 0
        while True:
            bundle_file = os.path.join(args.bundle_dir, f"{args.bundle_prefix}{bundle_id}.smi")
            if os.path.exists(bundle_file):
                bundle_ids.append(bundle_id)
                bundle_id += 1
            else:
                break
    
    print(f"Found {len(bundle_ids)} bundles to process: {bundle_ids}")
    
    # Process bundles
    all_results = {target.name: [] for target in targets}
    
    start_time = time.time()
    
    for bundle_id in bundle_ids:
        bundle_results = process_bundle(
            args.bundle_dir, args.bundle_prefix, bundle_id,
            targets, property_weights, args.sampling_fraction,
            args.max_decoys, args.charge_filter
        )
        
        # Merge results
        for target_name, decoys in bundle_results.items():
            all_results[target_name].extend(decoys)
    
    # Sort and limit final results
    for target_name in all_results:
        decoys = all_results[target_name]
        decoys.sort(key=lambda x: x.score)
        all_results[target_name] = decoys[:args.max_decoys]
    
    processing_time = time.time() - start_time
    
    # Save results
    print(f"\nSaving results to {args.output_dir}...")
    save_results(all_results, args.output_dir)
    
    # Summary
    total_decoys = sum(len(decoys) for decoys in all_results.values())
    print(f"\nSummary:")
    print(f"  Bundles processed: {len(bundle_ids)}")
    print(f"  Total processing time: {processing_time:.2f}s")
    print(f"  Total decoys found: {total_decoys}")
    print(f"  Average decoys per target: {total_decoys / len(targets):.1f}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()