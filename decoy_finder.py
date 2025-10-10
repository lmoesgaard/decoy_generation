#!/usr/bin/env python3
"""
Simple Decoy Finder - Clean Implementation with Sampling

Find property-matched molecular decoys from large enumerated libraries.
Processes 2B+ molecule SMILES bundles with precomputed properties using
multiprocessing, adaptive thresholding, and distributed sampling.
"""

import os
import sys
import json
import argparse
import multiprocessing as mp
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors, Crippen, Lipinski
from openbabel import openbabel as ob


def get_charge(mol: ob.OBMol) -> List[int]:
    has_positive = 0
    has_negative = 0

    for atom in ob.OBMolAtomIter(mol):
        formal_charge = atom.GetFormalCharge()
        atomic_num = atom.GetAtomicNum()

        # Skip nitro group charges (N+ with O- neighbors)
        if formal_charge > 0 and atomic_num == 7:
            is_nitro = any(
                n.GetAtomicNum() == 8 and n.GetFormalCharge() < 0
                for n in ob.OBAtomAtomIter(atom)
            )
            if not is_nitro:
                has_positive = 1
        elif formal_charge > 0:
            has_positive = 1
        elif formal_charge < 0 and atomic_num == 8:
            is_nitro = any(
                n.GetAtomicNum() == 7 and n.GetFormalCharge() > 0
                for n in ob.OBAtomAtomIter(atom)
            )
            if not is_nitro:
                has_negative = 1
        elif formal_charge < 0:
            has_negative = 1

    return [has_positive, has_negative]


def generate_sampling_mask(bundle_size: int, sample_fraction: float, bundle_id: int) -> np.ndarray:
    """Generate a distributed boolean sampling mask for a bundle."""
    if sample_fraction >= 1.0:
        return np.ones(bundle_size, dtype=bool)
    
    # Use bundle_id as seed for reproducible sampling
    np.random.seed(bundle_id * 42 + 123)
    
    # Generate random mask distributed across the bundle
    return np.random.random(bundle_size) < sample_fraction


class TargetMolecule:
    """Simple container for target molecule and its properties."""
    
    def __init__(self, smiles: str, name: str):
        self.smiles = smiles
        self.name = name
        self.mol = Chem.MolFromSmiles(smiles)
        
        if self.mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Calculate all properties as specified
        self.properties = np.array([
            rdMolDescriptors.CalcNumHeavyAtoms(self.mol),
            Descriptors.MolWt(self.mol), 
            Crippen.MolLogP(self.mol),
            Lipinski.NumHAcceptors(self.mol), 
            Lipinski.NumHDonors(self.mol),
            rdMolDescriptors.CalcTPSA(self.mol), 
            Descriptors.NumValenceElectrons(self.mol), 
            Lipinski.NumRotatableBonds(self.mol)
        ], dtype=np.float32)
        
        # Calculate charge state using OpenBabel
        self.charge = self._calculate_charge()
    
    def _calculate_charge(self) -> List[int]:
        """Calculate charge state [#positive, #negative] using OpenBabel."""
        try:
            conv = ob.OBConversion()
            conv.SetInAndOutFormats("smi", "smi")
            
            mol = ob.OBMol()
            if not conv.ReadString(mol, self.smiles):
                return [0, 0]
            
            return get_charge(mol)

        except Exception:
            return [0, 0]


def load_targets(smiles_file: str) -> List[TargetMolecule]:
    """Load target molecules from SMILES file."""
    targets = []
    
    with open(smiles_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) != 2:
                print(f"Warning: Skipping invalid line {line_num}: {line}")
                continue
            
            smiles, name = parts
            try:
                target = TargetMolecule(smiles, name)
                targets.append(target)
                print(f"Loaded target: {name} - {smiles}")
            except ValueError as e:
                print(f"Warning: Skipping invalid target {name}: {e}")
    
    print(f"Loaded {len(targets)} target molecules")
    return targets


def calculate_charge_openbabel(smiles: str) -> List[int]:
    """Calculate charge state for a SMILES using OpenBabel."""
    try:
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("smi", "smi")
        
        mol = ob.OBMol()
        if not conv.ReadString(mol, smiles):
            return [0, 0]

        return get_charge(mol)

    except Exception:
        return [0, 0]


def calculate_molecule_properties(smiles: str) -> Optional[np.ndarray]:
    """Calculate molecular properties for a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        return np.array([
            rdMolDescriptors.CalcNumHeavyAtoms(mol),
            Descriptors.MolWt(mol), 
            Crippen.MolLogP(mol),
            Lipinski.NumHAcceptors(mol), 
            Lipinski.NumHDonors(mol),
            rdMolDescriptors.CalcTPSA(mol), 
            Descriptors.NumValenceElectrons(mol), 
            Lipinski.NumRotatableBonds(mol)
        ], dtype=np.float32)
        
    except Exception:
        return None


def calculate_similarity_score(target_props: np.ndarray, candidate_props: np.ndarray, 
                             stds: np.ndarray, active_mask: np.ndarray) -> float:
    """Calculate similarity score between target and candidate properties."""
    # Only consider active properties
    target_active = target_props[active_mask]
    candidate_active = candidate_props[active_mask]
    stds_active = stds[active_mask]
    
    # Prevent division by zero
    stds_active = np.where(stds_active > 1e-10, stds_active, 1.0)
    
    # Calculate average normalized difference
    differences = np.abs(target_active - candidate_active) / stds_active
    return np.mean(differences)


def process_bundle_worker(args: Tuple) -> str:
    """
    Worker function to process a single bundle for all targets with sampling.
    
    Args:
        args: (bundle_id, smi_prefix, prop_prefix, smi_dir, prop_dir, 
               targets, config, stds, output_dir, batch_size, sample_fraction)
    
    Returns:
        str: Status message
    """
    (bundle_id, smi_prefix, prop_prefix, smi_dir, prop_dir, 
     targets, config, stds, output_dir, batch_size, sample_fraction) = args
    
    # File paths
    smi_file = Path(smi_dir) / f"{smi_prefix}{bundle_id}.smi"
    prop_file = Path(prop_dir) / f"{prop_prefix}{bundle_id}.npy"
    
    if not smi_file.exists() or not prop_file.exists():
        return f"Bundle {bundle_id}: Files not found"
    
    # Load properties with memory mapping
    try:
        prop_array = np.load(prop_file, mmap_mode='r')
    except Exception as e:
        return f"Bundle {bundle_id}: Error loading properties - {e}"
    
    # Generate sampling mask
    total_molecules = prop_array.shape[0]
    sampling_mask = generate_sampling_mask(total_molecules, sample_fraction, bundle_id)
    actual_samples = np.sum(sampling_mask)
    
    # Create active property mask
    property_names = ["NumHeavyAtoms", "MWt", "logP", "NumHAcceptors", 
                     "NumHDonors", "CalcTPSA", "NumValenceElectrons", "NumRotatableBonds"]
    
    active_mask = np.array([config.get(name, False) for name in property_names])
    use_charges = config.get("Charge", False)
    
    # Initialize results for each target
    target_results = {}
    for target in targets:
        target_results[target.name] = {
            'best_decoys': [],
            'current_threshold': config.get("Threshold", 0.2),
            'max_decoys': config.get("Max_decoys_per_ligand", 50),
            'oversample': config.get("oversample", 10)
        }
    
    molecules_processed = 0
    
    # Process SMILES file in batches with sampling
    with open(smi_file, 'r') as f:
        batch_smiles = []
        batch_names = []
        batch_indices = []
        
        for line_idx, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            
            # Check if this line should be sampled
            if line_idx < len(sampling_mask) and not sampling_mask[line_idx]:
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                smiles, name = parts[0], parts[1]
                batch_smiles.append(smiles)
                batch_names.append(name)
                batch_indices.append(line_idx)
            
            # Process batch when full
            if len(batch_smiles) >= batch_size:
                molecules_processed += _process_sampled_batch(
                    batch_smiles, batch_names, batch_indices,
                    prop_array, targets, target_results, active_mask, 
                    stds, use_charges, config
                )
                batch_smiles = []
                batch_names = []
                batch_indices = []
        
        # Process remaining molecules
        if batch_smiles:
            molecules_processed += _process_sampled_batch(
                batch_smiles, batch_names, batch_indices,
                prop_array, targets, target_results, active_mask, 
                stds, use_charges, config
            )
    
    # Write results to CSV files
    total_decoys_written = 0
    for target in targets:
        result = target_results[target.name]
        if result['best_decoys']:
            csv_file = Path(output_dir) / f"bundle_{bundle_id}_{target.name}_decoys.csv"
            _write_decoys_csv(result['best_decoys'], csv_file, property_names, active_mask)
            total_decoys_written += len(result['best_decoys'])
    
    return f"Bundle {bundle_id}: Sampled {actual_samples}/{total_molecules} molecules ({sample_fraction:.1%}), processed {molecules_processed}, wrote {total_decoys_written} decoys"


def _process_sampled_batch(batch_smiles: List[str], batch_names: List[str], batch_indices: List[int],
                          prop_array: np.ndarray, targets: List[TargetMolecule], 
                          target_results: Dict, active_mask: np.ndarray, stds: np.ndarray,
                          use_charges: bool, config: Dict) -> int:
    """Process a batch of sampled molecules against all targets."""
    processed = 0
    
    for i, (smiles, name, mol_idx) in enumerate(zip(batch_smiles, batch_names, batch_indices)):
        # Skip if index out of bounds
        if mol_idx >= prop_array.shape[0]:
            continue
        
        # Get precomputed properties
        candidate_props = prop_array[mol_idx]
        
        # Calculate charge if needed
        candidate_charge = None
        if use_charges:
            candidate_charge = calculate_charge_openbabel(smiles)
        
        # Check against all targets
        for target in targets:
            result = target_results[target.name]
            
            # Skip if we have enough high-quality decoys
            max_total = result['max_decoys'] * result['oversample']
            if len(result['best_decoys']) >= max_total:
                continue
            
            # Check charge compatibility if required
            if use_charges and candidate_charge != target.charge:
                continue
            
            # Calculate similarity score
            score = calculate_similarity_score(
                target.properties, candidate_props, stds, active_mask
            )
            
            # Check if this decoy meets the current threshold
            if score <= result['current_threshold']:
                # Add to best decoys
                decoy_entry = {
                    'smiles': smiles,
                    'name': name,
                    'score': score,
                    'properties': candidate_props.copy()
                }
                result['best_decoys'].append(decoy_entry)
                
                # Sort by score and keep only the best
                result['best_decoys'].sort(key=lambda x: x['score'])
                
                # If we have too many, remove worst and tighten threshold
                if len(result['best_decoys']) > max_total:
                    result['best_decoys'] = result['best_decoys'][:max_total]
                    # Tighten threshold to the score of the worst kept decoy
                    result['current_threshold'] = result['best_decoys'][-1]['score']
        
        processed += 1
    
    return processed


def _write_decoys_csv(decoys: List[Dict], csv_file: Path, property_names: List[str], 
                     active_mask: np.ndarray) -> None:
    """Write decoys to CSV file."""
    if not decoys:
        return
    
    # Create DataFrame
    rows = []
    for decoy in decoys:
        row = {
            'smi': decoy['smiles'],
            'name': decoy['name'],
            'score': decoy['score']
        }
        
        # Add active properties
        for i, prop_name in enumerate(property_names):
            if active_mask[i]:
                row[prop_name] = decoy['properties'][i]
        
        rows.append(row)
    
    df = pd.DataFrame(rows)
    df.to_csv(csv_file, index=False)


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Find property-matched molecular decoys from large enumerated libraries"
    )
    
    parser.add_argument("--targets", required=True, help="Target molecules SMILES file (SMILES NAME per line)")
    parser.add_argument("--smi-dir", required=True, help="Directory containing SMILES bundle files")
    parser.add_argument("--prop-dir", required=True, help="Directory containing property bundle files")
    parser.add_argument("--config", required=True, help="Configuration JSON file")
    parser.add_argument("--stds", required=True, help="Standard deviations numpy file")
    parser.add_argument("--smi-prefix", default="smiles_", help="Prefix for SMILES files")
    parser.add_argument("--prop-prefix", default="props_", help="Prefix for property files")
    parser.add_argument("--output-dir", default="decoy_results", help="Output directory")
    parser.add_argument("--batch-size", type=int, default=1000, help="Batch size for processing")
    parser.add_argument("--n-proc", type=int, default=mp.cpu_count(), help="Number of processes")
    parser.add_argument("--sample-fraction", type=float, default=1.0, 
                       help="Fraction of each bundle to sample (0.0-1.0, default: 1.0 = full bundle)")
    
    args = parser.parse_args()
    
    # Validate sample fraction
    if not 0.0 < args.sample_fraction <= 1.0:
        print("Error: sample-fraction must be between 0.0 and 1.0")
        return 1
    
    # Load configuration
    with open(args.config, 'r') as f:
        config = json.load(f)
    
    # Load standard deviations
    stds = np.load(args.stds)
    
    # Load target molecules
    targets = load_targets(args.targets)
    if not targets:
        print("No valid target molecules found!")
        return 1
    
    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Find available bundles
    available_bundles = []
    bundle_id = 0
    while True:
        smi_file = Path(args.smi_dir) / f"{args.smi_prefix}{bundle_id}.smi"
        prop_file = Path(args.prop_dir) / f"{args.prop_prefix}{bundle_id}.npy"
        
        if smi_file.exists() and prop_file.exists():
            available_bundles.append(bundle_id)
            bundle_id += 1
        else:
            break
    
    if not available_bundles:
        print("No bundle files found!")
        return 1
    
    print(f"Found {len(available_bundles)} bundles to process")
    print(f"Processing with {args.n_proc} processes")
    print(f"Sample fraction: {args.sample_fraction:.1%}")
    
    # Prepare arguments for multiprocessing
    worker_args = []
    for bundle_id in available_bundles:
        worker_args.append((
            bundle_id, args.smi_prefix, args.prop_prefix, 
            args.smi_dir, args.prop_dir, targets, config, stds,
            args.output_dir, args.batch_size, args.sample_fraction
        ))
    
    # Process bundles in parallel
    with mp.Pool(args.n_proc) as pool:
        results = pool.map(process_bundle_worker, worker_args)
    
    # Print results
    for result in results:
        print(result)
    
    print("Decoy finding complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())