"""
Integration layer for the new parallel bundle architecture.

This module provides a drop-in replacement for the existing DecoysGenerator
that uses the new parallel bundle processing architecture.
"""

import os
import json
import tempfile
import pandas as pd
from typing import Dict, Optional

from parallel_generator import ParallelBundleGenerator
from models import Molecule
from database import DatabaseProperties
from utils import get_charge, predict_charge


class ParallelDecoysGenerator:
    """
    Drop-in replacement for DecoysGenerator using parallel bundle architecture.
    
    This class maintains the same interface as the original DecoysGenerator
    but uses the new parallel processing implementation underneath.
    """
    
    def __init__(self, enamine_std_path: str, prop_bundle_dir: str, smi_bundle_dir: str,
                 prop_prefix: str, smi_prefix: str, config_file: str, ligand_smi_file: str,
                 oversample: int = 10, n_proc: int = 8, batch_size: int = 10000, 
                 verbose: bool = False, apply_similarity_filter: bool = True):
        """
        Initialize the parallel decoy generator.
        
        Args match the original DecoysGenerator for compatibility.
        """
        self.enamine_std_path = enamine_std_path
        self.prop_bundle_dir = prop_bundle_dir
        self.smi_bundle_dir = smi_bundle_dir
        self.prop_prefix = prop_prefix
        self.smi_prefix = smi_prefix
        self.config_file = config_file
        self.ligand_smi_file = ligand_smi_file
        self.oversample = oversample
        self.n_proc = n_proc
        self.batch_size = batch_size
        self.verbose = verbose
        self.apply_similarity_filter = apply_similarity_filter
        
        # Initialize database handler for compatibility
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
        
        # Prepare configuration for parallel processing
        self.bundle_config = {
            'enamine_std_path': enamine_std_path,
            'prop_bundle_dir': prop_bundle_dir,
            'smi_bundle_dir': smi_bundle_dir,
            'prop_prefix': prop_prefix,
            'smi_prefix': smi_prefix,
            'config_file': config_file  # Use the actual config file, not temp
        }
        
        self.sampling_config = {
            'database_fraction': self.config.get('Database_fraction', 1.0),
            'charge_filter': self.config.get('Charge', False),
            'protonate': self.config.get('Protonate', False),
            'max_decoys_per_target': self.n_oversampled_decoys,
            'similarity_cutoff': self.config.get('Max_tc', 0.35)
        }
    
    def _load_ligands(self) -> Dict[str, Molecule]:
        """Load target ligand molecules from file (same as original)."""
        if not os.path.exists(self.ligand_smi_file):
            raise FileNotFoundError(f"Ligand file not found: {self.ligand_smi_file}")
        
        actives = {}
        
        try:
            with open(self.ligand_smi_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split()
                    if len(parts) != 2:
                        raise ValueError(
                            f"Invalid format on line {line_num}. Expected: 'SMILES NAME'"
                        )
                    
                    smiles, name = parts
                    
                    # Calculate charge and properties
                    if self.config.get('Protonate', False):
                        charge = predict_charge(smiles)
                    else:
                        charge = get_charge(smiles)
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
    
    def generate_decoys(self) -> Dict[str, Molecule]:
        """
        Generate decoys using parallel bundle processing.
        
        Returns:
            Dict[str, Molecule]: Updated dictionary of molecules with decoys
        """
        if self.verbose:
            print("Using new parallel bundle architecture")
        
        # Create parallel processor
        parallel_processor = ParallelBundleGenerator(
            bundle_config=self.bundle_config,
            target_molecules=self.actives,
            sampling_config=self.sampling_config,
            n_processes=self.n_proc,
            verbose=self.verbose
        )
        
        # Process bundles in parallel
        updated_actives = parallel_processor.process_bundles_parallel()
        
        # Apply similarity filtering if enabled
        if self.apply_similarity_filter:
            if self.verbose:
                print(f"Applying similarity filtering (cutoff: {self.config['Max_tc']})...")
            self._apply_similarity_filtering(updated_actives)
        
        return updated_actives
    
    def _apply_similarity_filtering(self, actives: Dict[str, Molecule]):
        """Apply similarity filtering to molecules (same as original)."""
        from utils import sim_filter
        
        for name, molecule in actives.items():
            if molecule.num_decoys > 0 and 'smi' in molecule.decoys.columns:
                current_decoys = molecule.decoys
                smi_list = current_decoys['smi'].tolist()
                target_size = min(len(smi_list), self.config['Max_decoys_per_ligand'] * self.oversample)
                
                if len(smi_list) > 1:
                    if 'score' in current_decoys.columns:
                        sorted_decoys = current_decoys.sort_values('score')
                        sorted_smiles = sorted_decoys['smi'].tolist()
                    else:
                        sorted_smiles = smi_list
                    
                    keep_mask = sim_filter(sorted_smiles, target_size, self.config['Max_tc'])
                    
                    if 'score' in current_decoys.columns:
                        filtered_decoys = sorted_decoys.iloc[keep_mask].copy()
                    else:
                        filtered_indices = current_decoys.index[keep_mask]
                        filtered_decoys = current_decoys.loc[filtered_indices].copy()
                    
                    molecule.decoys = filtered_decoys
                    
                    if self.verbose and len(smi_list) != len(filtered_decoys):
                        print(f"    {name}: {len(smi_list)} -> {len(filtered_decoys)} decoys after similarity filtering")
    
    def finalize_decoys(self) -> Dict[str, Molecule]:
        """Finalize decoy selection (same as original)."""
        target_decoys = self.config['Max_decoys_per_ligand']
        
        for name, molecule in self.actives.items():
            if molecule.num_decoys > 0 and 'smi' in molecule.decoys.columns:
                current_decoys = molecule.decoys
                smi_list = current_decoys['smi'].tolist()
                
                if len(smi_list) > 1:
                    if 'score' in current_decoys.columns:
                        sorted_decoys = current_decoys.sort_values('score')
                        sorted_smiles = sorted_decoys['smi'].tolist()
                    else:
                        sorted_smiles = smi_list
                    
                    from utils import sim_filter
                    keep_mask = sim_filter(sorted_smiles, target_decoys, self.config['Max_tc'])
                    
                    if 'score' in current_decoys.columns:
                        filtered_decoys = sorted_decoys.iloc[keep_mask].copy()
                    else:
                        filtered_indices = current_decoys.index[keep_mask]
                        filtered_decoys = current_decoys.loc[filtered_indices].copy()
                    
                    molecule.decoys = filtered_decoys
                    
                    if self.verbose and len(smi_list) != len(filtered_decoys):
                        print(f"Final filtering {name}: {len(smi_list)} -> {len(filtered_decoys)} decoys")
        
        return self.actives
    
    def get_summary(self) -> Dict[str, any]:
        """Get summary statistics (same as original)."""
        summaries = {name: mol.get_summary() for name, mol in self.actives.items()}
        
        total_decoys = sum(mol.num_decoys for mol in self.actives.values())
        
        return {
            'num_target_molecules': len(self.actives),
            'total_decoys_generated': total_decoys,
            'average_decoys_per_molecule': total_decoys / len(self.actives) if self.actives else 0,
            'config_summary': self.db_props.get_config_summary(),
            'molecule_summaries': summaries
        }


# For backwards compatibility, allow importing the parallel version as the main generator
DecoysGenerator = ParallelDecoysGenerator