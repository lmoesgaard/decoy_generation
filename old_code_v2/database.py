"""
Database handler for molecular property databases.

This module handles all database-related operations including property
calculation, data loading, and database access for the decoy generation
process.
"""

import os
import json
import numpy as np
from rdkit import Chem
from numpy.lib.format import open_memmap
from typing import Dict, Any, Tuple

from utils import MOLECULAR_PROPERTIES


class DatabaseProperties:
    """
    Handles molecular property database operations and calculations.
    
    This class manages the molecular property database, including loading
    standardization data, calculating properties for molecules, and providing
    access to database bundles.
    """
    
    def __init__(self, enamine_std_path: str, prop_bundle_dir: str, smi_bundle_dir: str, 
                 prop_prefix: str, smi_prefix: str, config_file: str):
        """
        Initialize the database properties handler.
        
        Args:
            enamine_std_path: Path to the Enamine standard deviations file
            prop_bundle_dir: Directory containing property bundle files
            smi_bundle_dir: Directory containing SMILES bundle files
            prop_prefix: Prefix for property bundle filenames
            smi_prefix: Prefix for SMILES bundle filenames
            config_file: Path to the configuration file
        """
        self.enamine_std_path = enamine_std_path
        self.prop_bundle_dir = prop_bundle_dir
        self.smi_bundle_dir = smi_bundle_dir
        self.prop_prefix = prop_prefix
        self.smi_prefix = smi_prefix
        
        # Load configuration
        self.config = self._load_config(config_file)
        
        # Load standardization data
        self.enamine_std = self._load_enamine_std()
        
        # Create column mask based on configuration
        self.column_mask = self._create_column_mask()
        
        # Apply mask to standardization data
        self.queried_enamine_std = self.enamine_std[self.column_mask]
    
    def _load_config(self, config_file: str) -> Dict[str, Any]:
        """Load configuration from JSON file."""
        try:
            with open(config_file, "r") as f:
                config = json.load(f)
            return config
        except FileNotFoundError:
            raise FileNotFoundError(f"Configuration file not found: {config_file}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in configuration file: {e}")
    
    def _load_enamine_std(self) -> np.ndarray:
        """Load Enamine standard deviations from file."""
        try:
            return np.load(self.enamine_std_path)
        except FileNotFoundError:
            raise FileNotFoundError(f"Enamine std file not found: {self.enamine_std_path}")
    
    def _create_column_mask(self) -> np.ndarray:
        """
        Create a boolean mask for properties to be used based on configuration.
        
        Returns:
            np.ndarray: Boolean mask indicating which properties to use
        """
        mask = []
        for name, _ in MOLECULAR_PROPERTIES:
            if name in self.config:
                # Handle both nested dict format and simple boolean format
                if isinstance(self.config[name], dict):
                    mask.append(self.config[name].get('active', False))
                else:
                    mask.append(bool(self.config[name]))
            else:
                # Default to False if property not in config
                mask.append(False)
        
        return np.array(mask, dtype=bool)
    
    def calculate_all_properties(self, smiles: str) -> np.ndarray:
        """
        Calculate all molecular properties for a given SMILES string.
        
        Args:
            smiles: SMILES string representation of the molecule
            
        Returns:
            np.ndarray: Array of calculated property values
            
        Raises:
            ValueError: If the SMILES string is invalid
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        properties = []
        for _, calc_func in MOLECULAR_PROPERTIES:
            try:
                value = calc_func(mol)
                properties.append(value)
            except Exception as e:
                raise ValueError(f"Failed to calculate property with {calc_func.__name__}: {e}")
        
        return np.array(properties)
    
    def calculate_queried_properties(self, smiles: str) -> np.ndarray:
        """
        Calculate only the properties specified in the configuration.
        
        Args:
            smiles: SMILES string representation of the molecule
            
        Returns:
            np.ndarray: Array of calculated property values for queried properties only
        """
        all_properties = self.calculate_all_properties(smiles)
        return all_properties[self.column_mask]
    
    def get_bundle_filenames(self, bundle_index: int) -> Tuple[str, str]:
        """
        Get the filenames for SMILES and property bundles at a given index.
        
        Args:
            bundle_index: Index of the bundle
            
        Returns:
            Tuple[str, str]: (smiles_file, properties_file) paths
        """
        smi_file = os.path.join(self.smi_bundle_dir, f"{self.smi_prefix}{bundle_index}.smi")
        prop_file = os.path.join(self.prop_bundle_dir, f"{self.prop_prefix}{bundle_index}.npy")
        
        return smi_file, prop_file
    
    def validate_bundle_files(self, bundle_index: int) -> bool:
        """
        Check if bundle files exist for a given index.
        
        Args:
            bundle_index: Index of the bundle to validate
            
        Returns:
            bool: True if both files exist, False otherwise
        """
        smi_file, prop_file = self.get_bundle_filenames(bundle_index)
        return os.path.exists(smi_file) and os.path.exists(prop_file)
    
    def load_property_bundle(self, bundle_index: int) -> np.ndarray:
        """
        Load a property bundle as a memory-mapped array.
        
        Args:
            bundle_index: Index of the bundle to load
            
        Returns:
            np.ndarray: Memory-mapped array of properties
            
        Raises:
            FileNotFoundError: If the bundle file doesn't exist
        """
        _, prop_file = self.get_bundle_filenames(bundle_index)
        
        if not os.path.exists(prop_file):
            raise FileNotFoundError(f"Property bundle not found: {prop_file}")
        
        return open_memmap(prop_file, "r")
    
    def get_active_property_names(self) -> list:
        """
        Get the names of properties that are active in the current configuration.
        
        Returns:
            list: Names of active properties
        """
        return [name for (name, _), active in zip(MOLECULAR_PROPERTIES, self.column_mask) if active]
    
    def get_num_active_properties(self) -> int:
        """
        Get the number of active properties in the current configuration.
        
        Returns:
            int: Number of active properties
        """
        return int(self.column_mask.sum())
    
    def get_config_summary(self) -> Dict[str, Any]:
        """
        Get a summary of the current configuration.
        
        Returns:
            dict: Summary of configuration settings
        """
        return {
            'total_properties': len(MOLECULAR_PROPERTIES),
            'active_properties': self.get_num_active_properties(),
            'active_property_names': self.get_active_property_names(),
            'enamine_std_shape': self.enamine_std.shape,
            'queried_std_shape': self.queried_enamine_std.shape,
            'bundle_directories': {
                'properties': self.prop_bundle_dir,
                'smiles': self.smi_bundle_dir
            }
        }