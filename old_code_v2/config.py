"""
Configuration management for decoy generation.

This module handles all configuration-related functionality, including
loading, validation, and creation of configuration files for the decoy
generation process.
"""

import json
from dataclasses import dataclass, asdict
from typing import Dict, Any, Optional
from pathlib import Path


@dataclass
class DecoyGenerationConfig:
    """
    Configuration for decoy generation parameters.
    
    This class defines all the parameters needed for decoy generation,
    making it easy to modify and validate input parameters.
    """
    
    # Molecular property filters (True = use as filter, False = ignore)
    protonate: bool = False
    molecular_weight: bool = False  
    num_heavy_atoms: bool = True
    log_p: bool = True
    num_h_acceptors: bool = True
    num_h_donors: bool = True
    calc_tpsa: bool = False
    num_valence_electrons: bool = False
    num_rotatable_bonds: bool = False
    charge: bool = True
    
    # Similarity and count constraints
    max_tc: float = 0.35  # Maximum Tanimoto coefficient to ligand
    max_decoys_per_ligand: int = 50  # Maximum number of decoys per ligand
    threshold: float = 0.2  # Threshold value for similarity filtering
    database_fraction: float = 1.0  # Fraction of database to process (1.0 = full database)
    
    def __post_init__(self):
        """Validate parameters after initialization."""
        self.validate()
    
    def validate(self) -> None:
        """Validate that all parameters are within acceptable ranges."""
        if not (0.0 <= self.max_tc <= 1.0):
            raise ValueError(f"max_tc must be between 0.0 and 1.0, got {self.max_tc}")

        if self.max_decoys_per_ligand <= 0:
            raise ValueError(f"max_decoys_per_ligand must be positive, got {self.max_decoys_per_ligand}")
        
        if not (0.0 < self.threshold):
            raise ValueError(f"threshold must be positive, got {self.threshold}")
        
        if not (0.0 < self.database_fraction <= 1.0):
            raise ValueError(f"database_fraction must be between 0.0 and 1.0, got {self.database_fraction}")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the dataclass to a dictionary with original key names."""
        data = asdict(self)
        key_mapping = {
            'protonate': 'Protonate',
            'molecular_weight': 'MWt',
            'num_heavy_atoms': 'NumHeavyAtoms',
            'log_p': 'logP',
            'num_h_acceptors': 'NumHAcceptors',
            'num_h_donors': 'NumHDonors',
            'calc_tpsa': 'CalcTPSA',
            'num_valence_electrons': 'NumValenceElectrons',
            'num_rotatable_bonds': 'NumRotatableBonds',
            'charge': 'Charge',
            'max_tc': 'Max_tc',
            'max_decoys_per_ligand': 'Max_decoys_per_ligand',
            'threshold': 'Threshold',
            'database_fraction': 'Database_fraction'
        }
        return {key_mapping[k]: v for k, v in data.items()}
    
    def save_to_file(self, filename: str = "decoy_generation.json") -> None:
        """
        Save the configuration to a JSON file.
        
        Args:
            filename: Path to save the configuration file
        """
        filepath = Path(filename)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        with open(filepath, "w") as f:
            json.dump(self.to_dict(), f, indent=4)
    
    @classmethod
    def load_from_file(cls, filename: str = "decoy_generation.json") -> "DecoyGenerationConfig":
        """
        Load configuration from a JSON file.
        
        Args:
            filename: Path to the configuration file
            
        Returns:
            DecoyGenerationConfig: Loaded configuration object
            
        Raises:
            FileNotFoundError: If the configuration file doesn't exist
            ValueError: If the configuration file is invalid
        """
        filepath = Path(filename)
        if not filepath.exists():
            raise FileNotFoundError(f"Configuration file not found: {filename}")
        
        try:
            with open(filepath, "r") as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in configuration file: {e}")
        
        # Map from original names back to pythonic names
        reverse_mapping = {
            'Protonate': 'protonate',
            'MWt': 'molecular_weight',
            'NumHeavyAtoms': 'num_heavy_atoms',
            'logP': 'log_p',
            'NumHAcceptors': 'num_h_acceptors',
            'NumHDonors': 'num_h_donors',
            'CalcTPSA': 'calc_tpsa',
            'NumValenceElectrons': 'num_valence_electrons',
            'NumRotatableBonds': 'num_rotatable_bonds',
            'Charge': 'charge',
            'Max_tc': 'max_tc',
            'Max_decoys_per_ligand': 'max_decoys_per_ligand',
            'Threshold': 'threshold',
            'Database_fraction': 'database_fraction'
        }
        
        kwargs = {reverse_mapping[k]: v for k, v in data.items() if k in reverse_mapping}
        
        # Handle missing keys with defaults
        instance = cls()
        for key, value in kwargs.items():
            setattr(instance, key, value)
        
        instance.validate()
        return instance
    
    @classmethod
    def default_config(cls) -> "DecoyGenerationConfig":
        """Create a default configuration."""
        return cls()
    
    @classmethod
    def strict_config(cls) -> "DecoyGenerationConfig":
        """
        Create a strict configuration with more property filters enabled.
        
        Returns:
            DecoyGenerationConfig: Strict configuration
        """
        return cls(
            protonate=True,
            molecular_weight=True,
            num_heavy_atoms=True,
            log_p=True,
            num_h_acceptors=True,
            num_h_donors=True,
            calc_tpsa=True,
            num_valence_electrons=True,
            num_rotatable_bonds=True,
            charge=True,
            max_tc=0.25,
            max_decoys_per_ligand=30,
            threshold=0.15
        )
    
    @classmethod
    def permissive_config(cls) -> "DecoyGenerationConfig":
        """
        Create a permissive configuration with fewer filters.
        
        Returns:
            DecoyGenerationConfig: Permissive configuration
        """
        return cls(
            protonate=False,
            molecular_weight=False,
            num_heavy_atoms=True,
            log_p=True,
            num_h_acceptors=False,
            num_h_donors=False,
            calc_tpsa=False,
            num_valence_electrons=False,
            num_rotatable_bonds=False,
            charge=False,
            max_tc=0.5,
            max_decoys_per_ligand=100,
            threshold=0.3
        )
    
    def get_active_properties(self) -> Dict[str, bool]:
        """
        Get a dictionary of property names and their active status.
        
        Returns:
            Dict[str, bool]: Property names and whether they're active
        """
        return {
            'Protonate': self.protonate,
            'MWt': self.molecular_weight,
            'NumHeavyAtoms': self.num_heavy_atoms,
            'logP': self.log_p,
            'NumHAcceptors': self.num_h_acceptors,
            'NumHDonors': self.num_h_donors,
            'CalcTPSA': self.calc_tpsa,
            'NumValenceElectrons': self.num_valence_electrons,
            'NumRotatableBonds': self.num_rotatable_bonds,
            'Charge': self.charge
        }
    
    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of the current configuration.
        
        Returns:
            Dict[str, Any]: Configuration summary
        """
        active_props = self.get_active_properties()
        num_active = sum(active_props.values())
        
        return {
            'num_active_properties': num_active,
            'total_properties': len(active_props),
            'active_properties': [name for name, active in active_props.items() if active],
            'similarity_constraints': {
                'max_tc': self.max_tc,
                'threshold': self.threshold
            },
            'decoy_limits': {
                'max_decoys_per_ligand': self.max_decoys_per_ligand
            }
        }


def create_config_file(output_path: str = "decoy_generation.json", 
                      config_type: str = "default") -> DecoyGenerationConfig:
    """
    Create and save a configuration file.
    
    Args:
        output_path: Path to save the configuration file
        config_type: Type of configuration ("default", "strict", "permissive")
        
    Returns:
        DecoyGenerationConfig: Created configuration object
        
    Raises:
        ValueError: If config_type is not recognized
    """
    if config_type == "default":
        config = DecoyGenerationConfig.default_config()
    elif config_type == "strict":
        config = DecoyGenerationConfig.strict_config()
    elif config_type == "permissive":
        config = DecoyGenerationConfig.permissive_config()
    else:
        raise ValueError(f"Unknown config type: {config_type}. Use 'default', 'strict', or 'permissive'")
    
    config.save_to_file(output_path)
    return config


def load_or_create_config(config_path: str = "decoy_generation.json") -> DecoyGenerationConfig:
    """
    Load configuration from file, or create default if file doesn't exist.
    
    Args:
        config_path: Path to the configuration file
        
    Returns:
        DecoyGenerationConfig: Loaded or created configuration
    """
    try:
        return DecoyGenerationConfig.load_from_file(config_path)
    except FileNotFoundError:
        print(f"Configuration file not found at {config_path}. Creating default configuration.")
        config = DecoyGenerationConfig.default_config()
        config.save_to_file(config_path)
        return config


if __name__ == "__main__":
    """Create default configuration file when run as script."""
    config = create_config_file()
    print(f"Created default configuration file: decoy_generation.json")
    print(f"Configuration summary:")
    summary = config.get_summary()
    for key, value in summary.items():
        print(f"  {key}: {value}")