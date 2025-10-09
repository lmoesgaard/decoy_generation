"""
Data models for decoy generation.

This module defines the data structures used throughout the decoy generation
process, including the Molecule class and related data types.
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Optional


@dataclass
class Molecule:
    """
    Represents a molecule with its properties and associated decoys.
    
    This class stores all the information needed for a molecule during
    the decoy generation process, including its SMILES representation,
    calculated properties, and generated decoys.
    
    Attributes:
        smiles: SMILES string representation of the molecule
        charge: Formal charge of the molecule
        properties: Array of calculated molecular properties
        threshold: Threshold value for decoy similarity filtering
        decoys: DataFrame containing generated decoy molecules
    """
    smiles: str
    charge: int
    properties: np.ndarray
    threshold: float
    decoys: pd.DataFrame
    
    def __post_init__(self):
        """Validate the molecule data after initialization."""
        if not isinstance(self.smiles, str) or not self.smiles.strip():
            raise ValueError("SMILES must be a non-empty string")
        
        if not isinstance(self.properties, np.ndarray):
            raise ValueError("Properties must be a numpy array")
        
        if not isinstance(self.decoys, pd.DataFrame):
            raise ValueError("Decoys must be a pandas DataFrame")
    
    @property
    def num_decoys(self) -> int:
        """Get the current number of decoys for this molecule."""
        return len(self.decoys)
    
    def add_decoys(self, new_decoys: pd.DataFrame) -> None:
        """
        Add new decoys to the existing decoy set.
        
        Args:
            new_decoys: DataFrame containing new decoy molecules to add
        """
        if not isinstance(new_decoys, pd.DataFrame):
            raise ValueError("New decoys must be a pandas DataFrame")
        
        self.decoys = pd.concat([self.decoys, new_decoys], ignore_index=True)
    
    def filter_decoys_by_score(self, max_score: float) -> None:
        """
        Filter decoys to keep only those with score below the threshold.
        
        Args:
            max_score: Maximum allowed score for decoys
        """
        if 'score' not in self.decoys.columns:
            raise ValueError("Decoys DataFrame must contain a 'score' column")
        
        self.decoys = self.decoys[self.decoys['score'] < max_score].copy()
    
    def sample_decoys(self, n_samples: int, weight_column: str = 'score') -> None:
        """
        Sample a subset of decoys, optionally with weights.
        
        Args:
            n_samples: Number of decoys to sample
            weight_column: Column name to use for weighted sampling (optional)
        """
        if n_samples >= len(self.decoys):
            return  # No sampling needed
        
        if weight_column in self.decoys.columns:
            # Weighted sampling (higher scores get lower weights)
            weights = self.decoys[weight_column].max() - self.decoys[weight_column] + 1e-6
            self.decoys = self.decoys.sample(n=n_samples, weights=weights)
            # Update threshold to the maximum score in the sampled set
            self.threshold = self.decoys[weight_column].max()
        else:
            # Random sampling
            self.decoys = self.decoys.sample(n=n_samples)
    
    def get_summary(self) -> dict:
        """
        Get a summary of the molecule and its decoys.
        
        Returns:
            dict: Summary information about the molecule
        """
        summary = {
            'smiles': self.smiles,
            'charge': self.charge,
            'num_properties': len(self.properties),
            'threshold': self.threshold,
            'num_decoys': self.num_decoys
        }
        
        if self.num_decoys > 0 and 'score' in self.decoys.columns:
            summary.update({
                'min_score': self.decoys['score'].min(),
                'max_score': self.decoys['score'].max(),
                'mean_score': self.decoys['score'].mean()
            })
        
        return summary