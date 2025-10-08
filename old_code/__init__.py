"""
Decoy Generation Package

A package for generating molecular decoys based on property similarity.

This package provides tools for:
- Loading and processing molecular property databases
- Calculating molecular properties and similarities
- Generating decoy molecules for drug discovery applications
- Managing configuration and running the complete pipeline

Main modules:
- config: Configuration management
- database: Database handling and property calculations
- generator: Main decoy generation logic
- models: Data structures and models
- utils: Utility functions for molecular properties
- main: Command-line interface
"""

from config import DecoyGenerationConfig, load_or_create_config, create_config_file
from database import DatabaseProperties
from generator import DecoysGenerator
from models import Molecule
from utils import get_charge, predict_charge, MOLECULAR_PROPERTIES

__version__ = "1.0.0"
__author__ = "Decoy Generation Team"

__all__ = [
    "DecoyGenerationConfig",
    "load_or_create_config", 
    "create_config_file",
    "DatabaseProperties",
    "DecoysGenerator",
    "Molecule",
    "get_charge",
    "predict_charge",
    "MOLECULAR_PROPERTIES"
]