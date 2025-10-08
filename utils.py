"""
Molecular property calculation utilities for decoy generation.

This module contains utility functions for calculating molecular properties
and handling SMILES strings, including charge prediction and molecular
property definitions.
"""

from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski
from rdkit import DataStructs, Chem
from rdkit.Chem import rdFingerprintGenerator
from openbabel import openbabel as ob
from typing import Tuple, List, Callable
import numpy as np


def get_charge(smi: str) -> List[int]:
    """
    Get the charge distribution of a molecule from its SMILES string without pH correction.
    
    Args:
        smi: SMILES string representation of the molecule
        
    Returns:
        List[int]: [has_positive_charge, has_negative_charge]
                   [1,0] = cation only
                   [0,1] = anion only  
                   [1,1] = zwitterion
                   [0,0] = neutral
        
    Raises:
        ValueError: If the SMILES string is invalid
    """
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "smi")

    mol = ob.OBMol()
    if not conv.ReadString(mol, smi):
        raise ValueError(f"Invalid SMILES: {smi}")

    # No pH correction - use molecule as-is
    has_positive = 0
    has_negative = 0
    
    # Check each atom for formal charge, excluding nitro groups
    for atom in ob.OBMolAtomIter(mol):
        formal_charge = atom.GetFormalCharge()
        atomic_num = atom.GetAtomicNum()
        
        # Skip nitro group charges (N with +1 connected to O with -1)
        if formal_charge > 0 and atomic_num == 7:  # Nitrogen with positive charge
            # Check if this is part of a nitro group
            is_nitro = False
            for neighbor in ob.OBAtomAtomIter(atom):
                if (neighbor.GetAtomicNum() == 8 and  # Oxygen
                    neighbor.GetFormalCharge() < 0):  # With negative charge
                    is_nitro = True
                    break
            if not is_nitro:
                has_positive = 1
        elif formal_charge > 0:
            has_positive = 1
        elif formal_charge < 0 and atomic_num == 8:  # Oxygen with negative charge
            # Check if this is part of a nitro group
            is_nitro = False
            for neighbor in ob.OBAtomAtomIter(atom):
                if (neighbor.GetAtomicNum() == 7 and  # Nitrogen
                    neighbor.GetFormalCharge() > 0):  # With positive charge
                    is_nitro = True
                    break
            if not is_nitro:
                has_negative = 1
        elif formal_charge < 0:
            has_negative = 1
    
    return [has_positive, has_negative]


def predict_charge(smi: str, pH: float = 7.4) -> List[int]:
    """
    Predict the charge distribution of a molecule at a given pH.
    
    This function protonates a molecule at the specified pH and returns
    the presence of positive and negative charges as [pos, neg].
    Nitro groups are excluded as they have localized +/- charges.
    
    Args:
        smi: SMILES string representation of the molecule
        pH: pH value for protonation (default: 7.4, physiological pH)
        
    Returns:
        List[int]: [has_positive_charge, has_negative_charge]
                   [1,0] = cation only
                   [0,1] = anion only  
                   [1,1] = zwitterion
                   [0,0] = neutral
        
    Raises:
        ValueError: If the SMILES string is invalid
    """
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "smi")

    mol = ob.OBMol()
    if not conv.ReadString(mol, smi):
        raise ValueError(f"Invalid SMILES: {smi}")

    mol.CorrectForPH(pH)
    
    has_positive = 0
    has_negative = 0
    
    # Check each atom for formal charge, excluding nitro groups
    for atom in ob.OBMolAtomIter(mol):
        formal_charge = atom.GetFormalCharge()
        atomic_num = atom.GetAtomicNum()
        
        # Skip nitro group charges (N with +1 connected to O with -1)
        if formal_charge > 0 and atomic_num == 7:  # Nitrogen with positive charge
            # Check if this is part of a nitro group
            is_nitro = False
            for neighbor in ob.OBAtomAtomIter(atom):
                if (neighbor.GetAtomicNum() == 8 and  # Oxygen
                    neighbor.GetFormalCharge() < 0):  # With negative charge
                    is_nitro = True
                    break
            if not is_nitro:
                has_positive = 1
        elif formal_charge > 0:
            has_positive = 1
        elif formal_charge < 0 and atomic_num == 8:  # Oxygen with negative charge
            # Check if this is part of a nitro group
            is_nitro = False
            for neighbor in ob.OBAtomAtomIter(atom):
                if (neighbor.GetAtomicNum() == 7 and  # Nitrogen
                    neighbor.GetFormalCharge() > 0):  # With positive charge
                    is_nitro = True
                    break
            if not is_nitro:
                has_negative = 1
        elif formal_charge < 0:
            has_negative = 1
    
    return [has_positive, has_negative]


# Property definitions with their corresponding calculation functions
# Order is important and must match the database schema
MOLECULAR_PROPERTIES: List[Tuple[str, Callable]] = [
    ('NumHeavyAtoms', rdMolDescriptors.CalcNumHeavyAtoms),
    ('MWt', Descriptors.MolWt),
    ('logP', Crippen.MolLogP),
    ('NumHAcceptors', Lipinski.NumHAcceptors),
    ('NumHDonors', Lipinski.NumHDonors),
    ('CalcTPSA', rdMolDescriptors.CalcTPSA),
    ('NumValenceElectrons', Descriptors.NumValenceElectrons),
    ('NumRotatableBonds', Lipinski.NumRotatableBonds)
]


def get_property_names() -> List[str]:
    """
    Get the list of property names in the correct order.
    
    Returns:
        List[str]: Property names in database order
    """
    return [name for name, _ in MOLECULAR_PROPERTIES]


def get_property_calculators() -> List[Callable]:
    """
    Get the list of property calculation functions in the correct order.
    
    Returns:
        List[Callable]: Property calculation functions in database order
    """
    return [func for _, func in MOLECULAR_PROPERTIES]


def smi2fp(smi: str):
    """Convert SMILES to Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024).GetFingerprint(mol)
    return fp


def sim_filter(smi_lst: list, pop_size: int, cutoff: float = 0.35) -> np.ndarray:
    """
    Filter SMILES list to remove similar molecules.
    
    Args:
        smi_lst: List of SMILES strings
        pop_size: Maximum number of molecules to keep
        cutoff: Tanimoto similarity cutoff (molecules above this similarity are filtered out)
        
    Returns:
        np.ndarray: Boolean mask indicating which molecules to keep
    """
    if len(smi_lst) == 0:
        return np.array([], dtype=bool)
    
    if len(smi_lst) <= pop_size:
        return np.ones(len(smi_lst), dtype=bool)
    
    # Start with first molecule
    fps = []
    first_fp = smi2fp(smi_lst[0])
    if first_fp is not None:
        fps.append(first_fp)
    
    mask = np.zeros(len(smi_lst), dtype=bool)
    mask[0] = True
    ind = 1
    
    while mask.sum() < pop_size and ind < len(smi_lst):
        fp = smi2fp(smi_lst[ind])
        if fp is not None and len(fps) > 0:
            sim = np.max(DataStructs.BulkTanimotoSimilarity(fp, fps))
            if sim < cutoff:
                fps.append(fp)
                mask[ind] = True
        elif fp is not None and len(fps) == 0:
            # First valid fingerprint
            fps.append(fp)
            mask[ind] = True
        ind += 1
    
    return mask