# Simple Decoy Finder

A clean, minimal implementation for finding property-matched molecular decoys from large enumerated libraries.

## Overview

This tool processes 2B+ molecule SMILES bundles with precomputed properties to find decoys that match target molecules based on user-specified molecular properties.

## Key Features

- **Large-scale processing**: Handles billion-molecule SMILES bundles efficiently
- **Memory mapping**: Uses numpy memory mapping for large property files
- **Multiprocessing**: Parallel processing of multiple bundles
- **Adaptive thresholding**: Automatically tightens similarity thresholds as better decoys are found
- **Charge matching**: Uses OpenBabel for on-demand charge state calculation
- **Flexible property selection**: User configures which properties to match

## Usage

```bash
python decoy_finder.py \
    --targets targets.smi \
    --smi-dir /path/to/smiles/bundles \
    --prop-dir /path/to/property/bundles \
    --config config.json \
    --stds sampled_std.npy \
    --n-proc 8
```

## File Formats

### Target File (targets.smi)
```
CCO ethanol
c1ccccc1 benzene
CC(C)O isopropanol
```

### Configuration (config.json)
```json
{
    "NumHeavyAtoms": true,
    "MWt": false,
    "logP": true,
    "NumHAcceptors": true,
    "NumHDonors": true,
    "CalcTPSA": false,
    "NumValenceElectrons": false,
    "NumRotatableBonds": false,
    "Charge": true,
    "Max_tc": 0.35,
    "Max_decoys_per_ligand": 50,
    "Threshold": 0.2,
    "oversample": 10
}
```

### Bundle Files
- SMILES: `{smi_prefix}{i}.smi` (e.g., `smiles_0.smi`)
- Properties: `{prop_prefix}{i}.npy` (e.g., `props_0.npy`)

## Properties Calculated

1. NumHeavyAtoms
2. MolWt (Molecular Weight)
3. logP (Lipophilicity)
4. NumHAcceptors (H-bond acceptors)
5. NumHDonors (H-bond donors)
6. CalcTPSA (Topological Polar Surface Area)
7. NumValenceElectrons
8. NumRotatableBonds

## Output

Creates CSV files: `bundle_{id}_{target_name}_decoys.csv`

Each file contains:
- `smi`: SMILES string
- `name`: Molecule name
- `score`: Similarity score (lower = better)
- Active properties as specified in config

## Algorithm

1. **Target Loading**: Load target molecules and calculate all properties
2. **Bundle Processing**: For each bundle, process in parallel using multiprocessing
3. **Property Matching**: Calculate similarity scores using normalized property differences
4. **Charge Matching**: Use OpenBabel to match charge states if enabled
5. **Adaptive Thresholding**: Tighten thresholds as better decoys are found
6. **Output Generation**: Write best decoys to CSV files (one per bundle per target)

## Performance

- Designed for 2B+ molecule datasets
- Memory-efficient through numpy memory mapping
- Parallel processing across bundles
- Batch processing within bundles to manage memory