# Decoy Generation Package

A clean, modular Python package for generating molecular decoys based on property similarity.

## Overview

This package provides a complete pipeline for generating molecular decoys for drug discovery applications. 

## Project Structure

```
decoy_generation/
├── __init__.py          # Package initialization
├── config.py            # Configuration management
├── database.py          # Database handling and property calculations
├── generator.py         # Main decoy generation logic
├── models.py            # Data structures and models
├── utils.py             # Utility functions for molecular properties
├── main.py              # Command-line interface
├── sample_props.py      # Property sampling utility
└── README.md            # This file
```

## Installation

Ensure you have the required dependencies:

```bash
pip install rdkit pandas numpy openbabel matplotlib
```

**Note**: `matplotlib` is optional but recommended for generating property comparison reports.

## Usage

### Command Line Interface

The main script provides a comprehensive command-line interface:

```bash
python main.py \
    --enamine-std /path/to/sampled_std.npy \
    --prop-bundle-dir /path/to/property_bundles \
    --smi-bundle-dir /path/to/smiles_bundles \
    --ligands /path/to/ligands.smi \
    --config /path/to/config.json \
    --output-dir results \
    --n-proc 8 \
    --verbose
```

### Python API

You can also use the package programmatically:

```python
from decoy_generation import DecoysGenerator, DecoyGenerationConfig

# Create configuration
config = DecoyGenerationConfig.default_config()
config.save_to_file("my_config.json")

# Initialize generator
generator = DecoysGenerator(
    enamine_std_path="sampled_std.npy",
    prop_bundle_dir="property_bundles",
    smi_bundle_dir="smiles_bundles", 
    prop_prefix="properties_",
    smi_prefix="smiles_",
    config_file="my_config.json",
    ligand_smi_file="ligands.smi"
)

# Generate decoys
results = generator.generate_decoys()
final_results = generator.finalize_decoys()

# Access results
for name, molecule in final_results.items():
    print(f"{name}: {molecule.num_decoys} decoys generated")
    print(molecule.decoys.head())
```

### Configuration

Create different types of configurations:

```python
from decoy_generation import DecoyGenerationConfig

# Default configuration
config = DecoyGenerationConfig.default_config()

# Strict configuration (more filters)
strict_config = DecoyGenerationConfig.strict_config()

# Permissive configuration (fewer filters)
permissive_config = DecoyGenerationConfig.permissive_config()

# Custom configuration
custom_config = DecoyGenerationConfig(
    num_heavy_atoms=True,
    log_p=True,
    charge=True,
    max_decoys_per_ligand=100,
    threshold=0.25
)

# Save to file
config.save_to_file("config.json")

# Load from file
loaded_config = DecoyGenerationConfig.load_from_file("config.json")
```

### Property Sampling

Generate standard deviations for normalization:

```bash
python sample_props.py \
    --bundle-dir /path/to/property_bundles \
    --prop-prefix properties_ \
    --sample-fraction 1e-7 \
    --output sampled_std.npy \
    --database-fraction 0.1
```

## File Formats

### Ligand File Format
```
SMILES_STRING MOLECULE_NAME
CCO ethanol
CC(C)O isopropanol
```

### Configuration File Format
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
    "Threshold": 0.2
}
```

## Molecular Properties

The package calculates the following molecular properties:

1. **NumHeavyAtoms**: Number of heavy (non-hydrogen) atoms
2. **MWt**: Molecular weight
3. **logP**: Lipophilicity (partition coefficient)
4. **NumHAcceptors**: Number of hydrogen bond acceptors
5. **NumHDonors**: Number of hydrogen bond donors
6. **CalcTPSA**: Topological polar surface area
7. **NumValenceElectrons**: Number of valence electrons
8. **NumRotatableBonds**: Number of rotatable bonds

## Output Files

The decoy generation process creates several output files in the specified output directory:

### CSV Files
- **`{molecule_name}_decoys.csv`**: One file per target molecule containing:
  - `smi`: SMILES string of the decoy molecule
  - `name`: Name/ID of the decoy molecule
  - `score`: Similarity score (lower is better)
  - Active molecular properties (e.g., `NumHeavyAtoms`, `logP`, `NumHAcceptors`, `NumHDonors`)

### Property Comparison Report
- **`property_report.pdf`**: Multi-page PDF with detailed property analysis:
  - **Spider/Radar plots**: Visual comparison showing actual property values
  - **Property value table**: Exact numerical values for target vs. decoys
  - **Statistical summary**: Mean, standard deviation, and range of decoy properties
  - **Annotations**: Target (T) and Decoy mean (D) values shown directly on plots
  - One page per target molecule with both visualization and data table

### Summary
- **`summary.json`**: JSON file with statistics about the generation process

Use `--no-report` to disable PDF report generation if matplotlib is not available.

## Advanced Usage

### Custom Property Calculations

You can extend the property calculations by modifying `utils.py`:

```python
from rdkit.Chem import Descriptors

# Add to MOLECULAR_PROPERTIES list
MOLECULAR_PROPERTIES.append(('MyProperty', Descriptors.MyDescriptor))
```

### Custom Filtering

Implement custom filtering logic in the generator:

```python
def custom_filter(df, molecule, config):
    # Your custom filtering logic here
    return df[df.custom_score > threshold]
```

## Performance Considerations

- **Memory Usage**: The package uses memory-mapped arrays to handle large datasets efficiently
- **Parallel Processing**: Adjust `n_proc` parameter based on your system capabilities
- **Batch Size**: Larger batch sizes improve performance but use more memory

## Troubleshooting

### Debug Mode

Use verbose mode for detailed output:

```bash
python main.py --verbose [other options]
```
