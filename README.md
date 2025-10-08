# Decoy Generation Package

A clean, modular Python package for generating molecular decoys based on property similarity.

## Overview

This package provides a complete pipeline for generating molecular decoys for drug discovery applications. It features property-based filtering, structural similarity analysis, enhanced charge detection, and comprehensive PDF reporting.

## Features

- **Modular Design**: Clean separation of concerns with dedicated modules for configuration, database handling, and generation logic
- **Property-Based Filtering**: Filter molecules by molecular weight, LogP, hydrogen bonding, TPSA, and more
- **Structural Similarity**: Tanimoto coefficient-based filtering for chemical diversity
- **Enhanced Charge Detection**: Detects cations, anions, zwitterions, and neutral molecules with pH-dependent protonation
- **PDF Reports**: Automatic generation of spider plots and data tables showing property distributions
- **Configurable Pipeline**: JSON-based configuration for easy customization
- **Performance Optimized**: Multi-processing support and memory-efficient batch processing

## Installation

### Option 1: Install from GitHub
```bash
# Clone and install in one step
pip install git+https://github.com/lmoesgaard/decoy_generation.git

# Or clone first, then install
git clone https://github.com/lmoesgaard/decoy_generation.git
cd decoy_generation
pip install .
```

### Option 2: Development Installation
```bash
# Clone the repository
git clone https://github.com/lmoesgaard/decoy_generation.git
cd decoy_generation

# Install in development mode (recommended for contributors)
pip install -e .

# Or install with optional dependencies
pip install -e .[openbabel,dev]
```

### Option 3: Manual Dependency Installation
```bash
# If you prefer to install dependencies separately
pip install -r requirements.txt
python main.py --help  # Use directly without pip install
```

### Dependencies

**Core dependencies** (automatically installed):
- **rdkit**: Molecular informatics toolkit
- **pandas**: Data manipulation and analysis  
- **numpy**: Numerical computing
- **matplotlib**: Plotting and visualization
- **scipy**: Scientific computing

**Optional dependencies**:
- **openbabel-wheel**: Enhanced charge detection with pH-dependent protonation
  ```bash
  pip install openbabel-wheel
  # or install with: pip install .[openbabel]
  ```

### Verify Installation

```bash
# Test the installation
python -c "from main import main; print('Installation successful!')"

# Test with help
python main.py --help
```

## Usage

After installation, you can use the package in several ways:

### Option 1: Direct Script Execution (Recommended)

```bash
# Generate decoys
python main.py --ligands ligands.smi --config examples/default_config.json --output-dir results

# Sample properties 
python sample_props.py --bundle-dir /path/to/bundles --output sampled_std.npy

# Get help
python main.py --help
```

### Option 2: Python Import

```bash
# From the package directory
python -c "from main import main; main()" --ligands ligands.smi --config config.json
```

### Example Command

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
