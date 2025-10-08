# Decoy Generation Package Usage

## Quick Start

### Basic Command
```bash
python main.py \
    --enamine-std /path/to/sampled_std.npy \
    --prop-bundle-dir /path/to/property_bundles \
    --smi-bundle-dir /path/to/smiles_bundles \
    --ligands /path/to/ligands.smi \
    --config examples/default_config.json \
    --output-dir results \
    --verbose
```

### Help
```bash
python main.py --help
```

## Example Usage

### Full Database Run
```bash
python main.py \
    --enamine-std examples/sampled_std.npy \
    --prop-bundle-dir /path/to/property_bundles \
    --smi-bundle-dir /path/to/smiles_bundles \
    --ligands examples/ligands.smi \
    --config examples/default_config.json \
    --output-dir results \
    --n-proc 8 \
    --batch-size 50000 \
    --verbose
```

### Testing with Limited Data
```bash
python main.py \
    --enamine-std examples/sampled_std.npy \
    --prop-bundle-dir /path/to/property_bundles \
    --smi-bundle-dir /path/to/smiles_bundles \
    --ligands examples/ligands.smi \
    --config examples/default_config.json \
    --database-fraction 0.1 \
    --n-proc 4 \
    --verbose
```

### High Performance Run
```bash
python main.py \
    --enamine-std /path/to/sampled_std.npy \
    --prop-bundle-dir /path/to/property_bundles \
    --smi-bundle-dir /path/to/smiles_bundles \
    --ligands /path/to/ligands.smi \
    --config examples/default_config.json \
    --output-dir results \
    --oversample 20 \
    --n-proc 16 \
    --batch-size 100000 \
    --database-fraction 1.0 \
    --verbose
```

## Key Parameters

### Required Parameters
- `--enamine-std`: Path to standard deviations file (.npy)
- `--prop-bundle-dir`: Directory containing property bundle files
- `--smi-bundle-dir`: Directory containing SMILES bundle files  
- `--ligands`: File containing target ligand SMILES (format: "SMILES NAME")
- `--config`: Configuration file path (JSON format)

### Performance Parameters
- `--database-fraction`: Fraction of database to use (0.0-1.0, default: 1.0)
  - `0.1` = Use 10% of each bundle (fast testing)
  - `1.0` = Use full database (complete run)
- `--n-proc`: Number of parallel processes (default: 8)
- `--batch-size`: Molecules per batch (default: 10000)
- `--oversample`: Oversampling factor (default: 10)

### Output Parameters
- `--output-dir`: Output directory for results (default: "results")
- `--verbose`: Enable detailed progress output
- `--no-report`: Skip PDF report generation
- `--no-similarity-filter`: Disable structural similarity filtering

## Configuration File

The configuration file controls molecular property selection and filtering:

```json
{
    "Protonate": false,          // Apply pH 7.4 protonation to input molecules
    "MWt": false,               // Include molecular weight
    "NumHeavyAtoms": true,      // Include heavy atom count
    "logP": true,               // Include lipophilicity
    "NumHAcceptors": true,      // Include H-bond acceptors
    "NumHDonors": true,         // Include H-bond donors
    "CalcTPSA": false,          // Include topological polar surface area
    "NumValenceElectrons": false, // Include valence electrons
    "NumRotatableBonds": false,  // Include rotatable bonds
    "Charge": true,             // Filter by charge distribution [pos, neg]
    "Max_tc": 0.35,            // Structural similarity cutoff (Tanimoto)
    "Max_decoys_per_ligand": 50, // Final number of decoys per ligand
    "Threshold": 0.2           // Max property difference threshold
}
```

### Important Configuration Notes

- **"Protonate"**: When `true`, applies pH 7.4 protonation to determine ionization states
- **"Charge"**: Filters by charge distribution `[has_positive, has_negative]`:
  - `[1,0]` = Cations only
  - `[0,1]` = Anions only  
  - `[1,1]` = Zwitterions
  - `[0,0]` = Neutral molecules
- **"Max_tc"**: Tanimoto similarity cutoff for structural diversity (0.0-1.0)
- **"Threshold"**: Maximum allowed property difference (lower = stricter matching)

## Other Utilities

### Sample Properties
```bash
python sample_props.py \
    --bundle-dir /path/to/property_bundles \
    --database-fraction 0.1 \
    --output sampled_std.npy
```

### Create Configuration Examples  
```bash
python write_input.py
```

This creates example configuration files in the `examples/` directory.