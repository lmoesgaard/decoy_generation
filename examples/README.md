# Example Configuration Files

This directory contains example JSON configuration files for the decoy generation package.

## Configuration Files

### `default_config.json`
The default configuration that balances filtering and permissiveness:
- Uses basic molecular properties (heavy atoms, logP, H-bond donors/acceptors, charge)
- Moderate similarity thresholds
- Standard decoy count limits

### `strict_config.json` 
A more restrictive configuration that enables all property filters:
- All molecular properties enabled
- Lower similarity thresholds (more restrictive)
- Fewer decoys per ligand

### `permissive_config.json`
A more lenient configuration with minimal filtering:
- Only basic properties enabled (heavy atoms, logP)
- Higher similarity thresholds (less restrictive)
- More decoys per ligand allowed

## Configuration Parameters

### Molecular Property Filters
- `Protonate`: Apply protonation at physiological pH
- `MWt`: Filter by molecular weight  
- `NumHeavyAtoms`: Filter by number of heavy atoms
- `logP`: Filter by lipophilicity
- `NumHAcceptors`: Filter by hydrogen bond acceptors
- `NumHDonors`: Filter by hydrogen bond donors
- `CalcTPSA`: Filter by topological polar surface area
- `NumValenceElectrons`: Filter by valence electrons
- `NumRotatableBonds`: Filter by rotatable bonds
- `Charge`: Filter by charge distribution pattern (see below)

### Charge Filtering
The `Charge` parameter filters by ionization state pattern `[has_positive, has_negative]`:
- `[1,0]` = Cations only (molecules with positive charges)
- `[0,1]` = Anions only (molecules with negative charges)  
- `[1,1]` = Zwitterions (molecules with both positive and negative charges)
- `[0,0]` = Neutral molecules (no ionizable groups)

### Protonation Control
- `Protonate`: Apply pH 7.4 protonation to determine ionization states
  - `true`: Considers pH-dependent ionization (amines protonate, acids deprotonate)
  - `false`: Uses SMILES charge states as-is

### Similarity Constraints
- `Max_tc`: Structural similarity cutoff using Tanimoto coefficient (0.0-1.0)
  - Lower values = more diverse decoys, fewer total decoys
  - Higher values = less diverse decoys, more total decoys
- `Threshold`: Maximum property difference threshold 
  - Lower values = stricter property matching
  - Higher values = more permissive property matching
- `Max_decoys_per_ligand`: Final number of decoys to generate per ligand

## Usage

Copy one of these files and modify the parameters as needed:

```bash
cp examples/default_config.json my_config.json
# Edit my_config.json as needed
python main.py --config my_config.json [other options]
```

Or create programmatically:

```python
from config import DecoyGenerationConfig

# Load and modify
config = DecoyGenerationConfig.load_from_file("examples/default_config.json")
config.max_decoys_per_ligand = 75
config.save_to_file("my_config.json")
```