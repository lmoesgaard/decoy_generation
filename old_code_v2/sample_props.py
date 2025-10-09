"""
Sample properties from molecular property database to compute standard deviations.

This script samples a small fraction of molecules from each property bundle
to compute standard deviations for normalization purposes in decoy generation.
"""

import argparse
import numpy as np
from numpy.lib.format import open_memmap
from pathlib import Path
from typing import Optional


def sample_properties(bundle_dir: str, prop_prefix: str = "properties_", 
                     max_bundles: Optional[int] = None, 
                     sample_fraction: float = 1e-7,
                     output_file: str = "sampled_std.npy",
                     random_seed: Optional[int] = None) -> np.ndarray:
    """
    Sample properties from molecular property bundles and compute standard deviations.
    
    Args:
        bundle_dir: Directory containing property bundle files
        prop_prefix: Prefix for property bundle files
        max_bundles: Maximum number of bundles to process (None for all available)
        sample_fraction: Fraction of molecules to sample from each bundle
        output_file: Path to save the computed standard deviations
        random_seed: Random seed for reproducibility
        
    Returns:
        np.ndarray: Standard deviations for each property
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    
    bundle_path = Path(bundle_dir)
    if not bundle_path.exists():
        raise FileNotFoundError(f"Bundle directory not found: {bundle_dir}")
    
    combined_arrays = []
    bundle_count = 0
    bundle_index = 0
    
    print(f"Sampling properties from {bundle_dir}")
    print(f"Sample fraction: {sample_fraction:.2e}")
    if max_bundles:
        print(f"Max bundles: {max_bundles}")
    
    while max_bundles is None or bundle_count < max_bundles:
        prop_file = bundle_path / f"{prop_prefix}{bundle_index}.npy"
        
        if not prop_file.exists():
            if bundle_count == 0:
                print(f"No property files found with prefix '{prop_prefix}'")
                return np.array([])
            print(f"Bundle {bundle_index} not found, stopping.")
            break
        
        try:
            # Load property array as memory map
            arr = open_memmap(str(prop_file), mode="r")
            
            # Calculate number of samples
            n_samples = max(1, int(arr.shape[0] * sample_fraction))
            
            # Generate random indices
            idx = np.random.randint(0, arr.shape[0], size=n_samples, dtype=np.int64)
            
            # Sample and copy to avoid memory map issues
            sampled = arr[idx].copy()
            combined_arrays.append(sampled)
            
            print(f"Bundle {bundle_index}: sampled {n_samples}/{arr.shape[0]} molecules "
                  f"({arr.shape[1]} properties)")
            
            bundle_count += 1
            bundle_index += 1
            
        except Exception as e:
            print(f"Error processing bundle {bundle_index}: {e}")
            break
    
    if not combined_arrays:
        raise ValueError("No property data could be loaded")
    
    # Combine all samples
    print(f"Combining samples from {bundle_count} bundles...")
    combined_arr = np.vstack(combined_arrays)
    total_samples = combined_arr.shape[0]
    
    # Compute standard deviations
    print(f"Computing standard deviations for {total_samples} total samples...")
    combined_std = combined_arr.std(axis=0)
    
    # Save results
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(output_file, combined_std)
    
    print(f"Saved standard deviations to {output_file}")
    print(f"Standard deviations shape: {combined_std.shape}")
    print(f"Standard deviations: {combined_std}")
    
    return combined_std


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Sample molecular properties to compute standard deviations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--bundle-dir",
        default="../xreal_v0p0_2B_fplibrary",
        help="Directory containing property bundle files"
    )
    parser.add_argument(
        "--prop-prefix",
        default="properties_",
        help="Prefix for property bundle files"
    )
    parser.add_argument(
        "--max-bundles",
        type=int,
        help="Maximum number of bundles to process"
    )
    parser.add_argument(
        "--sample-fraction",
        type=float,
        default=1e-7,
        help="Fraction of molecules to sample from each bundle"
    )
    parser.add_argument(
        "--output",
        default="sampled_std.npy",
        help="Output file for standard deviations"
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    args = parser.parse_args()
    
    try:
        std_devs = sample_properties(
            bundle_dir=args.bundle_dir,
            prop_prefix=args.prop_prefix,
            max_bundles=args.max_bundles,
            sample_fraction=args.sample_fraction,
            output_file=args.output,
            random_seed=args.random_seed
        )
        
        if args.verbose:
            print(f"\nDetailed results:")
            for i, std in enumerate(std_devs):
                print(f"  Property {i}: {std:.6f}")
        
        print(f"\nSuccess! Standard deviations computed and saved.")
        
    except Exception as e:
        print(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())