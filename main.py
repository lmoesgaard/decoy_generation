"""
Main script for decoy generation.

This script provides a simple command-line interface for generating
molecular decoys using the decoy generation pipeline.
"""

import os
import sys
import argparse
import json
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from typing import Optional, Dict, List

from config import DecoyGenerationConfig, load_or_create_config
#from generator import DecoysGenerator
from parallel_integration import ParallelDecoysGenerator as DecoysGenerator
from models import Molecule


def save_results(results: dict, output_dir: str = "results") -> None:
    """
    Save decoy generation results to CSV files.
    
    Each CSV file contains:
    - smi: SMILES string of the decoy molecule
    - name: Name/ID of the decoy molecule
    - score: Similarity score (lower is better)
    - Active molecular properties (based on configuration)
    
    Args:
        results: Dictionary of molecule name -> Molecule objects
        output_dir: Directory to save results
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Save individual decoy files
    for name, molecule in results.items():
        if molecule.num_decoys > 0:
            decoy_file = output_path / f"{name}_decoys.csv"
            molecule.decoys.to_csv(decoy_file, index=False)
            print(f"Saved {molecule.num_decoys} decoys for {name} to {decoy_file}")
    
    # Save summary
    summary_data = {}
    for name, molecule in results.items():
        summary_data[name] = molecule.get_summary()
    
    summary_file = output_path / "summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary_data, f, indent=2)
    
    print(f"Saved summary to {summary_file}")


def create_spider_plot(ax, target_properties: np.ndarray, decoy_properties: np.ndarray, 
                      property_names: List[str], molecule_name: str) -> None:
    """
    Create a spider/radar plot comparing target and decoy properties using actual values.
    
    Args:
        ax: Matplotlib axis to plot on
        target_properties: Array of target molecule properties
        decoy_properties: 2D array of decoy properties (n_decoys x n_properties)
        property_names: List of property names
        molecule_name: Name of the target molecule
    """
    # Number of properties
    n_props = len(property_names)
    
    # Calculate angles for each property
    angles = np.linspace(0, 2 * np.pi, n_props, endpoint=False).tolist()
    angles += angles[:1]  # Complete the circle
    
    # Calculate statistics for each property
    property_stats = []
    for i in range(n_props):
        if decoy_properties.shape[0] > 0:
            decoy_col = decoy_properties[:, i].astype(float)  # Ensure float type
            stats = {
                'target': float(target_properties[i]),
                'decoy_mean': float(decoy_col.mean()),
                'decoy_std': float(decoy_col.std()),
                'decoy_min': float(decoy_col.min()),
                'decoy_max': float(decoy_col.max())
            }
        else:
            stats = {
                'target': float(target_properties[i]),
                'decoy_mean': float(target_properties[i]),
                'decoy_std': 0.0,
                'decoy_min': float(target_properties[i]),
                'decoy_max': float(target_properties[i])
            }
        property_stats.append(stats)
    
    # Determine scale for each property (use reasonable range around data)
    scales = []
    for i, stats in enumerate(property_stats):
        # Find the range that encompasses all data with some padding
        all_values = [stats['target'], stats['decoy_mean'], 
                     stats['decoy_mean'] - stats['decoy_std'],
                     stats['decoy_mean'] + stats['decoy_std']]
        
        data_min = min(all_values)
        data_max = max(all_values)
        
        # Add padding (20% on each side)
        data_range = data_max - data_min
        if data_range < 1e-10:  # Handle case where all values are the same
            data_range = abs(data_max) * 0.1 if data_max != 0 else 1.0
        
        scale_min = data_min - 0.2 * data_range
        scale_max = data_max + 0.2 * data_range
        
        # For some properties, ensure reasonable minimum values
        if property_names[i] in ['NumHeavyAtoms', 'NumHAcceptors', 'NumHDonors', 
                                'NumValenceElectrons', 'NumRotatableBonds']:
            scale_min = max(0, scale_min)  # These can't be negative
        
        scales.append((scale_min, scale_max))
    
    # Normalize values to 0-1 for plotting (but keep original values for labels)
    def normalize_value(value, scale_min, scale_max):
        # Ensure all values are floats to avoid type issues
        value = float(value)
        scale_min = float(scale_min)
        scale_max = float(scale_max)
        
        if scale_max - scale_min > 1e-10:
            return (value - scale_min) / (scale_max - scale_min)
        else:
            return 0.5
    
    # Prepare normalized data for plotting
    norm_target = []
    norm_decoy_mean = []
    norm_decoy_upper = []
    norm_decoy_lower = []
    
    for i, stats in enumerate(property_stats):
        scale_min, scale_max = scales[i]
        
        norm_target.append(normalize_value(stats['target'], scale_min, scale_max))
        norm_decoy_mean.append(normalize_value(stats['decoy_mean'], scale_min, scale_max))
        
        upper_val = stats['decoy_mean'] + stats['decoy_std']
        lower_val = stats['decoy_mean'] - stats['decoy_std']
        norm_decoy_upper.append(normalize_value(upper_val, scale_min, scale_max))
        norm_decoy_lower.append(normalize_value(lower_val, scale_min, scale_max))
    
    # Complete the circle for plotting
    norm_target.append(norm_target[0])
    norm_decoy_mean.append(norm_decoy_mean[0])
    norm_decoy_upper.append(norm_decoy_upper[0])
    norm_decoy_lower.append(norm_decoy_lower[0])
    
    # Create the plot
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    
    # Plot target molecule
    ax.plot(angles, norm_target, 'o-', linewidth=2, label=f'Target ({molecule_name})', 
           color='red', markersize=6)
    ax.fill(angles, norm_target, alpha=0.25, color='red')
    
    # Plot decoy mean
    ax.plot(angles, norm_decoy_mean, 'o-', linewidth=2, label='Decoy Mean', 
           color='blue', markersize=4)
    
    # Plot decoy std deviation as shaded area
    ax.fill_between(angles, norm_decoy_lower, norm_decoy_upper, alpha=0.2, color='blue', 
                   label='Decoy ±1 StdDev')
    
    # Customize the plot with actual values
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(property_names)
    ax.set_ylim(0, 1)
    
    # Create custom radial tick labels - show relative scale
    n_ticks = 5
    tick_positions = np.linspace(0, 1, n_ticks)
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(['Low', '', 'Med', '', 'High'])
    
    ax.grid(True)
    ax.legend(loc='lower left', bbox_to_anchor=(0.0, -0.1), ncol=3, fontsize=9)
    ax.set_title(f'{molecule_name} - Property Comparison', 
                 size=12, fontweight='bold', pad=20)


def generate_property_report(results: Dict[str, Molecule], output_dir: str = "results") -> None:
    """
    Generate a PDF report with spider plots comparing target and decoy properties.
    
    Args:
        results: Dictionary of molecule name -> Molecule objects
        output_dir: Directory to save the report
    """
    output_path = Path(output_dir)
    report_file = output_path / "property_report.pdf"
    
    # Skip if no molecules have decoys
    molecules_with_decoys = {name: mol for name, mol in results.items() if mol.num_decoys > 0}
    if not molecules_with_decoys:
        print("No molecules with decoys found. Skipping property report generation.")
        return
    
    # Get property names from the first molecule with decoys
    first_molecule = next(iter(molecules_with_decoys.values()))
    property_columns = [col for col in first_molecule.decoys.columns 
                       if col not in ['smi', 'smiles', 'name', 'score', 'target']]
    
    if not property_columns:
        print("No property columns found in decoy data. Skipping property report generation.")
        return
    
    try:
        with PdfPages(report_file) as pdf:
            for name, molecule in molecules_with_decoys.items():
                # Create figure with subplots: spider plot + table
                fig = plt.figure(figsize=(14, 10))
                
                # Spider plot on the left
                ax1 = plt.subplot(121, projection='polar')
                
                # Extract decoy properties
                decoy_data = molecule.decoys[property_columns].values.astype(float)
                
                # Get target properties (only the active ones)
                target_properties = np.array(molecule.properties, dtype=float)
                
                # Create spider plot
                create_spider_plot(ax1, target_properties, decoy_data, 
                                 property_columns, name)
                
                # Create summary table on the right
                ax2 = plt.subplot(122)
                ax2.axis('off')  # Hide axes for table
                
                # Calculate statistics for table
                table_data = []
                headers = ['Property', 'Target', 'Mean', 'StdDev', 'Range']
                
                for i, prop_name in enumerate(property_columns):
                    target_val = target_properties[i]
                    decoy_col = decoy_data[:, i]
                    decoy_mean = decoy_col.mean()
                    decoy_std = decoy_col.std()
                    decoy_min = decoy_col.min()
                    decoy_max = decoy_col.max()
                    
                    # Format values based on property type
                    if 'log' in prop_name.lower():
                        target_str = f"{target_val:.2f}"
                        mean_str = f"{decoy_mean:.2f}"
                        std_str = f"{decoy_std:.2f}"
                        range_str = f"{decoy_min:.2f}-{decoy_max:.2f}"
                    elif prop_name in ['NumHeavyAtoms', 'NumHAcceptors', 'NumHDonors', 
                                      'NumValenceElectrons', 'NumRotatableBonds']:
                        target_str = f"{target_val:.0f}"
                        mean_str = f"{decoy_mean:.1f}"
                        std_str = f"{decoy_std:.1f}"
                        range_str = f"{decoy_min:.0f}-{decoy_max:.0f}"
                    else:
                        target_str = f"{target_val:.1f}"
                        mean_str = f"{decoy_mean:.1f}"
                        std_str = f"{decoy_std:.1f}"
                        range_str = f"{decoy_min:.1f}-{decoy_max:.1f}"
                    
                    table_data.append([prop_name, target_str, mean_str, std_str, range_str])
                
                # Create table with better column widths
                table = ax2.table(cellText=table_data, colLabels=headers,
                                 cellLoc='center', loc='upper center',
                                 colWidths=[0.3, 0.15, 0.15, 0.15, 0.25])
                table.auto_set_font_size(False)
                table.set_fontsize(10)
                table.scale(1.0, 1.8)
                
                # Style the table
                for i in range(len(headers)):
                    table[(0, i)].set_facecolor('#40466e')
                    table[(0, i)].set_text_props(weight='bold', color='white')
                
                # Alternate row colors
                for i in range(1, len(table_data) + 1):
                    for j in range(len(headers)):
                        if i % 2 == 0:
                            table[(i, j)].set_facecolor('#f1f1f2')
                        else:
                            table[(i, j)].set_facecolor('white')
                
                ax2.set_title('Property Values Summary', fontsize=12, fontweight='bold', pad=5)
                
                # Add overall figure title and molecule info
                fig.suptitle(f'Property Analysis Report: {name}', fontsize=16, fontweight='bold')
                
                # Add molecule information at the bottom
                info_text = f"SMILES: {molecule.smiles}\n"
                info_text += f"Number of decoys: {molecule.num_decoys}\n"
                info_text += f"Max Property Difference: {molecule.threshold:.3f}"
                
                fig.text(0.5, 0.02, info_text, ha='center', va='bottom', fontsize=10,
                        bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
                
                plt.tight_layout()
                plt.subplots_adjust(top=0.9, bottom=0.15)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
        
        print(f"Generated property comparison report: {report_file}")
        
    except ImportError:
        print("Warning: matplotlib not available. Skipping property report generation.")
        print("Install matplotlib to generate property comparison plots: pip install matplotlib")
    except Exception as e:
        print(f"Error generating property report: {e}")


def validate_paths(args: argparse.Namespace) -> bool:
    """
    Validate that all required paths exist.
    
    Args:
        args: Parsed command line arguments
        
    Returns:
        bool: True if all paths are valid
    """
    required_paths = [
        (args.enamine_std, "Enamine standard deviations file"),
        (args.prop_bundle_dir, "Property bundle directory"),
        (args.smi_bundle_dir, "SMILES bundle directory"),
        (args.ligands, "Ligand SMILES file")
    ]
    
    all_valid = True
    for path, description in required_paths:
        if not os.path.exists(path):
            print(f"Error: {description} not found: {path}")
            all_valid = False
    
    return all_valid


def setup_config(config_path: Optional[str] = None) -> DecoyGenerationConfig:
    """
    Set up configuration for decoy generation.
    
    Args:
        config_path: Path to configuration file (optional)
        
    Returns:
        DecoyGenerationConfig: Configuration object
    """
    if config_path:
        try:
            config = DecoyGenerationConfig.load_from_file(config_path)
            print(f"Loaded configuration from {config_path}")
        except Exception as e:
            print(f"Error loading configuration: {e}")
            print("Using default configuration...")
            config = DecoyGenerationConfig.default_config()
    else:
        config = load_or_create_config()
    
    # Print configuration summary
    summary = config.get_summary()
    print("\nConfiguration Summary:")
    print(f"  Active properties: {summary['num_active_properties']}/{summary['total_properties']}")
    print(f"  Properties: {', '.join(summary['active_properties'])}")
    print(f"  Max TC: {summary['similarity_constraints']['max_tc']}")
    print(f"  Threshold: {summary['similarity_constraints']['threshold']}")
    print(f"  Max decoys per ligand: {summary['decoy_limits']['max_decoys_per_ligand']}")
    
    return config


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Generate molecular decoys for given ligands. "
                   "Outputs CSV files containing SMILES, names, similarity scores, "
                   "and molecular properties.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        "--enamine-std",
        required=True,
        help="Path to Enamine standard deviations file (.npy)"
    )
    parser.add_argument(
        "--prop-bundle-dir",
        required=True,
        help="Directory containing property bundle files"
    )
    parser.add_argument(
        "--smi-bundle-dir",
        required=True,
        help="Directory containing SMILES bundle files"
    )
    parser.add_argument(
        "--ligands",
        required=True,
        help="File containing ligand SMILES (format: SMILES NAME per line)"
    )
    
    # Optional arguments
    parser.add_argument(
        "--config",
        help="Path to configuration file (JSON format)"
    )
    parser.add_argument(
        "--prop-prefix",
        default="props_",
        help="Prefix for property bundle files"
    )
    parser.add_argument(
        "--smi-prefix", 
        default="smiles_",
        help="Prefix for SMILES bundle files"
    )
    parser.add_argument(
        "--output-dir",
        default="results",
        help="Directory to save results"
    )
    parser.add_argument(
        "--oversample",
        type=int,
        default=10,
        help="Oversampling factor for decoy generation"
    )
    parser.add_argument(
        "--n-proc",
        type=int,
        default=8,
        help="Number of parallel processes"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10000,
        help="Batch size for processing database molecules"
    )
    parser.add_argument(
        "--database-fraction",
        type=float,
        metavar="FRACTION",
        help="Fraction of database to use (0.0-1.0, overrides config)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    parser.add_argument(
        "--no-similarity-filter",
        action="store_true",
        help="Disable similarity filtering after each bundle (uses Max_tc from config)"
    )
    parser.add_argument(
        "--no-report",
        action="store_true",
        help="Disable generation of property comparison report (PDF)"
    )
    
    args = parser.parse_args()
    
    # Validate paths
    if not validate_paths(args):
        sys.exit(1)
    
    # Setup configuration
    if args.verbose:
        print("Setting up configuration...")
        setup_start = time.time()
    else:
        print("Setting up configuration...")
    
    try:
        # Create temporary config file for the generator
        config = setup_config(args.config)
        
        # Apply command line overrides
        if args.database_fraction is not None:
            if not (0.0 <= args.database_fraction <= 1.0):
                print("Error: database_fraction must be between 0.0 and 1.0")
                return
            config.database_fraction = args.database_fraction
        
        temp_config_file = "temp_config.json"
        config.save_to_file(temp_config_file)
        
        # Initialize generator
        if args.verbose:
            setup_time = time.time() - setup_start
            print(f"\nInitializing decoy generator... (config setup: {setup_time:.2f}s)")
            init_start = time.time()
        else:
            print("\nInitializing decoy generator...")
        
        generator = DecoysGenerator(
            enamine_std_path=args.enamine_std,
            prop_bundle_dir=args.prop_bundle_dir,
            smi_bundle_dir=args.smi_bundle_dir,
            prop_prefix=args.prop_prefix,
            smi_prefix=args.smi_prefix,
            config_file=temp_config_file,
            ligand_smi_file=args.ligands,
            oversample=args.oversample,
            n_proc=args.n_proc,
            batch_size=args.batch_size,
            verbose=args.verbose,
            apply_similarity_filter=not args.no_similarity_filter
        )
        
        # Don't clean up temporary config file yet - parallel generator needs it
        # It will be cleaned up at the end of processing
        temp_config_to_cleanup = temp_config_file
        
        if args.verbose:
            init_time = time.time() - init_start
            print(f"Loaded {len(generator.actives)} target molecules (initialization: {init_time:.2f}s)")
            for name, mol in generator.actives.items():
                print(f"  {name}: {mol.smiles} (charge: {mol.charge})")
        else:
            print(f"Loaded {len(generator.actives)} target molecules")
        
        # Generate decoys
        print(f"\nGenerating decoys...")
        if args.verbose:
            print(f"  Using {args.n_proc} processes")
            print(f"  Batch size: {args.batch_size}")
            print(f"  Similarity filtering: {'enabled' if not args.no_similarity_filter else 'disabled'}")
            if not args.no_similarity_filter:
                print(f"  Similarity cutoff: {config.max_tc} (from config)")
            print(f"  Property report: {'enabled' if not args.no_report else 'disabled'}")
            if args.database_fraction is not None:
                print(f"  Database fraction: {args.database_fraction:.2%}")
            elif config.database_fraction < 1.0:
                print(f"  Database fraction: {config.database_fraction:.2%} (from config)")
            generation_start = time.time()
        else:
            print(f"  Using {args.n_proc} processes")
            print(f"  Batch size: {args.batch_size}")
        
        results = generator.generate_decoys()
        
        # Finalize decoys
        if args.verbose:
            generation_time = time.time() - generation_start
            print(f"\nFinalizing decoy selection... (generation: {generation_time:.2f}s)")
            finalize_start = time.time()
        else:
            print("\nFinalizing decoy selection...")
        
        results = generator.finalize_decoys()
        
        # Print summary
        summary = generator.get_summary()
        
        if args.verbose:
            finalize_time = time.time() - finalize_start
            total_time = time.time() - setup_start
            print(f"\nDecoy Generation Complete! (finalization: {finalize_time:.2f}s, total: {total_time:.2f}s)")
        else:
            print(f"\nDecoy Generation Complete!")
        
        print(f"  Target molecules: {summary['num_target_molecules']}")
        print(f"  Total decoys: {summary['total_decoys_generated']}")
        print(f"  Average decoys per molecule: {summary['average_decoys_per_molecule']:.1f}")
        
        if args.verbose:
            print("\nPer-molecule summary:")
            for name, mol_summary in summary['molecule_summaries'].items():
                min_score = mol_summary.get('min_score')
                max_score = mol_summary.get('max_score')
                
                if min_score is not None and max_score is not None:
                    score_range = f"(score range: {min_score:.3f} - {max_score:.3f})"
                else:
                    score_range = "(no decoys found)"
                
                print(f"  {name}: {mol_summary['num_decoys']} decoys {score_range}")
        
        # Save results
        if args.verbose:
            print(f"\nSaving results to {args.output_dir}...")
            save_start = time.time()
        else:
            print(f"\nSaving results to {args.output_dir}...")
        
        save_results(results, args.output_dir)
        
        # Generate property comparison report
        if not args.no_report:
            if args.verbose:
                print("Generating property comparison report...")
                report_start = time.time()
            
            generate_property_report(results, args.output_dir)
            
            if args.verbose:
                save_time = time.time() - save_start
                report_time = time.time() - report_start
                total_save_time = time.time() - save_start
                print(f"Results saved (CSV files: {save_time - report_time:.2f}s, report: {report_time:.2f}s, total: {total_save_time:.2f}s)")
            else:
                print(f"Results and report saved to {args.output_dir}")
        else:
            if args.verbose:
                save_time = time.time() - save_start
                print(f"Results saved (saving: {save_time:.2f}s)")
            else:
                print(f"Results saved to {args.output_dir}")
        
        print("Done!")
        
        # Clean up temporary config file
        try:
            if 'temp_config_to_cleanup' in locals():
                os.remove(temp_config_to_cleanup)
        except:
            pass  # Ignore cleanup errors
        
    except Exception as e:
        # Clean up temporary config file on error too
        try:
            if 'temp_config_to_cleanup' in locals():
                os.remove(temp_config_to_cleanup)
        except:
            pass
            
        print(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()