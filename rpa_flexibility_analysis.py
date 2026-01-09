#!/usr/bin/env python3
"""
Structural Flexibility Analysis: Experimental B-factors vs. AlphaFold pLDDT
===========================================================================

This script compares the local flexibility of RPA70 observed in X-ray
crystallography (B-factors) against the predicted confidence by AlphaFold
(pLDDT scores).

Author: [Your Name]
Date: 2025-01-09
"""

import os
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from pathlib import Path

# Set style for publication-quality figures
plt.style.use('seaborn-v0_8-darkgrid')

# ============================================================================
# CONSTANTS
# ============================================================================

# File URLs
PDB_XRAY_URL = "https://files.rcsb.org/download/1JMC.pdb"
PDB_ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/files/AF-P27694-F1-model_v6.pdb"

# Local file names
XRAY_FILE = "1JMC.pdb"
ALPHAFOLD_FILE = "AF-P27694-F1-model_v6.pdb"
PYMOL_SCRIPT = "color_by_flexibility.pml"

# Chain ID for the experimental structure (Chain A of 1JMC)
TARGET_CHAIN = "A"


# ============================================================================
# STEP 1: DATA FETCHING
# ============================================================================

def download_pdb_file(url: str, filename: str, force_download: bool = False) -> str:
    """
    Download a PDB file from the given URL if it doesn't exist locally.

    Parameters
    ----------
    url : str
        URL to download the PDB file from
    filename : str
        Local filename to save the file
    force_download : bool, optional
        Force re-download even if file exists (default: False)

    Returns
    -------
    str
        Path to the downloaded file
    """
    if os.path.exists(filename) and not force_download:
        print(f"  [SKIP] {filename} already exists locally.")
        return filename

    print(f"  [DOWNLOAD] Fetching {filename} from {url}...")
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        with open(filename, 'wb') as f:
            f.write(response.content)

        print(f"  [SUCCESS] Downloaded {filename} ({len(response.content)} bytes)")
        return filename
    except requests.RequestException as e:
        raise RuntimeError(f"Failed to download {filename}: {e}")


def fetch_all_data(force_download: bool = False) -> tuple:
    """
    Download both experimental and AlphaFold PDB files.

    Parameters
    ----------
    force_download : bool, optional
        Force re-download even if files exist

    Returns
    -------
    tuple
        (xray_file_path, alphafold_file_path)
    """
    print("\n" + "="*60)
    print("STEP 1: FETCHING STRUCTURAL DATA")
    print("="*60)

    xray_path = download_pdb_file(PDB_XRAY_URL, XRAY_FILE, force_download)
    alphafold_path = download_pdb_file(PDB_ALPHAFOLD_URL, ALPHAFOLD_FILE, force_download)

    print("\n  [INFO] All files downloaded successfully.")
    return xray_path, alphafold_path


# ============================================================================
# STEP 2: PARSING & ALIGNMENT
# ============================================================================

def parse_xray_bfactors(pdb_file: str, chain_id: str) -> pd.DataFrame:
    """
    Parse X-ray structure and extract average C-alpha B-factors per residue.

    Parameters
    ----------
    pdb_file : str
        Path to the X-ray PDB file
    chain_id : str
        Chain ID to extract (e.g., "A")

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: [residue_number, bfactor]
    """
    print(f"\n  [PARSE] Reading X-ray structure from {pdb_file}...")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("xray", pdb_file)

    # Get the target chain
    chain = structure[0][chain_id]

    # Extract residue numbers and average B-factors
    data = []
    for residue in chain:
        if residue.id[0] != ' ':  # Skip heteroatoms/water
            continue

        res_num = residue.id[1]

        # Get C-alpha atom B-factor if available
        if 'CA' in residue:
            bfactor = residue['CA'].get_bfactor()
            data.append({'residue_number': res_num, 'bfactor': bfactor})

    df = pd.DataFrame(data)
    df = df.drop_duplicates(subset=['residue_number']).sort_values('residue_number')

    print(f"  [INFO] Extracted {len(df)} residues from Chain {chain_id}")
    print(f"  [INFO] Residue range: {df['residue_number'].min()} - {df['residue_number'].max()}")
    print(f"  [INFO] B-factor range: {df['bfactor'].min():.2f} - {df['bfactor'].max():.2f} A^2")

    return df


def parse_alphafold_plddt(pdb_file: str) -> pd.DataFrame:
    """
    Parse AlphaFold structure and extract pLDDT scores per residue.
    In AlphaFold PDB files, pLDDT is stored in the B-factor column.

    Parameters
    ----------
    pdb_file : str
        Path to the AlphaFold PDB file

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: [residue_number, plddt]
    """
    print(f"\n  [PARSE] Reading AlphaFold structure from {pdb_file}...")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("alphafold", pdb_file)

    # AlphaFold typically uses Chain A
    chain = structure[0]['A']

    # Extract residue numbers and pLDDT (stored in B-factor column)
    data = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != ' ':  # Skip heteroatoms
                    continue

                res_num = residue.id[1]

                # Get pLDDT from C-alpha atom
                if 'CA' in residue:
                    plddt = residue['CA'].get_bfactor()
                    data.append({'residue_number': res_num, 'plddt': plddt})

    df = pd.DataFrame(data)
    df = df.drop_duplicates(subset=['residue_number']).sort_values('residue_number')

    print(f"  [INFO] Extracted {len(df)} residues from AlphaFold prediction")
    print(f"  [INFO] Residue range: {df['residue_number'].min()} - {df['residue_number'].max()}")
    print(f"  [INFO] pLDDT range: {df['plddt'].min():.2f} - {df['plddt'].max():.2f}")

    return df


def align_residues(xray_df: pd.DataFrame, af_df: pd.DataFrame) -> pd.DataFrame:
    """
    Align the two dataframes to only include overlapping residues.

    Parameters
    ----------
    xray_df : pd.DataFrame
        DataFrame with X-ray B-factors
    af_df : pd.DataFrame
        DataFrame with AlphaFold pLDDT scores

    Returns
    -------
    pd.DataFrame
        Aligned DataFrame with columns:
        [residue_number, bfactor, plddt, predicted_error]
    """
    print("\n  [ALIGN] Aligning residues between structures...")

    # Merge on residue number (inner join = only overlapping residues)
    merged = pd.merge(
        xray_df,
        af_df,
        on='residue_number',
        how='inner'
    ).sort_values('residue_number')

    # Calculate predicted error from pLDDT
    merged['predicted_error'] = 100 - merged['plddt']

    print(f"  [INFO] Found {len(merged)} overlapping residues")
    print(f"  [INFO] Aligned residue range: {merged['residue_number'].min()} - {merged['residue_number'].max()}")

    return merged


# ============================================================================
# STEP 3: ANALYSIS
# ============================================================================

def min_max_normalize(series: pd.Series) -> pd.Series:
    """
    Normalize a series to 0-1 range using Min-Max normalization.

    Parameters
    ----------
    series : pd.Series
        Input series to normalize

    Returns
    -------
    pd.Series
        Normalized series (0-1 range)
    """
    min_val = series.min()
    max_val = series.max()

    if max_val == min_val:
        return pd.Series([0.5] * len(series), index=series.index)

    return (series - min_val) / (max_val - min_val)


def analyze_flexibility(merged_df: pd.DataFrame) -> pd.DataFrame:
    """
    Perform normalization and correlation analysis.

    Parameters
    ----------
    merged_df : pd.DataFrame
        Aligned DataFrame with bfactor, plddt, and predicted_error

    Returns
    -------
    pd.DataFrame
        DataFrame with added normalized columns and correlation info
    """
    print("\n" + "="*60)
    print("STEP 3: ANALYSIS")
    print("="*60)

    # Normalize metrics
    merged_df['bfactor_norm'] = min_max_normalize(merged_df['bfactor'])
    merged_df['predicted_error_norm'] = min_max_normalize(merged_df['predicted_error'])

    # Calculate Pearson correlation
    correlation = merged_df['bfactor_norm'].corr(merged_df['predicted_error_norm'])

    print(f"\n  [NORMALIZE] Applied Min-Max normalization (0-1 scale)")
    print(f"  [CORRELATION] Pearson r = {correlation:.4f}")
    print(f"  [INFO] Positive r = higher B-factor correlates with lower pLDDT (higher flexibility)")

    return merged_df, correlation


# ============================================================================
# STEP 4: VISUALIZATION
# ============================================================================

def create_dual_axis_plot(merged_df: pd.DataFrame, correlation: float,
                         output_file: str = "flexibility_comparison.png"):
    """
    Create a dual-axis plot comparing B-factors and pLDDT scores.

    Parameters
    ----------
    merged_df : pd.DataFrame
        DataFrame with normalized metrics
    correlation : float
        Pearson correlation coefficient
    output_file : str, optional
        Output filename for the figure
    """
    print("\n" + "="*60)
    print("STEP 4: VISUALIZATION")
    print("="*60)

    fig, ax1 = plt.subplots(figsize=(14, 6))

    residues = merged_df['residue_number'].values

    # Left Y-axis: Normalized B-factor (Experimental Flexibility)
    color1 = '#1f77b4'  # Blue
    ax1.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Normalized B-Factor\n(Experimental Flexibility)',
                   color=color1, fontsize=12, fontweight='bold')
    line1 = ax1.plot(residues, merged_df['bfactor_norm'],
                     color=color1, linewidth=2, label='B-factor (X-ray)', alpha=0.8)
    ax1.tick_params(axis='y', labelcolor=color1, labelsize=11)
    ax1.tick_params(axis='x', labelsize=11)
    ax1.set_ylim(0, 1)
    ax1.grid(True, alpha=0.3)

    # Right Y-axis: Normalized pLDDT (Predicted Confidence) - INVERTED
    color2 = '#ff7f0e'  # Orange
    ax2 = ax1.twinx()
    ax2.set_ylabel('Normalized AlphaFold Confidence\n(pLDDT, Inverted)',
                   color=color2, fontsize=12, fontweight='bold')

    # Invert pLDDT so high confidence = low value (matches low B-factor)
    inverted_plddt = 1 - merged_df['predicted_error_norm']
    line2 = ax2.plot(residues, inverted_plddt,
                     color=color2, linewidth=2, label='pLDDT (AlphaFold)', alpha=0.8)
    ax2.tick_params(axis='y', labelcolor=color2, labelsize=11)
    ax2.set_ylim(0, 1)

    # Highlight the 1JMC coverage region
    ax1.axvspan(merged_df['residue_number'].min(), merged_df['residue_number'].max(),
                alpha=0.15, color='gray', label='1JMC coverage')

    # Title with correlation coefficient
    title = f'RPA70 Flexibility Analysis: B-factor vs AlphaFold pLDDT\n'
    title += f'PDB: 1JMC (Chain A) | UniProt: P27694 | '
    title += f'Correlation: r = {correlation:.4f}'
    plt.title(title, fontsize=14, fontweight='bold', pad=20)

    # Combined legend - moved to upper right to avoid overlapping with text box
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper right', framealpha=0.9, fontsize=10)

    # Add text annotation
    info_text = f"n = {len(merged_df)} residues\n"
    info_text += f"X-ray range: {merged_df['residue_number'].min()}-{merged_df['residue_number'].max()}"
    ax1.text(0.02, 0.98, info_text, transform=ax1.transAxes,
             fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n  [SAVE] Figure saved to {output_file}")

    return fig


def create_scatter_plot(merged_df: pd.DataFrame, correlation: float,
                       output_file: str = "flexibility_scatter.png"):
    """
    Create a scatter plot showing the relationship between B-factors and pLDDT.

    Parameters
    ----------
    merged_df : pd.DataFrame
        DataFrame with normalized metrics
    correlation : float
        Pearson correlation coefficient
    output_file : str, optional
        Output filename for the figure
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    x = merged_df['bfactor_norm']
    y = merged_df['predicted_error_norm']

    # Scatter plot with color mapping by residue position
    scatter = ax.scatter(x, y, c=merged_df['residue_number'],
                        cmap='viridis', s=50, alpha=0.6, edgecolors='black', linewidth=0.5)

    # Add regression line
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    ax.plot(x, p(x), "r--", alpha=0.8, linewidth=2, label=f'Linear fit: y = {z[0]:.2f}x + {z[1]:.2f}')

    ax.set_xlabel('Normalized B-Factor (X-ray Flexibility)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Normalized Predicted Error\n(100 - pLDDT)', fontsize=12, fontweight='bold')
    ax.set_title(f'RPA70: B-factor vs AlphaFold Predicted Error\nCorrelation: r = {correlation:.4f}',
                fontsize=14, fontweight='bold')

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Residue Number', fontsize=11)

    ax.legend(fontsize=10, loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  [SAVE] Scatter plot saved to {output_file}")

    return fig


# ============================================================================
# STEP 5: PYMOL SCRIPT GENERATION
# ============================================================================

def generate_pymol_script(pdb_file: str, output_file: str = PYMOL_SCRIPT):
    """
    Generate a PyMOL script to color residues by B-factor.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file to load
    output_file : str, optional
        Output filename for the PyMOL script
    """
    print("\n" + "="*60)
    print("STEP 5: PYMOL SCRIPT GENERATION")
    print("="*60)

    script_content = f"""# PyMOL Script: Color RPA70 by Flexibility (B-factor)
# ===========================================================
# This script loads the 1JMC structure and colors residues
# from blue (rigid/low B-factor) to red (flexible/high B-factor)
#
# Usage: pymol {output_file}
# ===========================================================

# Load the structure
load {pdb_file}, rpa70

# Hide everything except protein
hide everything
show cartoon, rpa70

# Color by B-factor (spectrum coloring)
# Blue = low B-factor (rigid), Red = high B-factor (flexible)
spectrum b, blue_red, minimum=0, maximum=100, rpa70

# Set background
bg_color white

# Set representation
set cartoon_smooth_loops, 0
set cartoon_cylindrical_helices, 1
set cartoon_cylindrical_sheets, 1

# Show as cartoon with transparent surface
show surface, rpa70
set transparency, 0.5, rpa70

# Zoom to fit
zoom rpa70

# Add a color key (B-factor scale)
# Note: This creates a simple representation
# Blue  = Rigid (low B-factor, ~0-30 Å²)
# White = Medium (medium B-factor, ~30-60 Å²)
# Red   = Flexible (high B-factor, ~60+ Å²)

# Print statistics
print "=========================================="
print "RPA70 Flexibility Visualization"
print "=========================================="
print "Color scheme: B-factor (flexibility)"
print "  Blue  = Rigid (low B-factor)"
print "  White = Medium flexibility"
print "  Red   = Flexible (high B-factor)"
print "=========================================="
"""

    with open(output_file, 'w') as f:
        f.write(script_content)

    print(f"  [SAVE] PyMOL script saved to {output_file}")
    print(f"  [INFO] Run with: pymol {output_file}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function."""
    print("\n" + "="*70)
    print(" STRUCTURAL FLEXIBILITY ANALYSIS")
    print(" B-factors (X-ray) vs pLDDT (AlphaFold)")
    print(" Target: RPA70 (Human RPA1)")
    print("="*70)

    try:
        # Step 1: Fetch data
        xray_file, af_file = fetch_all_data()

        # Step 2: Parse structures
        print("\n" + "="*60)
        print("STEP 2: PARSING & ALIGNMENT")
        print("="*60)
        xray_df = parse_xray_bfactors(xray_file, TARGET_CHAIN)
        af_df = parse_alphafold_plddt(af_file)
        merged_df = align_residues(xray_df, af_df)

        # Step 3: Analyze
        merged_df, correlation = analyze_flexibility(merged_df)

        # Step 4: Visualize
        create_dual_axis_plot(merged_df, correlation)
        create_scatter_plot(merged_df, correlation)

        # Step 5: Generate PyMOL script
        generate_pymol_script(xray_file)

        # Summary
        print("\n" + "="*70)
        print(" ANALYSIS COMPLETE")
        print("="*70)
        print(f"\nGenerated files:")
        print(f"  - flexibility_comparison.png (Dual-axis plot)")
        print(f"  - flexibility_scatter.png (Scatter plot)")
        print(f"  - {PYMOL_SCRIPT} (PyMOL script)")
        print(f"\nKey findings:")
        print(f"  - Overlapping residues: {len(merged_df)}")
        print(f"  - Correlation (r): {correlation:.4f}")
        print(f"  - Residue range: {merged_df['residue_number'].min()}-{merged_df['residue_number'].max()}")

        # Save aligned data to CSV
        csv_file = "aligned_flexibility_data.csv"
        merged_df.to_csv(csv_file, index=False)
        print(f"\n  [SAVE] Aligned data saved to {csv_file}")

        print("\n" + "="*70)

    except Exception as e:
        print(f"\n[ERROR] {str(e)}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
