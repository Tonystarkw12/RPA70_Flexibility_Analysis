# RPA70 Structural Flexibility Analysis

## Overview
This script compares experimental flexibility (B-factors from X-ray crystallography) with predicted confidence (pLDDT from AlphaFold) for the human RPA70 protein (DNA-binding domains A and B).

## Quick Start

### 1. Install Dependencies
```bash
pip install -r requirements.txt
```

### 2. Run the Analysis
```bash
python rpa_flexibility_analysis.py
```

## Input Data

| Source | Type | Identifier | Description |
|--------|------|------------|-------------|
| X-ray | PDB | 1JMC | Human RPA70 DNA-binding domains A + B with ssDNA |
| AlphaFold | PDB | AF-P27694-F1 | Full-length human RPA1/RPA70 prediction |

## Output Files

| File | Description |
|------|-------------|
| `flexibility_comparison.png` | Dual-axis plot comparing normalized B-factors and pLDDT |
| `flexibility_scatter.png` | Scatter plot with correlation analysis |
| `color_by_flexibility.pml` | PyMOL script to color structure by B-factor |
| `aligned_flexibility_data.csv` | Raw aligned data for further analysis |

## Key Features

1. **Automatic Data Fetching**: Downloads PDB files directly from RCSB and AlphaFold
2. **Residue Alignment**: Only compares overlapping residues (1JMC covers ~181-422)
3. **Min-Max Normalization**: Makes metrics comparable on 0-1 scale
4. **Pearson Correlation**: Quantifies relationship between experimental and predicted flexibility
5. **Dual-Axis Visualization**: Inverted pLDDT axis for intuitive comparison

## Interpretation

- **High B-factor + Low pLDDT** = Experimentally flexible, predicted disordered
- **Low B-factor + High pLDDT** = Experimentally rigid, predicted confident
- **Positive correlation** indicates AlphaFold pLDDT captures experimental flexibility patterns
