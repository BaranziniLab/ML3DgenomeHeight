# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Project Overview

This repository contains code for "Machine learning reveals 3D regulatory mechanisms for height-associated haplotypes" - a research project that applies machine learning to analyze 3D genome structure in relation to genetic variants associated with human height.

## High-Level Architecture

The project is organized into three main components:

### 1. LDSC Enrichment Analysis (`LDSC Enrichment/`)
- **Purpose**: Performs linkage disequilibrium score regression enrichment analysis
- **Key Components**:
  - `make_ldsc_data/`: Scripts for generating LDSC-compatible annotations and GWAS data
  - `make_annot.R`: Creates annotation files from variant data
  - `munge_gwas.R`: Processes GWAS summary statistics for multiple populations
  - `gwas_raw/`: Contains GWAS summary statistics for different populations (AFR, AMR, EAS, EUR, SAS)
  - `akita_variants/`: Variant data from Akita model predictions

### 2. Hi-C Data Visualization (`Visualize HiC/`)
- **Purpose**: Creates publication-ready visualizations of Hi-C contact maps
- **Key Components**:
  - `visualization_hic.ipynb`: Jupyter notebook with specialized functions for Hi-C visualization
  - Includes diamond-rotated contact maps with genomic annotations
  - Functions for loading predictions, adding genomic tracks, and comparing contact maps

### 3. Project Documentation and Results
- **Figures/**: Contains publication figures and supplementary materials
- **ASHG Abstract/** and **ASHG Presentation/**: Conference materials
- **README.md**: Basic project description

## Common Development Commands

### LDSC Analysis
```bash
# Generate LDSC annotations (run from LDSC Enrichment/make_ldsc_data/)
Rscript make_annot.R

# Process GWAS data for multiple populations
Rscript munge_gwas.R
```

### Hi-C Visualization
```bash
# Launch Jupyter notebook for Hi-C visualization
jupyter notebook "Visualize HiC/visualization_hic.ipynb"
```

### Version Control
```bash
# Auto-sync scripts are provided for quick commits
./auto_sync.sh    # On Unix systems
auto_sync.bat     # On Windows
```

## Key Analysis Functions

### Hi-C Visualization (`visualization_hic.ipynb`)
- `from_upper_triu()`: Converts upper triangular vector to symmetric matrix
- `load_individual_map()`: Loads Hi-C contact maps from prediction files
- `annotations()`: Adds genomic annotations (genes, CTCF, conservation)
- `map_comparison_no_delta_with_annotations()`: Creates side-by-side contact map comparisons with genomic tracks

### Data Processing
- **Variant Data**: Located in `akita_variants/variants_in_divergent_windows_10AF.txt`
  - Format: chr, position, ref_allele, alt_allele, window_position
- **GWAS Data**: Population-specific files in `gwas_raw/`
  - Format: Standard GWAS summary statistics with hg38 coordinates

## Architecture Notes

### Data Flow
1. **Input**: Height GWAS variants and Hi-C predictions from Akita model
2. **Processing**: LDSC enrichment analysis and Hi-C visualization
3. **Output**: Publication figures and enrichment statistics

### Key Dependencies
- **R**: For LDSC data processing and statistical analysis
- **Python/Jupyter**: For Hi-C visualization with matplotlib, numpy, pandas
- **External Tools**: LDSC software suite (referenced but not included)

### File Organization
- Raw data files are organized by analysis type and population
- Visualization code is self-contained in Jupyter notebooks
- Auto-sync scripts maintain Git synchronization

## Development Notes

### Running Hi-C Visualizations
The visualization notebook expects specific file paths for annotation data:
- Gene annotations: `grch38_gene_annotations.bed`
- Conservation data: `phastConsElements100way_hg38.bed`
- CTCF binding sites: `ctcf_full_merged_hg38.bed`

### LDSC Data Processing
- GWAS files are expected in standard format with hg38 coordinates
- Variant annotations are derived from Akita model predictions
- Population-specific analysis is supported for AFR, AMR, EAS, EUR, and SAS populations

### Git Workflow
The project includes automated sync scripts that perform standard Git operations:
1. Stage all changes (`git add .`)
2. Commit with automatic message
3. Pull with unrelated histories flag
4. Push to main branch

This setup facilitates rapid iteration during active research phases.
