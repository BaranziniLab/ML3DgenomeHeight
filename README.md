# Machine Learning Reveals 3D Regulatory Mechanisms for Height-Associated Haplotypes

This repository contains the complete computational pipeline and analysis code for studying how genetic variants associated with human height affect 3D chromatin organization and gene regulation. The project combines machine learning approaches with genomics data to understand the regulatory mechanisms underlying height-associated genetic variants.

## Abstract

Genome-wide association studies (GWAS) have identified thousands of genetic variants associated with human height, but the mechanisms by which these variants affect gene regulation through 3D chromatin structure remain poorly understood. This study applies machine learning models to predict how height-associated variants alter 3D genome organization and gene expression patterns.

## Project Structure

### 📊 LDSC Enrichment Analysis (`LDSC Enrichment/`)

Linkage Disequilibrium Score Regression (LDSC) analysis to test for enrichment of height-associated variants in specific genomic regions.

#### Key Scripts:
- **`make_ldsc_data/make_annot.R`**: Generates LDSC-compatible annotation files from variant data
- **`make_ldsc_data/munge_gwas.R`**: Processes and formats GWAS summary statistics for multiple ancestries

#### Data Resources:
- **`akita_variants/variants_in_divergent_windows_10AF.txt`**: Height-associated variants identified from Akita model predictions
  - Format: chromosome, position, reference allele, alternative allele, window position
  - Contains variants showing divergent effects on chromatin organization

- **`gwas_raw/`**: Population-specific GWAS summary statistics for human height
  - `afr_hg38.txt`: African ancestry populations
  - `amr_hg38.txt`: Admixed American populations  
  - `eas_hg38.txt`: East Asian populations
  - `eur_hg38.txt`: European ancestry populations
  - `sas_hg38.txt`: South Asian populations
  - `all_hg38.txt`: Combined multi-ancestry analysis
  - All files use hg38 genomic coordinates with standard GWAS format (SNPID, RSID, CHR, POS, etc.)

### 🧬 Hi-C Data Visualization (`Visualize HiC/`)

Specialized tools for creating publication-ready visualizations of chromatin contact maps with genomic annotations.

#### Main Notebook:
- **`visualization_hic.ipynb`**: Comprehensive Jupyter notebook containing:
  - Custom functions for Hi-C data processing and visualization
  - Diamond-rotated contact map plotting
  - Integration of multiple genomic annotation tracks
  - Side-by-side comparison tools for different conditions

#### Key Functions:
- `from_upper_triu()`: Converts upper triangular vector format to symmetric contact matrices
- `load_individual_map()`: Loads and processes Hi-C predictions from model outputs
- `annotations()`: Adds genomic annotation tracks (genes, CTCF sites, conservation)
- `map_comparison_no_delta_with_annotations()`: Creates comparative visualizations with genomic context

#### Required Annotation Files (not included):
- `grch38_gene_annotations.bed`: Gene annotation data
- `phastConsElements100way_hg38.bed`: Evolutionary conservation tracks
- `ctcf_full_merged_hg38.bed`: CTCF binding sites

### 📄 Documentation and Results

- **`Figures/`**: Publication-quality figures and supplementary materials
  - `figure1.png`: Main results overview
  - `figure2.png`: LDSC enrichment analysis results
  - `figure3.png`: Hi-C visualization examples
  - `figure4.png`: Mechanistic model
  - `figure s1.png`: Supplementary analysis

- **`ASHG Abstract/`**: Conference materials
  - Abstract submission for the American Society of Human Genetics meeting

- **`ASHG Presentation/`**: Conference presentation materials
  - PowerPoint presentation with key findings and visualizations

- **`Height 3D Genome MS.pdf`**: Full manuscript (if available)

## Getting Started

### Prerequisites

- **R** (≥4.0) with packages for statistical analysis
- **Python** (≥3.8) with Jupyter notebook support
- Required Python packages: `matplotlib`, `numpy`, `pandas`, `scipy`
- Access to LDSC software suite (for enrichment analysis)

### Running the Analysis

#### 1. LDSC Enrichment Analysis
```bash
cd "LDSC Enrichment/make_ldsc_data/"

# Generate annotation files
Rscript make_annot.R

# Process GWAS summary statistics
Rscript munge_gwas.R
```

#### 2. Hi-C Visualizations
```bash
# Launch Jupyter notebook
jupyter notebook "Visualize HiC/visualization_hic.ipynb"

# Follow the notebook cells to:
# - Load Hi-C prediction data
# - Add genomic annotations
# - Generate publication figures
```

#### 3. Quick Git Synchronization
```bash
# Use provided auto-sync scripts for version control
./auto_sync.sh      # Unix/Mac
# OR
auto_sync.bat       # Windows
```

## Data Sources

- **Height GWAS**: Multi-ancestry genome-wide association study data
- **3D Genome Predictions**: Machine learning models (Akita) for chromatin organization
- **Genomic Annotations**: Gene locations, regulatory elements, conservation scores
- **Population Genetics**: Linkage disequilibrium patterns across human populations

## Key Findings

This analysis reveals:
1. Height-associated variants are enriched in regions with altered 3D chromatin structure
2. Population-specific differences in regulatory mechanisms
3. Novel regulatory targets for height-associated loci
4. Mechanistic insights into how genetic variants affect gene regulation through chromatin organization

## Citation

*Manuscript in preparation*

Preprint available soon

## Contact

For questions about the analysis or code, please open an issue in this repository.

## License

*MIT License*
