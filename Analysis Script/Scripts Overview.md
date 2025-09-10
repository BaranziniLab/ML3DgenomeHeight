# Height 3D-Genome Analysis Project

This repository contains R scripts for analyzing the relationship between human height GWAS variants and 3D genome architecture using the Akita deep learning model. The project investigates how genetic variants affecting human height alter chromatin interactions across multiple ancestries.

## Overview
This comprehensive analytical pipeline integrates genome-wide association study (GWAS) data with 3D chromatin interaction predictions to identify regulatory mechanisms underlying height variation. The analysis combines haplotype-based approaches with machine learning predictions to quantify how genetic variants alter 3D genome organization.

## Directory Structure and R Scripts

### 1. GIANT GWAS Data Processing (`giant_gwas/`)

Scripts for processing GIANT consortium height GWAS data:

- **`liftover.R`** - Converts genomic coordinates from hg19 to hg38 using liftOver chain files. Processes GWAS summary statistics for all ancestries (AFR, EUR, AMR, EAS, SAS).

- **`make_plink_files.R`** - Converts GWAS summary statistics to PLINK format for clumping analysis. Creates files with SNP ID, chromosome, position, and P-value columns.

- **`generate_hapmap_clump_runners.R`** - Generates SLURM job scripts for running PLINK clumping using HapMap reference panel. Creates ancestry-specific batch jobs with varying sample sizes.

- **`generate_tgp_clump_runner.R`** - Similar to HapMap version but uses 1000 Genomes Project (TGP) reference panel for clumping analysis.

- **`generate_topmed_clump_runner.R`** - Generates clumping scripts using TopMed reference panel for more comprehensive variant coverage.

### 2. LDSC Enrichment Analysis (`ldsc_enrichment/`)

Scripts for linkage disequilibrium score regression enrichment analysis:

#### `make_ldsc_data/`
- **`munge_gwas.R`** - Prepares GWAS summary statistics for LDSC analysis by formatting and filtering to HapMap3 SNPs.

- **`make_annot.R`** - Creates annotation files for LDSC by mapping Akita-identified 3D-genome-altering variants to genomic positions. Performs coordinate conversion from hg38 to hg19.

#### `regionr_enrichment/`
- **`calc_enrichment.R`** - Performs permutation testing using regioneR to test for enrichment of 3D-genome-altering variants in height GWAS loci.

### 3. Akita Model Analysis (`akita_rerun/`)

Scripts for analyzing Akita model predictions and 3D genome effects:

- **`cleanup.R`** - Main analysis script that processes divergence matrices, calculates maximum 3D genome divergence effects, and creates visualizations comparing haplotypic vs single-SNP effects.

- **`calc_pvalue.R`** - Statistical analysis of the relationship between effect sizes and 3D genome divergence. Creates scatter plots and regression analyses.

- **`make_table1.R`** - Generates summary table of top variants with significant 3D genome effects (divergence ≥ 0.01).

- **`make_scale_bar.R`** - Creates color scale bar for Hi-C contact map visualizations.

- **`snp_contribution.R`** - Analyzes individual SNP contributions to 3D genome effects at multi-SNP loci. Identifies cases where effects are distributed across multiple variants.

### 4. LCOR Locus Analysis (`lcor_chip/`)

- **`make_barplot.R`** - Analyzes ChIP-seq binding patterns at the LCOR locus (chr10:96789820), creating barplots of transcription factor binding frequencies.

### 5. GIANT GWAS Loci Analysis (`giant_gwas_loci/`)

Scripts for analyzing height-associated loci:

- **`integrate_results.R`** - Integrates clumping results with association statistics, applies FDR correction, and identifies replicated associations across ancestries.

- **`sanity_check.R`** - Quality control analysis comparing results between TopMed and published GWAS (Yengo et al.) for validation.

#### `clumped_result/`
- **`make_map.R`** - Creates genomic interval maps from clumped results for downstream analysis, defining ±1Mb windows around lead SNPs.

### 6. 1000 Genomes Project Analysis (`tgp/`)

Comprehensive haplotype analysis using 1000 Genomes data:

- **`make_haplotype.R`** - Main script for extracting haplotypes from 1000 Genomes VCF files. Processes phased genotype data and creates haplotype frequency tables across ancestries.

- **`make_sequence.R`** - Converts haplotype data to DNA sequences for Akita model input. Handles sequence padding and formatting for 1Mb fixed-length sequences.

- **`run1.R` & `run2.R`** - Batch processing scripts for generating haplotype sequences (run1: variants 1-5335, run2: variants 5336-10671).

- **`test_run.R`** - Test script for validating the sequence generation pipeline.

- **`investigate_lcor.R`** - Detailed analysis of the LCOR locus, including in silico mutagenesis and PCA of chromatin contact predictions.

- **`investigate_lcor_pca.R`** - Principal component analysis of LCOR locus 3D genome effects, examining individual SNP contributions to PC space.

- **`investigate_beta.R`** - Extracts GWAS effect sizes for specific loci across all ancestries for downstream analysis.

- **`visualize_prediction.R`** - Creates contact map visualizations and correlation analyses of Akita predictions across different haplotypes.

- **`make_qq.R`** - Generates Q-Q plots for assessing the distribution of maximum 3D genome divergence effects.

- **`analyse_divergence_mat.R`** - Comprehensive analysis of divergence matrices, identifying top divergent regions and generating visualization commands.

- **`insilico_mutagenesis.R`** - Performs systematic in silico mutagenesis by generating all possible haplotype combinations at a locus.

#### `calc_cormat/`
- **`calc_cormat.R`** - Calculates Spearman correlation matrices between Akita predictions for different haplotypes and computes divergence as 1-correlation.

### 7. TopMed Reference Processing (`topmed_ref_subsetting/`)

Scripts for processing TopMed reference panels:

- **`subset_bim.R`** - Subsets PLINK binary files to specific SNP sets (HapMap, 1000 Genomes, TopMed).

#### `HapMap/`
- **`liftover.R`** - Lifts over HapMap3 coordinates from hg18 to hg38.
- **`annotate.R`** - Annotates HapMap SNPs using Ensembl biomaRt for position verification.

### 8. TopMed Phenotype Processing (`topmed_3dgenome/`)

- **`get_pheno.R`** - Extracts and harmonizes phenotype data from multiple TopMed cohorts (ARIC, CARDIA, CHS, HCHS, MESA, WHI, FHS). Prepares phenotype files for association testing with height, lipids, and lung function.

## Data Sources

### GWAS Data
- **Source**: GIANT Consortium Height GWAS (Yengo et al., 2022)
- **Populations**: Multi-ancestry analysis including:
  - European (EUR)
  - African/African American (AFR) 
  - East Asian (EAS)
  - South Asian (SAS)
  - Hispanic/Latino (AMR)
  - Trans-ancestry meta-analysis (ALL)

### Reference Genomes and Annotations
- **Genome builds**: hg19 and hg38 (with liftOver conversion)
- **Population references**: 1000 Genomes Project, TOPMed
- **Functional annotations**: 
  - Gene annotations (GRCh38)
  - PhastCons conservation elements
  - CTCF binding sites
  - HapMap3 SNPs for LD Score Regression

## Analytical Methods

### 1. GWAS Data Processing
**Scripts**: `giant_gwas/download_giant.sh`, `giant_gwas/liftover.R`

- Downloaded multi-ancestry height GWAS summary statistics from GIANT consortium
- Performed genome coordinate liftover from hg19 to hg38 using `rtracklayer`
- Standardized variant identifiers across populations
- Generated population-specific datasets for downstream analysis

### 2. Genetic Association Analysis
**Scripts**: `giant_gwas_loci/height_assoc_runner.sh`, `giant_gwas_loci/integrate_results.R`

- **Software**: PLINK2
- **Covariates**: Age, sex, and first 5 principal components
- **Quality control**: 
  - Minor allele frequency (MAF) ≥ 0.001
  - Covariate variance standardization
- **Statistical approach**: Linear regression with population-specific analysis
- **Multiple testing correction**: False discovery rate (FDR) ≤ 0.05
- **LD clumping**: Population-specific clumping using 1000 Genomes and TOPMed references

### 3. 3D Genome Prediction and Analysis
**Scripts**: `akita_rerun/`, `akita_visualize/`

- **Machine learning model**: Akita deep learning model for 3D chromatin interaction prediction
- **Input**: DNA sequence variants and their impact on chromatin architecture
- **Output**: Hi-C contact map predictions showing 3D genome organization changes
- **Divergence calculation**: Spearman correlation-based divergence metrics to quantify 3D structural changes
- **Visualization**: 
  - Hi-C contact map visualization with rotated diamond plots
  - Comparative analysis between reference and alternate alleles
  - Integration with genomic annotations (genes, conservation, CTCF)

### 4. Enrichment Analysis
**Scripts**: `ldsc_enrichment/`, `ldsc_enrichment/regionr_enrichment/calc_enrichment.R`

#### LD Score Regression (LDSC)
- **Purpose**: Heritability enrichment analysis of 3D genome-altering variants
- **Input**: 
  - Munged GWAS summary statistics
  - Binary annotations for variants causing 3D genome divergence
  - HapMap3 SNP set for LD structure
- **Approach**: Stratified LD Score Regression to test enrichment of height heritability

#### RegionR Permutation Testing
- **Statistical framework**: Permutation-based enrichment testing
- **Null hypothesis**: Random genomic distribution
- **Test regions**: Height-associated loci (±100kb windows)
- **Query regions**: 3D genome-altering variants (divergence > 0.1)
- **Permutations**: 1,000 randomizations with `resampleRegions`
- **Universe**: Autosomal chromosomes 1-22

### 5. Statistical Analysis and Visualization
**Scripts**: `akita_rerun/calc_pvalue.R`, `akita_rerun/make_table1.R`, `akita_rerun/make_barplot.R`

- **Correlation analysis**: Relationship between effect size and 3D genome divergence
- **Statistical tests**: 
  - t-tests comparing haplotypic vs. non-haplotypic variants
  - Linear regression modeling effect size vs. divergence
- **Data filtering**: Variants with divergence > 0.001 for robust analysis
- **Visualization**: 
  - Scatter plots with linear regression fits
  - Bar plots for population-specific results
  - Multi-panel figures with statistical annotations

## Key Methodological Features

### Quality Control
- Population-specific analysis to account for genetic ancestry
- Stringent MAF thresholds and covariate adjustment
- FDR correction for multiple testing
- Systematic coordinate system conversion (hg19→hg38)

### Machine Learning Integration
- Deep learning-based prediction of 3D chromatin structure
- Quantitative assessment of variant impact on genome architecture
- High-resolution contact map analysis (2048bp resolution)

### Multi-Scale Analysis
- Variant-level effects on local chromatin structure
- Locus-level association with height phenotype  
- Genome-wide enrichment patterns
- Population-specific replication

### Statistical Rigor
- Multiple independent statistical approaches
- Permutation-based significance testing
- Correction for population structure and linkage disequilibrium
- Effect size quantification and correlation analysis

## Dependencies and Software

### R packages
- `data.table`, `dplyr` - Data manipulation
- `rtracklayer` - Genome coordinate conversion
- `GenomicRanges` - Genomic interval operations
- `regioneR` - Permutation testing
- `ggplot2`, `ggpubr` - Data visualization

### External software
- **PLINK2** - Genetic association analysis
- **Python libraries** - Akita model implementation (numpy, pandas, matplotlib, scipy)

### Data resources
- GIANT Height GWAS summary statistics
- 1000 Genomes Project reference panels
- TOPMed reference data
- UCSC Genome Browser annotations
- PhastCons conservation scores

## Output Files
- Population-specific association results (`*_results.tsv`)
- 3D genome divergence metrics (`max_divergence.csv`)
- Enrichment test results and visualizations
- Hi-C contact map visualizations and comparative plots

## Key Analysis Workflow

1. **Data Preparation**: GWAS summary statistics are lifted over to hg38 and formatted for analysis
2. **Variant Clumping**: Independent variants are identified using different reference panels
3. **Haplotype Extraction**: Phased haplotypes are extracted from 1000 Genomes data
4. **Sequence Generation**: Haplotype sequences are prepared for Akita model input
5. **3D Genome Prediction**: Akita model predicts chromatin contact maps for each haplotype
6. **Divergence Analysis**: Correlation between contact maps is calculated and divergence computed
7. **Statistical Analysis**: Relationship between GWAS effect sizes and 3D genome effects is assessed
8. **Enrichment Testing**: Permutation tests evaluate enrichment of 3D-altering variants in GWAS loci

## R Package Dependencies

The scripts require the following R packages:
- **Data manipulation**: `dplyr`, `data.table`, `purrr`
- **Genomics**: `rtracklayer`, `GenomicRanges`, `Biostrings`, `vcfR`
- **Statistics**: `regioneR`
- **Visualization**: `ggplot2`, `ggpubr`, `cowplot`, `gridExtra`, `patchwork`
- **File I/O**: `readxl`, `writexl`
- **Bioinformatics**: `biomaRt`
- **String operations**: `stringdist`
- **Scales**: `scales`

## External Tools

- **PLINK2**: For genetic association analysis and clumping
- **BCFtools**: For VCF file manipulation
- **Akita model**: For predicting 3D genome structure
- **LDSC**: For heritability enrichment analysis

This comprehensive analytical pipeline enables systematic investigation of how genetic variants affecting complex traits alter 3D genome architecture through integration of population genetics, machine learning, and 3D genomics approaches.
