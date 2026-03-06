# GWAS Catalog VCF Files

Variants from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/) converted to VCF
format (VCFv4.2, GRCh38/hg38), with REF/ALT alleles annotated from dbSNP via
multiple annotation pipelines.

---

## Output files

### `gwas_catalog.vcf.gz` + `gwas_catalog.vcf.gz.tbi`

| Property | Value |
|---|---|
| Format | VCFv4.2, bgzip-compressed, tabix-indexed |
| Genome build | GRCh38 / hg38 |
| Records | 886,598 variant–association rows |
| Unique rsIDs | 402,612 |
| Unique traits | 41,587 |
| Unique mapped genes | 42,732 |
| REF annotated | 98.6% |
| ALT annotated | 99.99% |
| REF == ALT (residual) | 0.06% |

All GWAS Catalog associations with a valid GRCh38 chromosome and position.
One row per variant–trait–study combination, so the same variant will appear
multiple times if it is associated with multiple traits or reported in multiple
studies.

### `gwas_catalog_significant.vcf.gz` + `gwas_catalog_significant.vcf.gz.tbi`

| Property | Value |
|---|---|
| Format | VCFv4.2, bgzip-compressed, tabix-indexed |
| Genome build | GRCh38 / hg38 |
| Records | 751,989 variant–association rows |
| Significance threshold | p < 5×10⁻⁸ (PVALUE_MLOG ≥ 7.301) |

Subset of `gwas_catalog.vcf.gz` containing only associations that reach the
conventional genome-wide significance threshold. All annotation columns are
identical to the full file.

---

## VCF structure

### Fixed columns

| Column | Description |
|---|---|
| `CHROM` | Chromosome in UCSC style (`chr1` … `chr22`, `chrX`, `chrY`) |
| `POS` | GRCh38 base-pair position (1-based) |
| `ID` | rsID from the GWAS Catalog `SNPS` column; falls back to `rs{SNP_ID_CURRENT}` |
| `REF` | Reference allele from dbSNP/GRCh38 (annotated post-hoc; `.` if unresolvable) |
| `ALT` | Alternate allele; the GWAS risk allele where it differs from REF, otherwise the dbSNP ALT |
| `QUAL` | `.` (not applicable for GWAS data) |
| `FILTER` | `.` (not applicable for GWAS data) |

### INFO field annotations

Every record carries the following INFO tags sourced from the GWAS Catalog
download files. Column header documentation:
<https://www.ebi.ac.uk/gwas/docs/fileheaders>

| Tag | Type | Source column | Description |
|---|---|---|---|
| `TRAIT` | String | `DISEASE/TRAIT` | Disease or trait associated with the variant |
| `MAPPED_GENE` | String | `MAPPED_GENE` | Gene(s) mapped by Ensembl; intergenic variants show upstream − downstream |
| `REPORTED_GENE` | String | `REPORTED GENE(S)` | Gene(s) reported by the original publication |
| `CONTEXT` | String | `CONTEXT` | Ensembl VEP functional consequence (e.g. `intergenic_variant`, `missense_variant`) |
| `RISK_ALLELE` | String | `STRONGEST SNP-RISK ALLELE` | Risk allele in `rsID-Allele` format |
| `RISK_ALLELE_FREQ` | Float | `RISK ALLELE FREQUENCY` | Risk allele frequency in the discovery sample |
| `PVALUE` | String | `P-VALUE` | Association p-value (e.g. `5E-8`) |
| `PVALUE_MLOG` | Float | `PVALUE_MLOG` | −log₁₀(p-value) |
| `PVALUE_TEXT` | String | `P-VALUE (TEXT)` | Qualifier for the p-value (e.g. age group, condition) |
| `OR_BETA` | String | `OR or BETA` | Odds ratio or beta coefficient |
| `CI_95` | String | `95% CI (TEXT)` | 95% confidence interval for OR/beta |
| `STUDY` | String | `STUDY` | Title of the GWAS publication |
| `PUBMEDID` | Integer | `PUBMEDID` | PubMed ID |
| `FIRST_AUTHOR` | String | `FIRST AUTHOR` | First author surname + initials |
| `JOURNAL` | String | `JOURNAL` | Abbreviated journal name |
| `DATE` | String | `DATE` | Publication date (YYYY-MM-DD) |
| `STUDY_ACCESSION` | String | ancestry file `STUDY ACCESSION` | GWAS Catalog accession(s) (GCSTXXXXXX) |
| `BROAD_ANCESTRY` | String | ancestry file `BROAD ANCESTRAL CATEGORY` | Broad ancestry of the study sample (e.g. `European`, `East_Asian`); comma-separated if multiple |
| `COUNTRY_OF_ORIGIN` | String | ancestry file `COUNTRY OF ORIGIN` | Country of origin of participants |
| `ADDITIONAL_ANCESTRY` | String | ancestry file `ADDITONAL ANCESTRY DESCRIPTION` | Additional ancestry detail beyond the broad category |
| `INITIAL_SAMPLE_SIZE` | String | `INITIAL SAMPLE SIZE` | Number and ancestry of the discovery sample |
| `REPLICATION_SAMPLE_SIZE` | String | `REPLICATION SAMPLE SIZE` | Number and ancestry of the replication sample |
| `PLATFORM` | String | `PLATFORM [SNPS PASSING QC]` | Genotyping platform and SNP count after QC |
| `SNP_ID_CURRENT` | String | `SNP_ID_CURRENT` | Current rsID numeric part after dbSNP merges |
| `INTERGENIC` | Integer | `INTERGENIC` | `1` if intergenic, `0` otherwise |
| `UPSTREAM_GENE_ID` | String | `UPSTREAM_GENE_ID` | Ensembl ID of the nearest upstream gene (intergenic only) |
| `DOWNSTREAM_GENE_ID` | String | `DOWNSTREAM_GENE_ID` | Ensembl ID of the nearest downstream gene (intergenic only) |
| `UPSTREAM_GENE_DIST` | Integer | `UPSTREAM_GENE_DISTANCE` | Distance in bp to the nearest upstream gene |
| `DOWNSTREAM_GENE_DIST` | Integer | `DOWNSTREAM_GENE_DISTANCE` | Distance in bp to the nearest downstream gene |
| `CNV` | String | `CNV` | Copy number variant flag (`Y`/`N`) |
| `REGION` | String | `REGION` | Cytogenetic region (e.g. `11p15.4`) |

> INFO field values have spaces replaced by underscores and special characters
> percent-encoded (`=` → `%3D`, `;` → `%3B`) to conform to the VCF spec.

---

## REF/ALT annotation notes

The GWAS Catalog download does not provide reference or alternate alleles.
Alleles were annotated in three stages:

1. **myvariant.info** (dbSNP/GRCh38) — batch rsID queries; 98.3% hit rate
2. **Ensembl REST API** `/variation/homo_sapiens` — fallback for rsIDs not in myvariant.info
3. **Ensembl BioMart** (`hsapiens_snp`) via R — secondary fallback

When the GWAS risk allele equals the GRCh38 reference allele (a known GWAS
Catalog quirk), the ALT column is set to the true dbSNP alternate allele so
that REF ≠ ALT. The original risk allele is always preserved in the
`RISK_ALLELE` INFO tag regardless.

Residual missing values (REF = `.` or ALT = `.`) represent variants where
none of the annotation sources returned a result — typically very rare or
recently-submitted variants not yet in dbSNP.

---

## Source files

| File | Description |
|---|---|
| `gwas-catalog-download-associations-v1.0-full.tsv` | Raw GWAS Catalog associations download (513 MB) |
| `gwas_catalog-ancestry_r2026-02-16.tsv` | GWAS Catalog ancestry/study metadata (44 MB) |

---

## Scripts

| Script | Language | Role |
|---|---|---|
| `make_vcf.py` | Python (uv) | Merge TSVs, build initial VCF with INFO fields |
| `annotate_alleles.py` | Python (uv) | Annotate REF/ALT via myvariant.info and Ensembl REST |
| `annotate_biomart.R` | R | BioMart fallback for rsIDs not resolved by Python steps |
| `fix_ref_eq_alt.py` | Python (uv) | Fix records where risk allele == REF using dbSNP ALT |
| `compress_vcf.py` | Python (uv) | bgzip-compress and tabix-index all VCF files |
| `validate_vcf.R` | R | QC summary statistics on the final VCF |

### Regeneration pipeline

```bash
# 1. Build initial VCF
uv run make_vcf.py

# 2. Annotate REF/ALT (myvariant.info + Ensembl REST fallback)
uv run annotate_alleles.py

# 3. biomaRt fallback for remaining unresolved rsIDs
Rscript annotate_biomart.R

# 4. Merge R annotations back
uv run annotate_alleles.py --apply-r-only

# 5. Fix records where risk allele == REF
uv run fix_ref_eq_alt.py

# 6. bgzip compress + tabix index
uv run compress_vcf.py

# 7. QC validation
Rscript validate_vcf.R
```

---

## Usage examples

Query a specific genomic region with tabix:
```bash
tabix gwas_catalog.vcf.gz chr6:26000000-34000000
```

Extract all variants associated with a specific trait (using bcftools):
```bash
bcftools view gwas_catalog.vcf.gz | grep "TRAIT=Type_2_diabetes"
```

Load in R with VariantAnnotation:
```r
library(VariantAnnotation)
vcf <- readVcf("gwas_catalog.vcf.gz", genome = "hg38")
```

Load in Python with pysam:
```python
import pysam
tbx = pysam.TabixFile("gwas_catalog.vcf.gz")
for row in tbx.fetch("chr6", 26000000, 34000000):
    print(row)
```
