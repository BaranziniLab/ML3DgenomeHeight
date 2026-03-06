#!/usr/bin/env python3
"""
Convert GWAS Catalog TSV files into a VCF (Variant Call Format) file.

Sources:
  - gwas-catalog-download-associations-v1.0-full.tsv  (associations)
  - gwas_catalog-ancestry_r2026-02-16.tsv             (ancestry/study metadata)

All genomic coordinates are reported in GRCh38 / hg38 as provided by the
GWAS Catalog (https://www.ebi.ac.uk/gwas/).

Column header documentation:
  https://www.ebi.ac.uk/gwas/docs/fileheaders
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

# ── paths ──────────────────────────────────────────────────────────────────────
HERE = Path(__file__).parent
ASSOC_FILE = HERE / "gwas-catalog-download-associations-v1.0-full.tsv"
ANCESTRY_FILE = HERE / "gwas_catalog-ancestry_r2026-02-16.tsv"
OUT_VCF = HERE / "gwas_catalog.vcf"


# ── helpers ────────────────────────────────────────────────────────────────────
def vcf_escape(value: str) -> str:
    """Escape special characters for VCF INFO field values."""
    if pd.isna(value):
        return "."
    s = str(value).strip()
    if not s:
        return "."
    # VCF INFO values: replace semicolons, equals-signs and spaces that would
    # break INFO parsing; keep commas (they separate multiple values).
    s = s.replace("%", "%25")
    s = s.replace(";", "%3B")
    s = s.replace("=", "%3D")
    s = s.replace(" ", "_")
    s = s.replace("\t", "%09")
    return s if s else "."


def extract_risk_allele(snp_risk: str) -> str:
    """Extract allele from 'rsXXXX-A' style strings; return '.' if absent."""
    if pd.isna(snp_risk):
        return "."
    m = re.search(r"-([ACGTacgt?]+)$", str(snp_risk).strip())
    return m.group(1).upper() if m else "."


def normalise_chrom(chrom) -> str | None:
    """Return UCSC-style chromosome string (chr1 … chrX/Y/MT) or None."""
    if pd.isna(chrom):
        return None
    c = str(chrom).strip().upper()
    if not c or c in ("NA", "NAN", "NONE", ""):
        return None
    c = c.replace("CHR", "")
    # keep only valid chromosomes
    valid = set([str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"])
    if c not in valid:
        return None
    if c == "M":
        c = "MT"
    return f"chr{c}"


def safe_int(val) -> int | None:
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return None


# ── load associations ──────────────────────────────────────────────────────────
print("Loading associations …", file=sys.stderr)
assoc = pd.read_csv(
    ASSOC_FILE,
    sep="\t",
    dtype=str,
    low_memory=False,
    na_values=["", "NA", "NaN", "nan", "N/A"],
)
assoc.columns = assoc.columns.str.strip()
print(f"  {len(assoc):,} association rows loaded", file=sys.stderr)

# ── load ancestry ──────────────────────────────────────────────────────────────
print("Loading ancestry …", file=sys.stderr)
ancestry = pd.read_csv(
    ANCESTRY_FILE,
    sep="\t",
    dtype=str,
    low_memory=False,
    na_values=["", "NA", "NaN", "nan", "N/A"],
)
ancestry.columns = ancestry.columns.str.strip()
print(f"  {len(ancestry):,} ancestry rows loaded", file=sys.stderr)

# ── collapse ancestry per (PUBMEDID, STAGE) → one row per study ───────────────
# We want broad ancestral categories (all stages, comma-joined if multiple)
ancestry_collapsed = (
    ancestry.groupby("PUBMEDID", dropna=False)
    .apply(
        lambda g: pd.Series(
            {
                "BROAD_ANCESTRY": ",".join(
                    sorted(set(g["BROAD ANCESTRAL CATEGORY"].dropna().tolist()))
                ),
                "COUNTRY_OF_ORIGIN": ",".join(
                    sorted(set(g["COUNTRY OF ORIGIN"].dropna().tolist()))
                ),
                "STUDY_ACCESSION": ",".join(
                    sorted(set(g["STUDY ACCESSION"].dropna().tolist()))
                ),
                "ADDITIONAL_ANCESTRY": ",".join(
                    sorted(
                        set(g["ADDITONAL ANCESTRY DESCRIPTION"].dropna().tolist())
                    )
                ),
            }
        ),
        include_groups=False,
    )
    .reset_index()
)

# ── merge ──────────────────────────────────────────────────────────────────────
print("Merging associations with ancestry …", file=sys.stderr)
df = assoc.merge(ancestry_collapsed, on="PUBMEDID", how="left")
print(f"  {len(df):,} rows after merge", file=sys.stderr)

# ── filter rows with usable chromosome + position ─────────────────────────────
df["_chrom"] = df["CHR_ID"].apply(normalise_chrom)
df["_pos"] = df["CHR_POS"].apply(safe_int)
df_valid = df[df["_chrom"].notna() & df["_pos"].notna()].copy()
print(
    f"  {len(df_valid):,} rows with valid GRCh38 coordinates "
    f"(dropped {len(df) - len(df_valid):,})",
    file=sys.stderr,
)

# ── deduplicate: keep one row per (chrom, pos, rsID, trait) ───────────────────
# When the same variant is reported in multiple studies for the same trait,
# keep all rows (they differ by study). Sort for deterministic output.
df_valid = df_valid.sort_values(["_chrom", "_pos", "SNPS", "DISEASE/TRAIT"])

# ── extract ALT allele from STRONGEST SNP-RISK ALLELE ─────────────────────────
df_valid["_alt"] = df_valid["STRONGEST SNP-RISK ALLELE"].apply(extract_risk_allele)


# ── build INFO string ──────────────────────────────────────────────────────────
def build_info(row) -> str:
    fields = []

    def add(key, value):
        v = vcf_escape(value)
        if v != ".":
            fields.append(f"{key}={v}")

    add("TRAIT", row.get("DISEASE/TRAIT"))
    add("MAPPED_GENE", row.get("MAPPED_GENE"))
    add("REPORTED_GENE", row.get("REPORTED GENE(S)"))
    add("CONTEXT", row.get("CONTEXT"))
    add("RISK_ALLELE", row.get("STRONGEST SNP-RISK ALLELE"))
    add("RISK_ALLELE_FREQ", row.get("RISK ALLELE FREQUENCY"))
    add("PVALUE", row.get("P-VALUE"))
    add("PVALUE_MLOG", row.get("PVALUE_MLOG"))
    add("PVALUE_TEXT", row.get("P-VALUE (TEXT)"))
    add("OR_BETA", row.get("OR or BETA"))
    add("CI_95", row.get("95% CI (TEXT)"))
    add("STUDY", row.get("STUDY"))
    add("PUBMEDID", row.get("PUBMEDID"))
    add("FIRST_AUTHOR", row.get("FIRST AUTHOR"))
    add("JOURNAL", row.get("JOURNAL"))
    add("DATE", row.get("DATE"))
    add("STUDY_ACCESSION", row.get("STUDY_ACCESSION"))
    add("BROAD_ANCESTRY", row.get("BROAD_ANCESTRY"))
    add("COUNTRY_OF_ORIGIN", row.get("COUNTRY_OF_ORIGIN"))
    add("ADDITIONAL_ANCESTRY", row.get("ADDITIONAL_ANCESTRY"))
    add("INITIAL_SAMPLE_SIZE", row.get("INITIAL SAMPLE SIZE"))
    add("REPLICATION_SAMPLE_SIZE", row.get("REPLICATION SAMPLE SIZE"))
    add("PLATFORM", row.get("PLATFORM [SNPS PASSING QC]"))
    add("SNP_ID_CURRENT", row.get("SNP_ID_CURRENT"))
    add("INTERGENIC", row.get("INTERGENIC"))
    add("UPSTREAM_GENE_ID", row.get("UPSTREAM_GENE_ID"))
    add("DOWNSTREAM_GENE_ID", row.get("DOWNSTREAM_GENE_ID"))
    add("UPSTREAM_GENE_DIST", row.get("UPSTREAM_GENE_DISTANCE"))
    add("DOWNSTREAM_GENE_DIST", row.get("DOWNSTREAM_GENE_DISTANCE"))
    add("CNV", row.get("CNV"))
    add("REGION", row.get("REGION"))

    return ";".join(fields) if fields else "."


print("Building INFO fields …", file=sys.stderr)
df_valid["_info"] = df_valid.apply(build_info, axis=1)

# ── VCF ID column: prefer rsID from SNPS column, fallback to SNP_ID_CURRENT ──
def vcf_id(row) -> str:
    snp = str(row.get("SNPS", "")).strip()
    cur = str(row.get("SNP_ID_CURRENT", "")).strip()
    if snp and snp not in ("nan", "NA", ""):
        return snp
    if cur and cur not in ("nan", "NA", ""):
        return f"rs{cur}"
    return "."


df_valid["_id"] = df_valid.apply(vcf_id, axis=1)


# ── write VCF ─────────────────────────────────────────────────────────────────
print(f"Writing {OUT_VCF} …", file=sys.stderr)

# Chromosome sort order
CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrMT": 25})
df_valid["_chrom_sort"] = df_valid["_chrom"].map(CHROM_ORDER).fillna(99)
df_valid = df_valid.sort_values(["_chrom_sort", "_pos"])

with open(OUT_VCF, "w") as fh:
    # ── meta-information lines ──
    fh.write("##fileformat=VCFv4.2\n")
    fh.write("##fileDate=20260304\n")
    fh.write("##source=GWAS_Catalog_v1.0_full_download\n")
    fh.write("##reference=GRCh38/hg38\n")
    fh.write(
        "##INFO=<ID=TRAIT,Number=.,Type=String,"
        'Description="Disease or trait associated with the variant '
        "(GWAS Catalog 'DISEASE/TRAIT' column). "
        "Spaces replaced by underscores."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=MAPPED_GENE,Number=.,Type=String,"
        'Description="Gene(s) mapped to the variant by Ensembl. '
        "Intergenic variants show upstream/downstream genes separated by ' - '. "
        "Source: GWAS Catalog 'MAPPED_GENE' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=REPORTED_GENE,Number=.,Type=String,"
        'Description="Gene(s) reported by the original publication '
        "(GWAS Catalog 'REPORTED GENE(S)' column)."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=CONTEXT,Number=1,Type=String,"
        'Description="Variant functional consequence/context annotated by Ensembl VEP '
        "(e.g. intergenic_variant, missense_variant). "
        "Source: GWAS Catalog 'CONTEXT' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=RISK_ALLELE,Number=1,Type=String,"
        'Description="Strongest SNP risk allele in the format rsID-Allele. '
        "Source: GWAS Catalog 'STRONGEST SNP-RISK ALLELE' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=RISK_ALLELE_FREQ,Number=1,Type=Float,"
        'Description="Reported risk allele frequency in the discovery sample. '
        "Source: GWAS Catalog 'RISK ALLELE FREQUENCY' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=PVALUE,Number=1,Type=String,"
        'Description="p-value of the most significant SNP reported for the association. '
        "Source: GWAS Catalog 'P-VALUE' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=PVALUE_MLOG,Number=1,Type=Float,"
        'Description="-log10(p-value). '
        "Source: GWAS Catalog 'PVALUE_MLOG' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=PVALUE_TEXT,Number=.,Type=String,"
        'Description="Qualifier for the p-value (e.g. condition, age group). '
        "Source: GWAS Catalog 'P-VALUE (TEXT)' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=OR_BETA,Number=1,Type=String,"
        'Description="Odds ratio or beta coefficient for the association. '
        "Source: GWAS Catalog 'OR or BETA' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=CI_95,Number=.,Type=String,"
        'Description="95%% confidence interval for the OR/Beta. '
        "Source: GWAS Catalog '95%% CI (TEXT)' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=STUDY,Number=.,Type=String,"
        'Description="Title of the GWAS publication. '
        "Source: GWAS Catalog 'STUDY' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=PUBMEDID,Number=1,Type=Integer,"
        'Description="PubMed ID of the publication. '
        "Source: GWAS Catalog 'PUBMEDID' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=FIRST_AUTHOR,Number=1,Type=String,"
        'Description="First author surname + initials. '
        "Source: GWAS Catalog 'FIRST AUTHOR' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=JOURNAL,Number=1,Type=String,"
        'Description="Abbreviated journal name. '
        "Source: GWAS Catalog 'JOURNAL' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=DATE,Number=1,Type=String,"
        'Description="Publication date (YYYY-MM-DD). '
        "Source: GWAS Catalog 'DATE' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=STUDY_ACCESSION,Number=.,Type=String,"
        'Description="GWAS Catalog study accession(s) (GCSTXXXXXX). '
        "Source: GWAS Catalog ancestry file 'STUDY ACCESSION' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=BROAD_ANCESTRY,Number=.,Type=String,"
        'Description="Broad ancestral category of the study sample '
        "(e.g. European, East_Asian). Comma-separated if multiple ancestries. "
        "Source: GWAS Catalog ancestry file 'BROAD ANCESTRAL CATEGORY' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=COUNTRY_OF_ORIGIN,Number=.,Type=String,"
        'Description="Country of origin of study participants. '
        "Source: GWAS Catalog ancestry file 'COUNTRY OF ORIGIN' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=ADDITIONAL_ANCESTRY,Number=.,Type=String,"
        'Description="Additional ancestry description beyond the broad category. '
        "Source: GWAS Catalog ancestry file 'ADDITONAL ANCESTRY DESCRIPTION' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=INITIAL_SAMPLE_SIZE,Number=.,Type=String,"
        'Description="Number and ancestry of initial (GWAS) sample. '
        "Source: GWAS Catalog 'INITIAL SAMPLE SIZE' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=REPLICATION_SAMPLE_SIZE,Number=.,Type=String,"
        'Description="Number and ancestry of replication sample. '
        "Source: GWAS Catalog 'REPLICATION SAMPLE SIZE' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=PLATFORM,Number=.,Type=String,"
        'Description="Genotyping platform and number of SNPs passing QC. '
        "Source: GWAS Catalog 'PLATFORM [SNPS PASSING QC]' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=SNP_ID_CURRENT,Number=1,Type=String,"
        'Description="Current rsID (numeric part only) after dbSNP merges. '
        "Source: GWAS Catalog 'SNP_ID_CURRENT' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=INTERGENIC,Number=1,Type=Integer,"
        'Description="Flag: 1 if the variant is intergenic, 0 otherwise. '
        "Source: GWAS Catalog 'INTERGENIC' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=UPSTREAM_GENE_ID,Number=1,Type=String,"
        'Description="Ensembl gene ID of the nearest upstream gene for intergenic variants. '
        "Source: GWAS Catalog 'UPSTREAM_GENE_ID' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=DOWNSTREAM_GENE_ID,Number=1,Type=String,"
        'Description="Ensembl gene ID of the nearest downstream gene for intergenic variants. '
        "Source: GWAS Catalog 'DOWNSTREAM_GENE_ID' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=UPSTREAM_GENE_DIST,Number=1,Type=Integer,"
        'Description="Distance in bp to the nearest upstream gene (intergenic only). '
        "Source: GWAS Catalog 'UPSTREAM_GENE_DISTANCE' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=DOWNSTREAM_GENE_DIST,Number=1,Type=Integer,"
        'Description="Distance in bp to the nearest downstream gene (intergenic only). '
        "Source: GWAS Catalog 'DOWNSTREAM_GENE_DISTANCE' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=CNV,Number=1,Type=String,"
        'Description="Copy number variant flag (Y/N). '
        "Source: GWAS Catalog 'CNV' column."
        '">\n'
    )
    fh.write(
        "##INFO=<ID=REGION,Number=1,Type=String,"
        'Description="Cytogenetic region of the variant. '
        "Source: GWAS Catalog 'REGION' column."
        '">\n'
    )
    fh.write(
        "##NOTE=Coordinates are in GRCh38/hg38 as provided by the GWAS Catalog. "
        "REF alleles are not available from the GWAS Catalog download; REF is set to '.' "
        "for all records. ALT is set to the extracted risk allele where determinable, "
        "otherwise '.'. INFO field values have spaces replaced by underscores and "
        "special characters percent-encoded per VCF spec.\n"
    )
    fh.write(
        "##GWAS_Catalog_URL=https://www.ebi.ac.uk/gwas/\n"
    )
    fh.write(
        "##GWAS_Catalog_fileheaders_docs=https://www.ebi.ac.uk/gwas/docs/fileheaders\n"
    )

    # ── column header ──
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # ── data rows ──
    written = 0
    for _, row in df_valid.iterrows():
        chrom = row["_chrom"]
        pos = int(row["_pos"])
        vid = row["_id"]
        alt = row["_alt"]
        info = row["_info"]
        fh.write(f"{chrom}\t{pos}\t{vid}\t.\t{alt}\t.\t.\t{info}\n")
        written += 1

print(f"Done — wrote {written:,} variant records to {OUT_VCF}", file=sys.stderr)
