#!/usr/bin/env python3
"""
compress_vcf.py

bgzip-compress and tabix-index all *.vcf files in the script's directory.
Uses pysam (bundles htslib bgzip + tabix) — no system-level htslib needed.

Output: <name>.vcf.gz + <name>.vcf.gz.tbi for each input VCF.
"""

import sys
import logging
from pathlib import Path

import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

HERE = Path(__file__).parent
vcf_files = sorted(HERE.glob("*.vcf"))

if not vcf_files:
    log.info("No *.vcf files found.")
    sys.exit(0)

for vcf in vcf_files:
    gz_path = vcf.with_suffix(".vcf.gz")
    tbi_path = Path(str(gz_path) + ".tbi")

    # ── bgzip compress ────────────────────────────────────────────────────────
    log.info(f"Compressing {vcf.name} → {gz_path.name} …")
    pysam.tabix_compress(
        str(vcf),
        str(gz_path),
        force=True,   # overwrite if already exists
    )
    log.info(f"  done  ({gz_path.stat().st_size / 1e9:.2f} GB)")

    # ── tabix index ───────────────────────────────────────────────────────────
    log.info(f"Indexing {gz_path.name} …")
    pysam.tabix_index(
        str(gz_path),
        preset="vcf",   # tells tabix: CHROM=col1, POS=col2, skip #-lines
        force=True,
    )
    log.info(f"  index written: {tbi_path.name}")

log.info("All VCFs compressed and indexed.")
