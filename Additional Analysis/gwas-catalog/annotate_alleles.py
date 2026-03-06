#!/usr/bin/env python3
"""
annotate_alleles.py

Adds proper REF/ALT alleles to gwas_catalog.vcf and writes:
  1. gwas_catalog.vcf          (all records, REF/ALT filled where possible)
  2. gwas_catalog_significant.vcf  (p < 5e-8 only, same annotation)

Annotation strategy (in order of priority):
  1. myvariant.info  — batch rsID lookup against dbSNP/GRCh38 (primary)
  2. Ensembl REST API POST /variation — batch fallback for rsIDs not in myvariant
  3. annotate_biomart.R — biomaRt queries run externally; results read back
     from r_annotations.tsv if that file exists

REF/ALT assignment rules:
  • REF = reference allele from dbSNP/GRCh38
  • ALT = alternate allele(s); if the GWAS risk allele is already the REF,
    ALT is set to the other known allele (or "." if unknown)
  • Positions and coordinates remain GRCh38 as-is

Usage:
  uv run annotate_alleles.py [--skip-myvariant] [--skip-ensembl]
"""

import re
import sys
import time
import json
import logging
import argparse
from pathlib import Path
from collections import defaultdict

import requests
import myvariant

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

HERE = Path(__file__).parent
VCF_IN       = HERE / "gwas_catalog.vcf"
VCF_OUT      = HERE / "gwas_catalog.vcf"          # overwrite
VCF_SIG_OUT  = HERE / "gwas_catalog_significant.vcf"
MISSING_TSV  = HERE / "missing_annotations.tsv"    # written for R fallback
R_ANNOT_TSV  = HERE / "r_annotations.tsv"          # read back from R

P_SIG_MLOG = 7.301029995663981   # -log10(5e-8)

VALID_BASES = set("ACGTacgt")

# ── helpers ────────────────────────────────────────────────────────────────────

def parse_pvalue_mlog(info: str) -> float | None:
    """Extract PVALUE_MLOG from INFO string."""
    m = re.search(r"PVALUE_MLOG=([^;]+)", info)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    # fallback: try to parse raw PVALUE
    m2 = re.search(r"PVALUE=([^;]+)", info)
    if m2:
        try:
            p = float(m2.group(1))
            if p > 0:
                import math
                return -math.log10(p)
        except ValueError:
            pass
    return None


def is_valid_allele(a: str) -> bool:
    return bool(a) and a != "." and all(c in VALID_BASES for c in a.upper())


def clean_allele(a) -> str:
    """Normalise an allele string to uppercase; return '.' if invalid."""
    if not a or str(a).strip() in ("", ".", "nan", "None"):
        return "."
    a = str(a).strip().upper()
    # drop ambiguous IUPAC / multi-allelic separators for simple representation
    # keep the longest unambiguous allele (first one if comma-separated)
    if "," in a:
        a = a.split(",")[0].strip()
    if all(c in VALID_BASES for c in a):
        return a
    return "."


# ── phase 1 — scan VCF ────────────────────────────────────────────────────────

def scan_vcf(path: Path) -> tuple[int, dict[str, tuple[str, int, str]]]:
    """
    Return (meta_line_count, {rsid: (chrom, pos, current_alt)}).
    Only collects rsIDs whose current REF == '.'.
    """
    meta_count = 0
    needs: dict[str, tuple[str, int, str]] = {}

    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                meta_count += 1
                continue
            if line.startswith("#CHROM"):
                meta_count += 1   # count the header line too
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos_s, vid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            if ref != ".":
                continue  # already annotated
            if vid and vid != ".":
                try:
                    pos = int(pos_s)
                except ValueError:
                    continue
                # For repeated rsIDs (same variant, many traits), store once
                if vid not in needs:
                    needs[vid] = (chrom, pos, alt)
    return meta_count, needs


# ── phase 2 — myvariant.info ──────────────────────────────────────────────────

def query_myvariant(rsids: list[str]) -> dict[str, tuple[str, str]]:
    """
    Batch-query myvariant.info.
    Returns dict: rsid → (REF, ALT).  Both strings are uppercase ATCG or '.'.
    """
    mv = myvariant.MyVariantInfo()
    annot: dict[str, tuple[str, str]] = {}

    batch_size = 1000
    total = len(rsids)
    log.info(f"myvariant.info: querying {total:,} rsIDs in batches of {batch_size}")

    for start in range(0, total, batch_size):
        batch = rsids[start : start + batch_size]
        if start % 10000 == 0:
            log.info(f"  myvariant progress: {start:,}/{total:,}  annotated so far: {len(annot):,}")

        try:
            results = mv.querymany(
                batch,
                scopes="dbsnp.rsid",
                fields="dbsnp.ref,dbsnp.alt,dbsnp.alleles",
                as_dataframe=False,
                verbose=False,
            )
        except Exception as e:
            log.warning(f"  myvariant batch error at {start}: {e}")
            time.sleep(2)
            continue

        for r in results:
            if r.get("notfound") or "error" in r:
                continue
            query_id = r.get("query", "")
            dbsnp = r.get("dbsnp")
            if not dbsnp:
                continue

            ref = clean_allele(dbsnp.get("ref"))
            alt_raw = dbsnp.get("alt")

            # alt can be a string, list, or missing
            if isinstance(alt_raw, list):
                alt = clean_allele(alt_raw[0]) if alt_raw else "."
            elif alt_raw:
                alt = clean_allele(str(alt_raw))
            else:
                # try alleles array
                alleles = dbsnp.get("alleles")
                if isinstance(alleles, list):
                    # pick first allele that differs from REF
                    alt = "."
                    for a_obj in alleles:
                        a = clean_allele(a_obj.get("allele", ".") if isinstance(a_obj, dict) else str(a_obj))
                        if a != ref and is_valid_allele(a):
                            alt = a
                            break
                else:
                    alt = "."

            if is_valid_allele(ref) and query_id:
                annot[query_id] = (ref, alt)

    log.info(f"myvariant.info: annotated {len(annot):,}/{total:,} rsIDs")
    return annot


# ── phase 3 — Ensembl REST fallback ──────────────────────────────────────────

ENSEMBL_VAR_POST = "https://rest.ensembl.org/variation/homo_sapiens"
ENSEMBL_HEADERS  = {"Content-Type": "application/json", "Accept": "application/json"}
ENSEMBL_BATCH    = 200   # POST endpoint accepts up to 200 IDs


def query_ensembl_batch(rsids: list[str]) -> dict[str, tuple[str, str]]:
    """Single POST to Ensembl variation endpoint for up to 200 rsIDs."""
    annot: dict[str, tuple[str, str]] = {}
    payload = {"ids": rsids}
    try:
        r = requests.post(
            ENSEMBL_VAR_POST,
            headers=ENSEMBL_HEADERS,
            data=json.dumps(payload),
            timeout=30,
        )
        if r.status_code == 429:
            retry_after = int(r.headers.get("Retry-After", 5))
            log.warning(f"Ensembl rate-limited; sleeping {retry_after}s")
            time.sleep(retry_after)
            return annot
        if r.status_code != 200:
            return annot
        data = r.json()
        for rsid, info in data.items():
            for mapping in info.get("mappings", []):
                allele_str = mapping.get("allele_string", "")
                if "/" in allele_str:
                    parts = allele_str.split("/")
                    ref = clean_allele(parts[0])
                    alt = clean_allele(parts[1]) if len(parts) > 1 else "."
                    if is_valid_allele(ref):
                        annot[rsid] = (ref, alt)
                        break
    except Exception as e:
        log.debug(f"Ensembl batch error: {e}")
    return annot


def query_ensembl(rsids: list[str]) -> dict[str, tuple[str, str]]:
    """Batch-query Ensembl REST API for rsIDs not found in myvariant."""
    annot: dict[str, tuple[str, str]] = {}
    total = len(rsids)
    log.info(f"Ensembl REST: querying {total:,} remaining rsIDs")

    for start in range(0, total, ENSEMBL_BATCH):
        batch = rsids[start : start + ENSEMBL_BATCH]
        annot.update(query_ensembl_batch(batch))
        # Ensembl asks for ≤15 requests/second
        time.sleep(0.07)
        if start % 2000 == 0 and start > 0:
            log.info(f"  Ensembl progress: {start:,}/{total:,}  annotated so far: {len(annot):,}")

    log.info(f"Ensembl REST: annotated {len(annot):,}/{total:,} rsIDs")
    return annot


# ── phase 4 — read R / biomaRt annotations if present ────────────────────────

def read_r_annotations(path: Path) -> dict[str, tuple[str, str]]:
    """Read r_annotations.tsv produced by annotate_biomart.R."""
    annot: dict[str, tuple[str, str]] = {}
    if not path.exists():
        return annot
    log.info(f"Reading R biomaRt annotations from {path.name}")
    with open(path) as fh:
        header = fh.readline()  # rsid\tref\talt
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            rsid, ref, alt = parts[0].strip(), parts[1].strip(), parts[2].strip()
            ref = clean_allele(ref)
            alt = clean_allele(alt)
            if is_valid_allele(ref) and rsid:
                annot[rsid] = (ref, alt)
    log.info(f"  R annotations loaded: {len(annot):,}")
    return annot


# ── phase 5 — write missing list for R ────────────────────────────────────────

def write_missing_for_r(
    needs: dict[str, tuple[str, int, str]],
    already_annotated: set[str],
    path: Path,
):
    still_missing = {k: v for k, v in needs.items() if k not in already_annotated}
    if not still_missing:
        log.info("No variants need R/biomaRt fallback.")
        return 0
    log.info(f"Writing {len(still_missing):,} rsIDs to {path.name} for biomaRt fallback")
    with open(path, "w") as fh:
        fh.write("rsid\tchrom\tpos\n")
        for rsid, (chrom, pos, _) in still_missing.items():
            fh.write(f"{rsid}\t{chrom}\t{pos}\n")
    return len(still_missing)


# ── phase 6 — rewrite VCF ────────────────────────────────────────────────────

def rewrite_vcf(
    in_path: Path,
    out_path: Path,
    sig_path: Path,
    annotations: dict[str, tuple[str, str]],
    meta_line_count: int,
):
    """
    Stream through the input VCF, apply REF/ALT annotations, write:
      - out_path  : all records
      - sig_path  : only records with PVALUE_MLOG >= P_SIG_MLOG (p < 5e-8)
    """
    # Build a temp path so we can overwrite in-place safely
    tmp_all = out_path.with_suffix(".vcf.tmp")
    tmp_sig = sig_path.with_suffix(".vcf.tmp")

    total = sig_written = annotated = 0
    not_annotated = 0

    with (
        open(in_path) as fh_in,
        open(tmp_all, "w") as fh_all,
        open(tmp_sig, "w") as fh_sig,
    ):
        meta_seen = 0
        for line in fh_in:
            # ── meta-information lines ──
            if line.startswith("##"):
                fh_all.write(line)
                fh_sig.write(line)
                meta_seen += 1
                continue

            # ── column header ──
            if line.startswith("#CHROM"):
                fh_all.write(line)
                fh_sig.write(line)
                continue

            # ── data line ──
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                fh_all.write(line)
                continue

            vid = cols[2]
            info = cols[7]

            # apply annotation
            if vid in annotations:
                ref, alt = annotations[vid]
                if is_valid_allele(ref):
                    cols[3] = ref
                    annotated += 1
                if is_valid_allele(alt) and cols[4] == ".":
                    cols[4] = alt
            else:
                not_annotated += 1

            new_line = "\t".join(cols) + "\n"
            fh_all.write(new_line)
            total += 1

            # significance filter
            pmlog = parse_pvalue_mlog(info)
            if pmlog is not None and pmlog >= P_SIG_MLOG:
                fh_sig.write(new_line)
                sig_written += 1

    # atomically replace
    tmp_all.replace(out_path)
    tmp_sig.replace(sig_path)

    log.info(
        f"Wrote {total:,} records to {out_path.name}  "
        f"(REF annotated: {annotated:,}, still missing: {not_annotated:,})"
    )
    log.info(
        f"Wrote {sig_written:,} significant records (p<5e-8) to {sig_path.name}"
    )
    return total, sig_written, annotated, not_annotated


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-myvariant", action="store_true")
    parser.add_argument("--skip-ensembl",   action="store_true")
    parser.add_argument("--apply-r-only",   action="store_true",
                        help="Skip API queries; only read r_annotations.tsv and rewrite VCFs")
    args = parser.parse_args()

    # ── scan ──
    log.info("Phase 1: scanning VCF for variants needing annotation …")
    meta_lines, needs = scan_vcf(VCF_IN)
    log.info(f"  meta lines: {meta_lines}  unique rsIDs needing annotation: {len(needs):,}")

    combined: dict[str, tuple[str, str]] = {}

    if not args.apply_r_only:
        # ── myvariant.info ──
        if not args.skip_myvariant:
            mv_annot = query_myvariant(list(needs.keys()))
            combined.update(mv_annot)
        else:
            log.info("Skipping myvariant.info (--skip-myvariant)")

        # ── Ensembl fallback ──
        remaining = [r for r in needs if r not in combined]
        if remaining and not args.skip_ensembl:
            ens_annot = query_ensembl(remaining)
            combined.update(ens_annot)
        else:
            if args.skip_ensembl:
                log.info("Skipping Ensembl REST (--skip-ensembl)")
            elif not remaining:
                log.info("No rsIDs remaining after myvariant.info — skipping Ensembl")

    # ── R / biomaRt annotations (read if present) ──
    r_annot = read_r_annotations(R_ANNOT_TSV)
    combined.update(r_annot)

    # ── write missing list for R ──
    n_missing = write_missing_for_r(needs, set(combined.keys()), MISSING_TSV)
    if n_missing:
        log.info(
            f"  Run annotate_biomart.R to annotate {n_missing:,} remaining variants, "
            f"then re-run this script with --apply-r-only"
        )

    # ── rewrite VCFs ──
    log.info("Phase 6: rewriting VCFs with annotations …")
    rewrite_vcf(VCF_IN, VCF_OUT, VCF_SIG_OUT, combined, meta_lines)

    log.info("Done.")


if __name__ == "__main__":
    main()
