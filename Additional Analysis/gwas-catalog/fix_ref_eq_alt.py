#!/usr/bin/env python3
"""
fix_ref_eq_alt.py

Some GWAS Catalog risk alleles ARE the reference allele; in such cases the
initial annotation step wrote REF == ALT (e.g., both "A"), which is not a
valid VCF record.

This script:
  1. Scans both VCFs for records where REF == ALT and REF != "."
  2. Re-queries myvariant.info for the corresponding rsIDs to get the
     true dbSNP ALT allele
  3. Rewrites gwas_catalog.vcf and gwas_catalog_significant.vcf with
     corrected ALT values
"""

import re
import sys
import time
import logging
from pathlib import Path

import myvariant

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

HERE = Path(__file__).parent
VCF_FULL = HERE / "gwas_catalog.vcf"
VCF_SIG  = HERE / "gwas_catalog_significant.vcf"

VALID_BASES = set("ACGTacgt")


def is_valid(a: str) -> bool:
    return bool(a) and a != "." and all(c in VALID_BASES for c in a.upper())


def clean(a) -> str:
    if not a or str(a).strip() in ("", ".", "nan", "None"):
        return "."
    s = str(a).strip().upper()
    if "," in s:
        s = s.split(",")[0].strip()
    return s if all(c in VALID_BASES for c in s) else "."


# ── phase 1: collect rsIDs needing fix ───────────────────────────────────────

def collect_bad_ids(path: Path) -> set[str]:
    bad = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 5:
                continue
            ref, alt, vid = c[3], c[4], c[2]
            # Bad = REF is a real base but ALT == REF (risk allele was REF)
            if is_valid(ref) and is_valid(alt) and ref == alt:
                if vid and vid != ".":
                    bad.add(vid)
    return bad


log.info("Collecting rsIDs where REF == ALT …")
bad_ids = collect_bad_ids(VCF_FULL)
log.info(f"  {len(bad_ids):,} unique rsIDs need ALT correction")


# ── phase 2: query myvariant.info for correct ALT ────────────────────────────

mv = myvariant.MyVariantInfo()
alt_lookup: dict[str, str] = {}   # rsid → correct ALT
rsid_list = list(bad_ids)
batch_size = 1000

log.info(f"Querying myvariant.info for {len(rsid_list):,} rsIDs …")
for start in range(0, len(rsid_list), batch_size):
    batch = rsid_list[start : start + batch_size]
    if start % 10000 == 0 and start > 0:
        log.info(f"  progress: {start:,}/{len(rsid_list):,}  fixed: {len(alt_lookup):,}")
    try:
        results = mv.querymany(
            batch,
            scopes="dbsnp.rsid",
            fields="dbsnp.ref,dbsnp.alt,dbsnp.alleles",
            as_dataframe=False,
            verbose=False,
        )
    except Exception as e:
        log.warning(f"  myvariant error at {start}: {e}")
        time.sleep(2)
        continue

    for r in results:
        if r.get("notfound") or "error" in r:
            continue
        query_id = r.get("query", "")
        dbsnp = r.get("dbsnp", {})
        if not dbsnp:
            continue

        ref = clean(dbsnp.get("ref"))
        alt_raw = dbsnp.get("alt")

        # Build a list of known ALT alleles
        alts: list[str] = []
        if isinstance(alt_raw, list):
            alts = [clean(a) for a in alt_raw]
        elif alt_raw:
            alts = [clean(str(alt_raw))]

        # Fallback: parse alleles array
        if not any(is_valid(a) for a in alts):
            alleles = dbsnp.get("alleles")
            if isinstance(alleles, list):
                for a_obj in alleles:
                    a = clean(a_obj.get("allele", ".") if isinstance(a_obj, dict) else str(a_obj))
                    if a != ref and is_valid(a):
                        alts.append(a)
                        break

        # Pick first alt that differs from REF
        chosen_alt = "."
        for a in alts:
            if is_valid(a) and a != ref:
                chosen_alt = a
                break

        if is_valid(chosen_alt) and query_id:
            alt_lookup[query_id] = chosen_alt

log.info(f"myvariant.info: resolved ALT for {len(alt_lookup):,}/{len(bad_ids):,} rsIDs")


# ── phase 3: rewrite both VCFs ───────────────────────────────────────────────

def rewrite(in_path: Path):
    tmp = in_path.with_suffix(".vcf.tmp")
    fixed = skipped = 0
    with open(in_path) as fi, open(tmp, "w") as fo:
        for line in fi:
            if line.startswith("#"):
                fo.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                fo.write(line)
                continue
            ref, alt, vid = cols[3], cols[4], cols[2]
            # Fix only records where REF == ALT AND we have a better ALT
            if is_valid(ref) and is_valid(alt) and ref == alt and vid in alt_lookup:
                cols[4] = alt_lookup[vid]
                fixed += 1
            else:
                skipped += 1
            fo.write("\t".join(cols) + "\n")
    tmp.replace(in_path)
    log.info(f"  {in_path.name}: fixed {fixed:,} records, left unchanged {skipped:,}")


log.info("Rewriting VCFs …")
rewrite(VCF_FULL)
rewrite(VCF_SIG)


# ── phase 4: final QC ─────────────────────────────────────────────────────────

def count_issues(path: Path):
    missing_ref = missing_alt = ref_eq_alt = total = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            total += 1
            ref, alt = c[3], c[4]
            if ref == ".":
                missing_ref += 1
            if alt == ".":
                missing_alt += 1
            if is_valid(ref) and is_valid(alt) and ref == alt:
                ref_eq_alt += 1
    log.info(
        f"  {path.name}: {total:,} records  "
        f"missing_ref={missing_ref:,} ({100*missing_ref/total:.2f}%)  "
        f"missing_alt={missing_alt:,} ({100*missing_alt/total:.2f}%)  "
        f"REF==ALT={ref_eq_alt:,} ({100*ref_eq_alt/total:.2f}%)"
    )


log.info("Final QC:")
count_issues(VCF_FULL)
count_issues(VCF_SIG)
log.info("Done.")
