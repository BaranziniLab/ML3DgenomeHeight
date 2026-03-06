#!/usr/bin/env Rscript
# annotate_biomart.R
#
# Fallback annotation for rsIDs that myvariant.info + Ensembl REST could not resolve.
# Reads:   missing_annotations.tsv  (written by annotate_alleles.py)
# Queries: Ensembl BioMart (hsapiens_snp) via biomaRt
# Writes:  r_annotations.tsv  → read back by annotate_alleles.py --apply-r-only
#
# Run AFTER annotate_alleles.py (first pass without --apply-r-only), then:
#   uv run annotate_alleles.py --apply-r-only
#
# Requirements: biomaRt (already installed)

suppressPackageStartupMessages({
  library(biomaRt)
  library(data.table)
})

missing_path <- "missing_annotations.tsv"
out_path     <- "r_annotations.tsv"

if (!file.exists(missing_path)) {
  cat("No missing_annotations.tsv found — nothing to do.\n")
  quit(status = 0)
}

missing <- fread(missing_path, colClasses = "character")
cat(sprintf("Variants to annotate via biomaRt: %s\n",
            format(nrow(missing), big.mark = ",")))

if (nrow(missing) == 0L) {
  cat("File is empty — writing empty r_annotations.tsv\n")
  fwrite(data.table(rsid = character(), ref = character(), alt = character()),
         out_path, sep = "\t")
  quit(status = 0)
}

rsids <- missing$rsid

# ── connect to Ensembl SNP mart ───────────────────────────────────────────────
cat("Connecting to Ensembl SNP mart (hsapiens_snp) …\n")
snp_mart <- tryCatch(
  useMart("ENSEMBL_MART_SNP",
          dataset  = "hsapiens_snp",
          host     = "https://www.ensembl.org"),
  error = function(e) {
    cat("WARNING: could not connect to Ensembl BioMart:", conditionMessage(e), "\n")
    NULL
  }
)

if (is.null(snp_mart)) {
  cat("BioMart unavailable — falling back to Ensembl REST API for remaining variants\n")
  # Use Ensembl REST API in R as last resort
  library(httr)
  library(jsonlite)

  ENSEMBL_POST <- "https://rest.ensembl.org/variation/homo_sapiens"
  batch_size   <- 200L
  results      <- list()

  for (start in seq(1, length(rsids), by = batch_size)) {
    batch <- rsids[start:min(start + batch_size - 1L, length(rsids))]
    tryCatch({
      resp <- POST(
        ENSEMBL_POST,
        add_headers("Content-Type" = "application/json",
                    "Accept"       = "application/json"),
        body   = toJSON(list(ids = batch), auto_unbox = FALSE),
        encode = "raw",
        timeout(30)
      )
      if (status_code(resp) == 200) {
        data <- fromJSON(content(resp, "text", encoding = "UTF-8"),
                         simplifyDataFrame = FALSE)
        for (rsid in names(data)) {
          for (mapping in data[[rsid]][["mappings"]]) {
            allele_str <- mapping[["allele_string"]]
            if (!is.null(allele_str) && grepl("/", allele_str)) {
              parts <- strsplit(allele_str, "/")[[1]]
              if (nchar(parts[1]) > 0) {
                results[[rsid]] <- list(ref = parts[1],
                                        alt = if (length(parts) > 1) parts[2] else ".")
                break
              }
            }
          }
        }
      }
    }, error = function(e) invisible(NULL))
    Sys.sleep(0.07)  # Ensembl rate limit
    if ((start %% 2000) == 1 && start > 1)
      cat(sprintf("  REST progress: %d/%d\n", start, length(rsids)))
  }

  out <- rbindlist(lapply(names(results), function(id) {
    data.table(rsid = id, ref = results[[id]]$ref, alt = results[[id]]$alt)
  }))
  cat(sprintf("REST fallback annotated %d/%d\n", nrow(out), length(rsids)))
  fwrite(out, out_path, sep = "\t")
  cat("Written:", out_path, "\n")
  quit(status = 0)
}

# ── BioMart query ─────────────────────────────────────────────────────────────
# BioMart returns allele info; we query in chunks to avoid timeouts
batch_size <- 2000L
all_results <- vector("list", ceiling(length(rsids) / batch_size))

cat(sprintf("Querying BioMart in batches of %d …\n", batch_size))
for (i in seq_along(all_results)) {
  start <- (i - 1L) * batch_size + 1L
  end   <- min(i * batch_size, length(rsids))
  batch <- rsids[start:end]

  cat(sprintf("  batch %d/%d  (%d–%d)\n",
              i, length(all_results), start, end))

  res <- tryCatch(
    getBM(
      attributes = c("refsnp_id", "allele"),
      filters    = "snp_filter",
      values     = batch,
      mart       = snp_mart
    ),
    error = function(e) {
      cat("  WARNING: BioMart batch failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  all_results[[i]] <- res
  Sys.sleep(0.5)   # polite pause between BioMart requests
}

bm_df <- rbindlist(Filter(Negate(is.null), all_results))
cat(sprintf("BioMart returned %s rows for %s rsIDs\n",
            format(nrow(bm_df), big.mark = ","),
            format(length(rsids), big.mark = ",")))

if (nrow(bm_df) == 0L) {
  cat("No results from BioMart — writing empty r_annotations.tsv\n")
  fwrite(data.table(rsid = character(), ref = character(), alt = character()),
         out_path, sep = "\t")
  quit(status = 0)
}

# BioMart 'allele' field: "A/T" or "A/T/C" — first token is REF
bm_df[, c("ref", "alt") := {
  parts <- strsplit(allele, "/")
  list(
    sapply(parts, `[`, 1L),
    sapply(parts, function(p) if (length(p) > 1L) p[2L] else ".")
  )
}]

# Keep one entry per rsID (take the first usable REF)
valid_bases <- "^[ACGTacgt]+$"
bm_clean <- bm_df[grepl(valid_bases, ref)][
  !duplicated(refsnp_id)
][, .(rsid = refsnp_id, ref, alt)]

cat(sprintf("After filtering: %s annotated rsIDs\n",
            format(nrow(bm_clean), big.mark = ",")))
fwrite(bm_clean, out_path, sep = "\t")
cat("Written:", out_path, "\n")
