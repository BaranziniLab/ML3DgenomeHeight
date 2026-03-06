#!/usr/bin/env Rscript
# validate_vcf.R
# Reads gwas_catalog.vcf and produces a quick QC summary.
# Run after make_vcf.py.

suppressPackageStartupMessages({
  library(data.table)
})

vcf_path <- "gwas_catalog.vcf"

cat("Reading VCF:", vcf_path, "\n")

# ── count meta-information lines ─────────────────────────────────────────────
con <- file(vcf_path, "r")
meta_lines <- 0L
repeat {
  ln <- readLines(con, n = 1L, warn = FALSE)
  if (length(ln) == 0L || !startsWith(ln, "##")) break
  meta_lines <- meta_lines + 1L
}
close(con)
cat(sprintf("Meta-information lines (##): %d\n", meta_lines))

# ── load data lines ───────────────────────────────────────────────────────────
vcf <- fread(
  vcf_path,
  sep        = "\t",
  header     = TRUE,
  skip       = meta_lines,
  quote      = "",
  colClasses = "character",
  col.names  = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
)

cat(sprintf("Total variant records    : %s\n", format(nrow(vcf), big.mark = ",")))
cat(sprintf("Unique chromosomes       : %s\n", paste(sort(unique(vcf$CHROM)), collapse = ", ")))
cat(sprintf("Unique rsIDs (non-.)     : %s\n", format(sum(vcf$ID != "."), big.mark = ",")))

# ── parse a few INFO tags for summaries ──────────────────────────────────────
parse_info_tag <- function(info_vec, tag) {
  pattern <- paste0("(?<=", tag, "=)[^;]+")
  regmatches(info_vec, regexpr(pattern, info_vec, perl = TRUE))
}

traits      <- parse_info_tag(vcf$INFO, "TRAIT")
genes       <- parse_info_tag(vcf$INFO, "MAPPED_GENE")
pvals_mlog  <- suppressWarnings(
  as.numeric(parse_info_tag(vcf$INFO, "PVALUE_MLOG"))
)
ancestry    <- parse_info_tag(vcf$INFO, "BROAD_ANCESTRY")

cat(sprintf("Unique traits            : %s\n",
            format(length(unique(traits)), big.mark = ",")))
cat(sprintf("Unique mapped genes      : %s\n",
            format(length(unique(genes[nzchar(genes)])), big.mark = ",")))
cat(sprintf("PVALUE_MLOG range        : %.2f – %.2f (median %.2f)\n",
            min(pvals_mlog, na.rm = TRUE),
            max(pvals_mlog, na.rm = TRUE),
            median(pvals_mlog, na.rm = TRUE)))

# Broad ancestry breakdown
ancestry_all <- unlist(strsplit(ancestry, ","))
ancestry_tab <- sort(table(trimws(ancestry_all)), decreasing = TRUE)
cat("\nBroad ancestry breakdown (top 10):\n")
print(head(ancestry_tab, 10))

# Per-chromosome record counts
chrom_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrMT")
chrom_tab    <- table(factor(vcf$CHROM, levels = chrom_levels))
cat("\nVariant records per chromosome:\n")
print(chrom_tab[chrom_tab > 0])

cat("\nValidation complete.\n")
