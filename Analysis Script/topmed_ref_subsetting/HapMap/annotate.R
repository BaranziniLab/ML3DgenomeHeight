library(rtracklayer)
library(data.table)
library(dplyr)
library(biomaRt)

snpMart = useEnsembl(biomart = "snps", dataset = "hsapiens_snp")

d = fread("hapmap.map")

rs_annotation = getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "allele"),
                       filters = "snp_filter",
                       values = d$SNPID,
                       mart = snpMart)
