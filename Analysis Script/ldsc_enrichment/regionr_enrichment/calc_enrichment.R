library(dplyr)
library(data.table)
library(regioneR)
library(GenomicRanges)

divergence = fread("divergent_windows_10AF_variants.txt")
height_gwas = fread("eur_1000_tgp_clumping_result.clumps")

divergence_filtered = divergence |> 
  filter(`1-spearman` > 0.1)

gr_divergence = GRanges(
  seqnames = divergence_filtered$chr,
  ranges = IRanges(
    start = divergence_filtered$pos,
    width = 1  # or set a window size if needed
  )
)

height_gwas_gr = GRanges(
  seqnames = paste0("chr", height_gwas$`#CHROM`),
  ranges = IRanges(
    start = height_gwas$POS - 100000,
    end = height_gwas$POS + 100000
  )
)

universe = data.frame(
  chr = paste0("chr", 1:22),
  length = rep(1.5e8, 22)  # adjust as needed for your genome
) |>
  toGRanges()



hits = findOverlaps(gr_divergence, height_gwas_gr)


set.seed(42)

pt = permTest(
  A = gr_divergence,
  B = height_gwas_gr,
  ntimes = 1000,
  randomize.function = resampleRegions,
  evaluate.function = numOverlaps,
  universe = universe,
  verbose = TRUE
)
