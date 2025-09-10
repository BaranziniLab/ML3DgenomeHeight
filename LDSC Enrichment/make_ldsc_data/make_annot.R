library(data.table)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

# Loading AKITA 3d-genome-altering variants
akita_variants = fread("akita_variants/variants_in_divergent_windows_10AF.txt")
names(akita_variants) = c("chr", "pos", "ref", "alt", "cm")
akita_variants = akita_variants[,-5]
akita_variants$chr = as.numeric(gsub(pattern = "chr", replacement = "", akita_variants$chr))

gr_hg38 = GRanges(
  seqnames = paste0("chr", akita_variants$chr),
  ranges = IRanges(start = akita_variants$pos, width = 1),
  ref = akita_variants$ref,
  alt = akita_variants$alt
)

chain_file = file.path("liftover_chain", "hg38ToHg19.over.chain")
chain = import.chain(chain_file)
gr_hg19 =liftOver(gr_hg38, chain)
gr_hg19_unlisted = unlist(gr_hg19)

akita_variants_hg19 = data.frame(
  chr = as.numeric(gsub("chr", "", seqnames(gr_hg19_unlisted))),
  pos = start(gr_hg19_unlisted),
  genome_3d_altering = 1
)

# Generate annot for all 22 chromosomes
for(chromosome in 1:22){
  cat(paste0("Annotating chr", chromosome, "...\n"))
  chr_path = file.path("ld_annot_template", paste0("nean.", chromosome, ".annot"))
  chrx = fread(chr_path)
  
  chrx = left_join(chrx, akita_variants_hg19, by = join_by("CHR" == "chr", "BP" == "pos"))
  chrx$genome_3d_altering = ifelse(is.na(chrx$genome_3d_altering), 0, chrx$genome_3d_altering)
  chrx = chrx |>
    dplyr::select(CHR, BP, SNP, CM, base, Coding_UCSC, genome_3d_altering) |>
    dplyr::distinct(BP, .keep_all = TRUE)
  fwrite(chrx, file = file.path("ld_annot", paste0("genome_3d.", chromosome, ".annot")), sep = "\t")
}

