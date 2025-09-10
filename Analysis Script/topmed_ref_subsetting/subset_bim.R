library(dplyr)
library(data.table)

hapmap = fread("HapMap/hapmap.map")
tgp = fread("TGP/1KG_ALL.20190425_NYGC_GATK_hg38_nodups_unrelated.bim")
names(tgp) = c("chr", "id", "distance", "pos", "ref", "alt")
topmed_chr = fread("TopMed/topmed_ref_chr7.bim")
names(topmed_chr) = c("chr", "id", "distance", "pos", "ref", "alt")

