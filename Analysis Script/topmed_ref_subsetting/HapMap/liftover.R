library(rtracklayer)
library(data.table)
library(dplyr)

chain36to38 = import.chain(con ="hg18ToHg38.over.chain")

d = fread("hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map")
names(d) = c("CHR", "SNPID", "DISTANCE", "POS")

d = d %>%
  mutate(CHR = ifelse(CHR == 23, "X", CHR)) %>%
  mutate(CHR = ifelse(CHR == 24, "Y", CHR))

gr = as(data.frame(chrom=paste0("chr", d$CHR), 
                   start = d$POS, 
                   end = d$POS,
                   SNPID = d$SNPID), 
        "GRanges")

gr = as.data.frame(liftOver(gr, chain36to38))

gr = gr %>%
  dplyr::select(seqnames, start, SNPID)

names(gr) = c("CHR", "POS", "SNPID")

fwrite(gr, file = "hapmap.map", sep = "\t")