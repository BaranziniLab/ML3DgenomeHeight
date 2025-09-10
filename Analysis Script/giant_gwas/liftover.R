library(rtracklayer)
library(data.table)
library(dplyr)

chain19to38 = import.chain(con ="hg19ToHg38.over.chain")
chain38to19 = import.chain(con ="hg38ToHg19.over.chain")

all = fread("raw_data/all.txt")
afr = fread("raw_data/afr.txt")
eur = fread("raw_data/eur.txt")
amr = fread("raw_data/amr.txt")
eas = fread("raw_data/eas.txt")
sas = fread("raw_data/sas.txt")

convert_build = function(d){
  gr = as(data.frame(chrom=paste0("chr", d$CHR), 
                     start = d$POS, 
                     end = d$POS,
                     SNPID = d$SNPID), 
          "GRanges")
  
  gr = as.data.frame(liftOver(gr, chain19to38))
  
  liftover_d = data.frame(
    HG38_POS = gr$start,
    SNPID = gr$SNPID
  )
  
  d = left_join(d, liftover_d, by = "SNPID")
  d = dplyr::select(d, SNPID, RSID, CHR, POS, HG38_POS,
                    everything()) %>%
    mutate(SNPID = paste(CHR, HG38_POS, OTHER_ALLELE, EFFECT_ALLELE, 
                         sep = ":"))
  return(d)
}

all_hg38 = convert_build(all)
afr_hg38 = convert_build(afr)
eur_hg38 = convert_build(eur)
eas_hg38 = convert_build(eas)
sas_hg38 = convert_build(sas)
amr_hg38 = convert_build(amr)


fwrite(all_hg38, file = "all_hg38.txt", sep = "\t")
fwrite(afr_hg38, file = "afr_hg38.txt", sep = "\t")
fwrite(eur_hg38, file = "eur_hg38.txt", sep = "\t")
fwrite(eas_hg38, file = "eas_hg38.txt", sep = "\t")
fwrite(sas_hg38, file = "sas_hg38.txt", sep = "\t")
fwrite(amr_hg38, file = "amr_hg38.txt", sep = "\t")







