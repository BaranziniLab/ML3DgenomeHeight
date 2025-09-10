library(data.table)
library(dplyr)

all = fread("hg38_data/all_hg38.txt")
afr = fread("hg38_data/afr_hg38.txt")
eur = fread("hg38_data/eur_hg38.txt")
amr = fread("hg38_data/amr_hg38.txt")
eas = fread("hg38_data/eas_hg38.txt")
sas = fread("hg38_data/sas_hg38.txt")

make_plink_pvalue_file = function(d, file_name = "all_plink.txt"){
  d = dplyr::select(d, SNPID, CHR, HG38_POS, P)
  names(d) = c("SNP", "CHR", "BP", "P")
  d$P = as.numeric(d$P)
  d = na.omit(d)
  fwrite(d, file = file_name, sep = "\t")
}

make_plink_pvalue_file(all, file_name = "plink_data/all_plink.txt")
make_plink_pvalue_file(afr, file_name = "plink_data/afr_plink.txt")
make_plink_pvalue_file(eur, file_name = "plink_data/eur_plink.txt")
make_plink_pvalue_file(eas, file_name = "plink_data/eas_plink.txt")
make_plink_pvalue_file(sas, file_name = "plink_data/sas_plink.txt")
make_plink_pvalue_file(amr, file_name = "plink_data/amr_plink.txt")
