library(data.table)
library(dplyr)

afr = fread("gwas_raw/afr_hg38.txt")
eur = fread("gwas_raw/eur_hg38.txt")
eas = fread("gwas_raw/eas_hg38.txt")
sas = fread("gwas_raw/sas_hg38.txt")
amr = fread("gwas_raw/amr_hg38.txt")
all = fread("gwas_raw/all_hg38.txt")
hm3 = unlist(fread("hm_snp/w_hm3.snp", header = FALSE))
names(hm3) = NULL
hm3 = data.frame(SNP = hm3)

munge_gwas = function(gwas = afr, gwas_name = "afr"){
  gwas = gwas |>
    select(RSID, EFFECT_ALLELE, OTHER_ALLELE, BETA, N)
  names(gwas) = c("SNP", "A1", "A2", "Z", "N")
  gwas = left_join(hm3, gwas, by = "SNP")
  fwrite(gwas, file = file.path("gwas_munged", paste0(gwas_name, "_height.sumstats")), sep = "\t")
}

munge_gwas(afr, gwas_name = "afr")
munge_gwas(eur, gwas_name = "eur")
munge_gwas(eas, gwas_name = "eas")
munge_gwas(sas, gwas_name = "sas")
munge_gwas(amr, gwas_name = "amr")
munge_gwas(all, gwas_name = "all")





