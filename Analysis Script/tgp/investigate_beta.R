library(dplyr)
library(data.table)
library(writexl)

load("haplotyping_results.rda")

get_gwas_beta = function(leading_marker){
  
  subset_sumstat = function(sumstat, gene_locus){
    sumstat = sumstat[sumstat$SNPID %in% gene_locus$fix$ID,]
    return(sumstat)
  }
  
  gene_locus = haplotyping_results[[which(names(haplotyping_results) == leading_marker)]]
  
  all = fread("yengo_results_raw/all_hg38.txt") %>%
    subset_sumstat(gene_locus)
  afr = fread("yengo_results_raw/afr_hg38.txt") %>%
    subset_sumstat(gene_locus) %>%
    mutate(AFR_AF = EFFECT_ALLELE_FREQ) %>%
    mutate(AFR_BETA = BETA) %>%
    mutate(AFR_P = P) %>%
    mutate(AFR_N = N)
  
  amr = fread("yengo_results_raw/amr_hg38.txt") %>%
    subset_sumstat(gene_locus) %>%
    mutate(AMR_AF = EFFECT_ALLELE_FREQ) %>%
    mutate(AMR_BETA = BETA) %>%
    mutate(AMR_P = P) %>%
    mutate(AMR_N = N)
  
  eas = fread("yengo_results_raw/eas_hg38.txt") %>%
    subset_sumstat(gene_locus) %>%
    mutate(EAS_AF = EFFECT_ALLELE_FREQ) %>%
    mutate(EAS_BETA = BETA) %>%
    mutate(EAS_P = P) %>%
    mutate(EAS_N = N)
  
  eur = fread("yengo_results_raw/eur_hg38.txt") %>%
    subset_sumstat(gene_locus) %>%
    mutate(EUR_AF = EFFECT_ALLELE_FREQ) %>%
    mutate(EUR_BETA = BETA) %>%
    mutate(EUR_P = P) %>%
    mutate(EUR_N = N)
  sas = fread("yengo_results_raw/sas_hg38.txt") %>%
    subset_sumstat(gene_locus) %>%
    mutate(SAS_AF = EFFECT_ALLELE_FREQ) %>%
    mutate(SAS_BETA = BETA) %>%
    mutate(SAS_P = P) %>%
    mutate(SAS_N = N)
  
  gwas_beta = cbind.data.frame(
    all,
    afr[,13:16],
    amr[,13:16],
    eas[,13:16],
    eur[,13:16],
    sas[,13:16]
  )
  

  return(gwas_beta)
}

lcor_gwas_beta = get_gwas_beta(leading_marker = "10:96789820:T:C")
slc41a2_gwas_beta = get_gwas_beta(leading_marker = "12:104877307:G:T")

# gene_locus = haplotyping_results[[which(names(haplotyping_results) == "12:104877307:G:T")]]

# write_xlsx(slc41a2_gwas_beta, path = "slc41a2_gwas_beta.xlsx")


