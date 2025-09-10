library(dplyr)
library(data.table)
library(vcfR)

afr_fam = fread("tgp_fam/1KG_AFR.20190425_NYGC_GATK_hg38_nodups.fam")
amr_fam = fread("tgp_fam/1KG_AMR.20190425_NYGC_GATK_hg38_nodups.fam")
eur_fam = fread("tgp_fam/1KG_EUR.20190425_NYGC_GATK_hg38_nodups.fam")
eas_fam = fread("tgp_fam/1KG_EAS.20190425_NYGC_GATK_hg38_nodups.fam")
sas_fam = fread("tgp_fam/1KG_SAS.20190425_NYGC_GATK_hg38_nodups.fam")

afr_results = fread("yengo_results/afr_1000_tgp_clumping_result.clumps") %>%
  mutate(ANCESTRY = "AFR")
amr_results = fread("yengo_results/amr_1000_tgp_clumping_result.clumps") %>%
  mutate(ANCESTRY = "AMR")
eur_results = fread("yengo_results/eur_1000_tgp_clumping_result.clumps") %>%
  mutate(ANCESTRY = "EUR") 
eas_results = fread("yengo_results/eas_893_tgp_clumping_result.clumps") %>%
  mutate(ANCESTRY = "EAS")
sas_results = fread("yengo_results/sas_103_tgp_clumping_result.clumps") %>%
  mutate(ANCESTRY = "SAS")

yengo_results = rbind.data.frame(afr_results,
                                 amr_results, 
                                 eur_results, 
                                 eas_results,
                                 sas_results) %>%
  arrange(`#CHROM`, P)

haplotyping_results = as.list(yengo_results[["ID"]])
names(haplotyping_results) = yengo_results[["ID"]]
for(i in 1:dim(yengo_results)[1]){
  try({
    chromosome = yengo_results[i,][["#CHROM"]]
    p_value = yengo_results[i,][["P"]]
    ancestry = yengo_results[i,][["ANCESTRY"]]
    leading_marker = yengo_results[i,][["ID"]]
    tagging_markers = unlist(strsplit(yengo_results[i,][["SP2"]], split = ","))
    all_markers = c(leading_marker, tagging_markers)
    
    system(paste0(
      'bcftools-1.18/bcftools view -r ',
      paste0(chromosome,":",
             strsplit(all_markers, split = "[:]") %>%
               lapply(function(x){x[2]}) %>%
               unlist() %>%
               as.numeric() %>% 
               range() %>%
               paste(collapse = "-")
      ),
      " -o temp.vcf.gz -O z ",
      gsub(pattern = "XXX", replacement = chromosome, 
           "tgp_vcf/ALL.chrXXX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")
    ))
    
    phased_vcf = read.vcfR(file = "temp.vcf.gz")
    system("rm temp.vcf.gz")
    
    fix = as.data.frame(phased_vcf@fix) %>%
      mutate(ID = paste(CHROM, POS, REF, ALT, sep = ":")) %>%
      select(CHROM, POS, ID, REF, ALT)
    gt = as.data.frame(phased_vcf@gt)[,-1]
    subjects = names(gt)
    marker_index = which(fix$ID %in% all_markers)
    fix = fix[marker_index,]
    gt = gt[marker_index,]
    
    afr_index = which(subjects %in% unlist(afr_fam[,1]))
    eur_index = which(subjects %in% unlist(eur_fam[,1]))
    amr_index = which(subjects %in% unlist(amr_fam[,1]))
    eas_index = which(subjects %in% unlist(eas_fam[,1]))
    sas_index = which(subjects %in% unlist(sas_fam[,1]))
    
    for(j in 1:dim(fix)[1]){
      ref_allele = fix[["REF"]][j]
      alt_allele = fix[["ALT"]][j]
      
      gt[j,] = gsub(pattern = "0", replacement = ref_allele, gt[j,]) %>%
        gsub(pattern = "1", replacement = alt_allele)
    }
    
    gt_afr = gt[,afr_index]
    gt_amr = gt[,amr_index]
    gt_eur = gt[,eur_index]
    gt_eas = gt[,eas_index]
    gt_sas = gt[,sas_index]
    
    get_maternal = function(x){
      x = apply(x, MARGIN = c(1,2), 
                function(x){
                  unlist(strsplit(x, split = ""))[1]
                })
      colnames(x) = paste0(colnames(x), "_maternal")
      x = as.data.frame(x)
      return(x)
    }
    
    get_paternal = function(x){
      x = apply(x, MARGIN = c(1,2), 
                function(x){
                  unlist(strsplit(x, split = ""))[3]
                })
      colnames(x) = paste0(colnames(x), "_paternal")
      x = as.data.frame(x)
      return(x)
    }
    
    gt_afr_haplotype = cbind.data.frame(get_maternal(gt_afr),
                                        get_paternal(gt_afr)) %>%
      apply(MARGIN = 2, function(x){
        paste(x, collapse = "")
      })
    gt_amr_haplotype = cbind.data.frame(get_maternal(gt_amr),
                                        get_paternal(gt_amr)) %>%
      apply(MARGIN = 2, function(x){
        paste(x, collapse = "")
      })
    gt_eur_haplotype = cbind.data.frame(get_maternal(gt_eur),
                                        get_paternal(gt_eur)) %>%
      apply(MARGIN = 2, function(x){
        paste(x, collapse = "")
      })
    gt_eas_haplotype = cbind.data.frame(get_maternal(gt_eas),
                                        get_paternal(gt_eas)) %>%
      apply(MARGIN = 2, function(x){
        paste(x, collapse = "")
      })
    gt_sas_haplotype = cbind.data.frame(get_maternal(gt_sas),
                                        get_paternal(gt_sas)) %>%
      apply(MARGIN = 2, function(x){
        paste(x, collapse = "")
      })
    
    make_haplotype_table = function(x, ancestry){
      x = table(x)
      x = tibble(
        haplotype = names(x),
        count = as.vector(x)
      )
      names(x)[2] = ancestry
      return(x)
    }
    
    haplotype_table = make_haplotype_table(gt_afr_haplotype, "afr") %>%
      full_join(make_haplotype_table(gt_amr_haplotype, "amr"), by = "haplotype") %>%
      full_join(make_haplotype_table(gt_eur_haplotype, "eur"), by = "haplotype") %>%
      full_join(make_haplotype_table(gt_eas_haplotype, "eas"), by = "haplotype") %>%
      full_join(make_haplotype_table(gt_sas_haplotype, "sas"), by = "haplotype")
    haplotype_table[is.na(haplotype_table)] = 0
    
    haplotyping_results[[i]] = list(
      ancestry = ancestry,
      leading_marker = leading_marker,
      p_value = p_value,
      fix = fix,
      haplotype_table = haplotype_table
    )
  }, silent = TRUE)
}

save(haplotyping_results, file = "haplotyping_results.rda")

# fasta_sequences = readDNAStringSet("ref_sequence/hg38_reference.fa", use.names = FALSE)



