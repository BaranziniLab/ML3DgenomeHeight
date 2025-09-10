library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

afr_topmed_result = fread("plink_result/height_afr_result.height.glm.linear")
afr_yengo_result = fread("yengo_result/afr_hg38.txt")

eur_topmed_result = fread("plink_result/height_eur_result.height.glm.linear")
eur_yengo_result = fread("yengo_result/eur_hg38.txt")

amr_topmed_result = fread("plink_result/height_amr_result.height.glm.linear")
amr_yengo_result = fread("yengo_result/amr_hg38.txt")

eas_topmed_result = fread("plink_result/height_eas_result.height.glm.linear")
eas_yengo_result = fread("yengo_result/eas_hg38.txt")

make_nlogp_fig = function(topmed_result, yengo_result, plot_title = "AFR Results"){
  topmed_result = topmed_result %>%
    select(ID, BETA, P) %>%
    mutate(BETA = as.numeric(BETA)) %>%
    mutate(P = as.numeric(P))
  names(topmed_result) = c("id", "topmed_beta", "topmed_p")
  
  yengo_result = yengo_result %>%
    select(SNPID, BETA, P) %>%
    mutate(BETA = as.numeric(BETA)) %>%
    mutate(P = as.numeric(P))
  names(yengo_result) = c("id", "yengo_beta", "yengo_p")
  
  combined_result = inner_join(topmed_result, yengo_result, by = "id") %>%
    mutate(topmed_nlogp = -log10(topmed_p)) %>%
    mutate(yengo_nlogp = -log10(yengo_p))
  
  plt = ggplot(data = combined_result, aes(x = topmed_nlogp, y = yengo_nlogp)) + 
    geom_point() + 
    ggtitle(label = plot_title) + 
    xlab("-log(P) from TopMed") + 
    ylab("-log(P) from Yengo et al.") +
    theme_pubr() + 
    theme(text = element_text(size = 15))
  
  ggsave(filename = paste0(plot_title, ".png"), dpi = 1200, device = "png",
         plot = plt, height = 6, width = 6)
  
  return(plt)
}


make_beta_fig = function(topmed_result, yengo_result, plot_title = "AFR Effect Size"){
  topmed_result = topmed_result %>%
    select(ID, BETA, P) %>%
    mutate(BETA = as.numeric(BETA)) %>%
    mutate(P = as.numeric(P)) %>%
    na.omit()
  topmed_result$BETA = (topmed_result$BETA - mean(topmed_result$BETA))/sd(topmed_result$BETA)
  
  names(topmed_result) = c("id", "topmed_beta", "topmed_p")
  
  yengo_result = yengo_result %>%
    select(SNPID, BETA, P) %>%
    mutate(BETA = as.numeric(BETA)) %>%
    mutate(P = as.numeric(P)) %>%
    na.omit()
  yengo_result$BETA = (yengo_result$BETA - mean(yengo_result$BETA))/sd(yengo_result$BETA)
  
  names(yengo_result) = c("id", "yengo_beta", "yengo_p")
  
  combined_result = inner_join(topmed_result, yengo_result, by = "id")
  
  plt = ggplot(data = combined_result, aes(x = topmed_beta, y = yengo_beta)) + 
    geom_point() + 
    ggtitle(label = plot_title) + 
    xlab("Effect size from TopMed") + 
    ylab("Effect size from Yengo et al.") +
    theme_pubr() + 
    theme(text = element_text(size = 15))
  
  ggsave(filename = paste0(plot_title, ".png"), dpi = 1200, device = "png",
         plot = plt, height = 6, width = 6)
  
  return(plt)
}

plt_afr = make_nlogp_fig(afr_topmed_result, afr_yengo_result, "AFR Results")
plt_eur = make_nlogp_fig(eur_topmed_result, eur_yengo_result, "EUR Results")
plt_amr = make_nlogp_fig(amr_topmed_result, amr_yengo_result, "AMR Results")
plt_eas = make_nlogp_fig(eas_topmed_result, eas_yengo_result, "EAS Results")


plt_afr = make_beta_fig(afr_topmed_result, afr_yengo_result, "AFR Effect Size")
plt_eur = make_beta_fig(eur_topmed_result, eur_yengo_result, "EUR Effect Size")
plt_amr = make_beta_fig(amr_topmed_result, amr_yengo_result, "AMR Effect Size")
plt_eas = make_beta_fig(eas_topmed_result, eas_yengo_result, "EAS Effect Size")
