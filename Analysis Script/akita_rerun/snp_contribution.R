library(dplyr)
library(data.table)
library(purrr)
library(glue)
library(ggplot2)
library(scales)
library(patchwork)

load("artificial_haplotypes.rda")

afr = fread("yengo_result/afr_hg38.txt") %>%
  mutate(snp_name = paste("AFR", SNPID, sep = "_")) %>%
  mutate(snp_name = gsub(pattern = "[:]", replacement = "_", snp_name))
eur = fread("yengo_result/eur_hg38.txt") %>%
  mutate(snp_name = paste("EUR", SNPID, sep = "_")) %>%
  mutate(snp_name = gsub(pattern = "[:]", replacement = "_", snp_name))
amr = fread("yengo_result/amr_hg38.txt") %>%
  mutate(snp_name = paste("AMR", SNPID, sep = "_")) %>%
  mutate(snp_name = gsub(pattern = "[:]", replacement = "_", snp_name))
sas = fread("yengo_result/sas_hg38.txt") %>%
  mutate(snp_name = paste("SAS", SNPID, sep = "_")) %>%
  mutate(snp_name = gsub(pattern = "[:]", replacement = "_", snp_name))
eas = fread("yengo_result/eas_hg38.txt") %>%
  mutate(snp_name = paste("EAS", SNPID, sep = "_")) %>%
  mutate(snp_name = gsub(pattern = "[:]", replacement = "_", snp_name))

gwas = rbind.data.frame(
  afr, eur, amr, sas, eas, eas
) %>%
  distinct(snp_name, .keep_all = TRUE) %>%
  select(snp_name, BETA, SE, N, RSID, EFFECT_ALLELE_FREQ)

names(gwas) = c("snp_name", "beta", "se", "sample_size", "rs_id", "eaf")

dmat_files = file.path("divergence_mat", list.files(path = "divergence_mat")) %>%
  as.list()

vs_ref = map(dmat_files, function(x){
  
  snp_name = gsub(pattern = ".csv", replacement = "",
                  unlist(strsplit(x, split = "dmat_"))[2])
  
  cat(glue("Reading {snp_name} ...\n\n"))
  
  d = fread(x) %>%
    as.data.frame()
  names(d)[1] = "haplotype"
  
  comp_to_ref = select(d, haplotype, REF) %>%
    filter(haplotype != "REF") %>%
    filter(haplotype != "ALTX1")
  
  altx_index = grepl(pattern = "ALTX", comp_to_ref[["haplotype"]])
  alt_index = !altx_index
  
  natural_haplotype_divergence = comp_to_ref[alt_index,]
  artificial_haplotype_divergence = comp_to_ref[altx_index,]
  
  return(list(
    snp_name = snp_name,
    natural_haplotype_divergence = natural_haplotype_divergence,
    artificial_haplotype_divergence = artificial_haplotype_divergence
  ))
})


snp_contribution = map(vs_ref, function(x){
  divergence = x$artificial_haplotype_divergence$REF
  top_snp = max(divergence)
  all_snp = sum(divergence)
  data.frame(
    snp_name = x$snp_name,
    top_snp = top_snp,
    all_snp = all_snp,
    top_percentage = top_snp/all_snp
  )
}) %>%
  bind_rows() %>%
  arrange(desc(top_snp)) %>%
  filter(top_snp >= 0.001)

filter(snp_contribution, top_percentage <= 0.6)

"
             snp_name     top_snp     all_snp top_percentage
1 EUR_4_122914501_C_A 0.007948376 0.013815098      0.5753398
2 EUR_12_57852573_C_A 0.001243476 0.002400103      0.5180929
3  EUR_7_27870446_C_T 0.001132177 0.001998972      0.5663795
"

plot_snp_contribution = function(snp_id1 = "EUR_4_122914501_C_A", 
                                 snp_id2 = "4:122914501:C:A",
                                 plot_title = "SNP Contribution at the chr4:122914501 Locus"){
  top_snp = snp_id1
  top_snp_divergence = vs_ref[[which((map(vs_ref, function(x){
    x$snp_name
  }) %>%
    unlist()) == top_snp)]]$artificial_haplotype_divergence$REF
  
  top_snp_info = artificial_haplotyping_results[[snp_id2]]$fix$ID
  top_snp_info = paste0("chr", top_snp_info)
  
  top_snp_divergence = data.frame(
    snp_id = top_snp_info,
    divergence = top_snp_divergence
  ) %>%
    filter(divergence > 1e-6) %>%
    arrange(desc(divergence))
  
  top_snp_divergence$snp_id = factor(top_snp_divergence$snp_id, 
                                     levels = top_snp_divergence$snp_id)
  
  plt = ggplot(top_snp_divergence, aes(x = snp_id, y = divergence*1e8)) +
    geom_bar(stat = "identity", fill = "gray", color = "black") +
    labs(x = "",
         y = "3D-Genome Divergence") +
    ggtitle(label = plot_title) + 
    scale_y_log10(breaks = c(0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8), 
                  labels = c("1e-8", "1e-7", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1")) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0), text = element_text(size = 12))
  
  return(plt)
}

plt1 = plot_snp_contribution("EUR_4_122914501_C_A", "4:122914501:C:A", 
                             plot_title = "chr4:122914501 Locus")
plt2 = plot_snp_contribution("EUR_12_57852573_C_A", "12:57852573:C:A", 
                             plot_title = "chr12:57852573 Locus")
plt3 = plot_snp_contribution("EUR_7_27870446_C_T", "7:27870446:C:T", 
                             plot_title = "chr7:27870446 Locus")

scatter_plot = ggplot(data = snp_contribution, aes(x = log(top_snp, base = 10), 
                                    y = top_percentage)) + 
  geom_point(size = 3, shape = 21, color = "black", fill = "gray") + 
  xlab(expression(log[10]~"(3D-Genome Divergence)")) +
  ylab("%Contribution of the Top SNP") + 
  scale_y_continuous(labels = percent_format(accuracy = 1)) + 
  theme_bw() + 
  theme(text = element_text(size = 12))

density_plt = ggplot(data = snp_contribution, 
       aes(x = top_percentage)) +
  geom_density(fill = "gray", alpha = 0.5) +
  theme_void() +
  coord_flip() +  # Flip coordinates to align it along the y-axis
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        text = element_text(size = 12))

plt4 = scatter_plot + density_plt +  plot_layout(widths = c(4, 1))

ggsave(filename = "chr4_122914501 Locus SNP Contribution.png", dpi = 1200, device = "png",
       height = 4.5, width = 6, plot = plt1)

ggsave(filename = "chr12_57852573 Locus SNP Contribution.png", dpi = 1200, device = "png",
       height = 4.5, width = 5, plot = plt2)

ggsave(filename = "chr7_27870446 Locus SNP Contribution.png", dpi = 1200, device = "png",
       height = 4.5, width = 5, plot = plt3)

ggsave(filename = "SNP Contribution Density Plot.png", dpi = 1200, device = "png",
       height = 4.5, width = 5, plot = plt4)
