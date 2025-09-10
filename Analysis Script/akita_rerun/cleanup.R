library(dplyr)
library(data.table)
library(purrr)
library(glue)
library(ggplot2)

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

natural_haplotype_divergence = map(vs_ref, function(x){
  divergence = x$natural_haplotype_divergence$REF
  divergence = divergence[divergence>1e-10]
  return(divergence)
}) %>% 
  unlist() %>%
  log(base = 10)

max_natural_divergence = map(vs_ref, function(x){
  divergence = max(x$natural_haplotype_divergence$REF)
  haplotype = x$natural_haplotype_divergence$haplotype[which.max(x$natural_haplotype_divergence$REF)]
  snp_name = x$snp_name
  
  divergence = data.frame(snp_name, haplotype, divergence)
  divergence = divergence %>%
    filter(divergence>1e-10) %>%
    mutate(log_divergence = log(divergence, base = 10))
  return(divergence)
}) %>% 
  bind_rows()

max_artificial_divergence = map(vs_ref, function(x){
  divergence = max(x$artificial_haplotype_divergence$REF)
  haplotype = x$artificial_haplotype_divergence$haplotype[which.max(x$artificial_haplotype_divergence$REF)]
  snp_name = x$snp_name
  
  divergence = data.frame(snp_name, haplotype, divergence)
  divergence = divergence %>%
    filter(divergence>1e-10) %>%
    mutate(log_divergence = log(divergence, base = 10))
  return(divergence)
}) %>% 
  bind_rows()


max_divergence = rbind.data.frame(
  max_natural_divergence,
  max_artificial_divergence
) %>%
  cbind.data.frame(
    group = c(rep("Haplotypic", dim(max_natural_divergence)[1]), 
              rep("Single-SNP", dim(max_artificial_divergence)[1]))
  )

p_values = map(artificial_haplotyping_results, function(x){
  if(length(x)==5){
    ancestry = x[[1]]
    marker = x[[2]]
    linker = paste(ancestry, marker, sep = "_")
    linker = gsub(pattern = "[:]", replacement = "_", linker)
    p_value = x[[3]]
    return(data.frame(snp_name = linker, p_value = p_value))
  }else{
    return(NULL)
  }
}) %>%
  bind_rows()

max_divergence = max_divergence %>% 
  left_join(p_values, by = "snp_name") %>%
  left_join(gwas, by = "snp_name")

write.csv(max_divergence, file = "max_divergence.csv",
          quote = FALSE, row.names = FALSE)

plt1 = ggplot(max_divergence, aes(x = log_divergence, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.1, color = "black") +
  labs(x = expression(log[10]~"(3D-Genome Divergence)"), 
       y = "Frequency", 
       title = "Haplotypic and single-SNP 3D-genome effects",
       fill = "") +
  scale_fill_manual(values = c("Single-SNP" = "#0072B2",
                               "Haplotypic" = "#D55E00")) +
  theme_minimal() + 
  theme(text = element_text(size = 12))

ggsave(filename = "3d Genome Effects.png", dpi = 1200, device = "png",
       height = 4, width = 5.5, plot = plt1)


natural_divergence = filter(max_divergence, group == "Haplotypic")
artificial_divergence = filter(max_divergence, group == "Single-SNP")

png("qqplot_natural_divergence.png", 
    width = 4, height = 4, units = "in", res = 1200)
qqnorm(natural_divergence$log_divergence, 
       main = expression("Q-Q Plot for 3D-genome Divergence"),
       ylab = expression(log[10]~"(3D-Genome Divergence)"), 
       bg = "#E0703080", pch = 21, cex = 1, ylim = c(-8.3, -0.5))
qqline(natural_divergence$log_divergence, col = "black", lty = 1, lwd = 1)
dev.off()

png("qqplot_artificial_divergence.png", 
    width = 4, height = 4, units = "in", res = 1200)
qqnorm(artificial_divergence$log_divergence, 
       main = expression("Q-Q Plot for 3D-genome Divergence"),
       ylab = expression(log[10]~"(3D-Genome Divergence)"), 
       bg = "#0072B280", pch = 21, cex = 1, ylim = c(-8.3, -0.5))
qqline(artificial_divergence$log_divergence, col = "black", lty = 1, lwd = 1)
dev.off()

calc_lambda = function(observed_values){
  observed_chi_squared = qnorm(ppoints(length(observed_values)))^2
  observed_median = median(observed_chi_squared, na.rm = TRUE)
  expected_median = qchisq(0.5, df = 1)
  lambda = observed_median / expected_median
  return(lambda)
}

calc_lambda(artificial_divergence$log_divergence)
calc_lambda(natural_divergence$log_divergence)

top_snp = "EUR_10_96789820_T_C"
top_snp_divergence = vs_ref[[which((map(vs_ref, function(x){
  x$snp_name
}) %>%
  unlist()) == top_snp)]]$artificial_haplotype_divergence$REF

top_snp_info = artificial_haplotyping_results[["10:96789820:T:C"]]$fix$ID
top_snp_info = paste0("chr", top_snp_info)

top_snp_divergence = data.frame(
  snp_id = top_snp_info,
  divergence = top_snp_divergence
) %>%
  filter(divergence > 1e-8) %>%
  arrange(desc(divergence))

top_snp_divergence$snp_id = factor(top_snp_divergence$snp_id, 
                                   levels = top_snp_divergence$snp_id)

plt2 = ggplot(top_snp_divergence, aes(x = snp_id, y = divergence*1e8)) +
  geom_bar(stat = "identity", fill = "gray", color = "black") +
  labs(x = "",
       y = "3D-Genome Divergence") +
  ggtitle(label = "chr10:96789820 Locus") + 
  scale_y_log10(breaks = c(0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8), 
                labels = c("1e-8", "1e-7", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0), text = element_text(size = 12))

ggsave(filename = "LCOR Locus SNP Contribution.png", dpi = 1200, device = "png",
       height = 4.5, width = 4, plot = plt2)


