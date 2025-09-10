library(dplyr)
library(purrr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)

d = fread("max_divergence.csv")

d = d %>%
   arrange(desc(divergence))

t.test(log_divergence ~ group, data = d)

d_haplotype = d %>%
  filter(group == "Haplotypic") %>%
  arrange(snp_name, desc(log_divergence)) %>%
  group_by(snp_name) %>%
  top_n(log_divergence, 1) %>%
  arrange(desc(log_divergence))

lm(data = d_haplotype,  log_divergence ~ abs(beta)) %>%
  summary()

lm(data = filter(d_haplotype, divergence > 0.001), log_divergence ~ abs(beta)) %>%
  summary()
          

plt1 = ggplot(data = d_haplotype, aes(x = abs(beta), y = log_divergence)) + 
  geom_point(size = 2, shape = 21, color = "black", fill = "gray") +
  geom_smooth(method = "lm", color = "black", se = FALSE, formula = 'y ~ x') + 
  xlab("Standardized Haplotypic Effect Size on Height") + 
  ylab(expression(log[10]~"(3D-Genome Divergence)")) +
  theme_bw() + 
  theme(text = element_text(size = 12)); plt1

plt2 = ggplot(data =  filter(d_haplotype, divergence > 0.001), 
              aes(x = abs(beta), y = log_divergence)) + 
  geom_point(size = 2, shape = 21, color = "black", fill = "gray") +
  geom_smooth(method = "lm", color = "black", se = FALSE, formula = 'y ~ x') + 
  xlab("Standardized Haplotypic Effect Size on Height") + 
  ylab(expression(log[10]~"(3D-Genome Divergence)")) +
  theme_bw() + 
  theme(text = element_text(size = 12)); plt2

plt = grid.arrange(plt1, plt2)

ggsave(filename = "figure s1.png", plot = plt, dpi = 1200, 
       device = "png", height = 6, width = 4)
          