library(dplyr)
library(data.table)
library(purrr)
library(readxl)
library(writexl)

d = fread("max_divergence.csv")

d = d %>%
  filter(group == "Haplotypic") %>%
  filter(divergence >= 0.01) %>%
  arrange(snp_name, desc(divergence)) %>%
  group_by(snp_name) %>%
  top_n(1, divergence) %>%
  ungroup() %>%
  arrange(desc(divergence))

write_xlsx(d, path = "table1.xlsx")