library(dplyr)
library(data.table)

afr_glm = fread("plink_result/height_afr_result.height.glm.linear")
amr_glm = fread("plink_result/height_amr_result.height.glm.linear")
eur_glm = fread("plink_result/height_eur_result.height.glm.linear")
eas_glm = fread("plink_result/height_eas_result.height.glm.linear")

afr_clump = fread("clumped_result/afr_1000_tgp_clumping_result.clumps")
amr_clump = fread("clumped_result/amr_1000_tgp_clumping_result.clumps")
eur_clump = fread("clumped_result/eur_1000_tgp_clumping_result.clumps")
eas_clump = fread("clumped_result/eas_893_tgp_clumping_result.clumps")

afr_replicated = afr_clump %>%
  select(ID, SP2, P) %>%
  rename("P_ORG" = "P") %>%
  left_join(afr_glm, by = "ID") %>%
  arrange(P) %>%
  mutate(FDR = p.adjust(P, method = "fdr")) %>%
  filter(FDR <= 0.05) %>%
  rename("CHR" = "#CHROM") %>%
  select(CHR, POS, REF, ALT, A1_FREQ, BETA, SE, P_ORG, P, FDR, SP2)

amr_replicated = amr_clump %>%
  select(ID, SP2, P) %>%
  rename("P_ORG" = "P") %>%
  left_join(amr_glm, by = "ID") %>%
  arrange(P) %>%
  mutate(FDR = p.adjust(P, method = "fdr")) %>%
  filter(FDR <= 0.05) %>%
  rename("CHR" = "#CHROM") %>%
  select(CHR, POS, REF, ALT, A1_FREQ, BETA, SE, P_ORG, P, FDR, SP2)

eur_replicated = eur_clump %>%
  select(ID, SP2, P) %>%
  rename("P_ORG" = "P") %>%
  left_join(eur_glm, by = "ID") %>%
  arrange(P) %>%
  mutate(FDR = p.adjust(P, method = "fdr")) %>%
  filter(FDR <= 0.05) %>%
  rename("CHR" = "#CHROM") %>%
  select(CHR, POS, REF, ALT, A1_FREQ, BETA, SE, P_ORG, P, FDR, SP2)

eas_replicated = eas_clump %>%
  select(ID, SP2, P) %>%
  rename("P_ORG" = "P") %>%
  left_join(eas_glm, by = "ID") %>%
  arrange(P) %>%
  mutate(FDR = p.adjust(P, method = "fdr")) %>%
  filter(FDR <= 0.05) %>%
  rename("CHR" = "#CHROM") %>%
  select(CHR, POS, REF, ALT, A1_FREQ, BETA, SE, P_ORG, P, FDR, SP2)

fwrite(afr_replicated, file = "afr_results.tsv", sep = "\t")
fwrite(eur_replicated, file = "eur_results.tsv", sep = "\t")
fwrite(amr_replicated, file = "amr_results.tsv", sep = "\t")
