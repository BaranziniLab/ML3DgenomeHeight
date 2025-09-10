library(dplyr)
library(data.table)
library(writexl)
library(ggplot2)
library(ggpubr)

lcor_pca = fread("lcor_akita_pca.csv")

plt1 = ggplot(data = lcor_pca, aes(x = pc1, y = pc2)) + 
  geom_point(aes(fill = chr10_96789820_T_C),
             size = 3, color = "black", shape = 21) + 
  labs(fill = "chr10:96789820") + 
  xlab("PC1") + 
  ylab("PC2") + 
  theme_pubr() + 
  theme(text = element_text(size = 15)); plt1

ggsave(filename = "chr10_96789820_T_C_pca12.png", device = "png", dpi = 1200,
       width = 3.5, height = 3.5, plot = plt1)

plt2 = ggplot(data = lcor_pca, aes(x = pc3, y = pc4)) + 
  geom_point(aes(fill = chr10_96789820_T_C),
             size = 3, color = "black", shape = 21) + 
  labs(fill = "chr10:96789820") + 
  xlab("PC3") + 
  ylab("PC4") + 
  theme_pubr() + 
  theme(text = element_text(size = 15)); plt2

ggsave(filename = "chr10_96789820_T_C_pca34.png", device = "png", dpi = 1200,
       width = 3.5, height = 3.5, plot = plt2)


effect_by_pc = cbind.data.frame(
  pc1 = summary(lm(data = lcor_pca, formula = pc1 ~ chr10_96748564_G_A + 
                     chr10_96749834_G_A + chr10_96754057_T_C + chr10_96773208_T_G + 
                     chr10_96774459_C_A + chr10_96774562_T_G + chr10_96789820_T_C))$coefficient[,1],
  pc2 = summary(lm(data = lcor_pca, formula = pc2 ~ chr10_96748564_G_A + 
                     chr10_96749834_G_A + chr10_96754057_T_C + chr10_96773208_T_G + 
                     chr10_96774459_C_A + chr10_96774562_T_G + chr10_96789820_T_C))$coefficient[,1],
  pc3 = summary(lm(data = lcor_pca, formula = pc3 ~ chr10_96748564_G_A + 
                     chr10_96749834_G_A + chr10_96754057_T_C + chr10_96773208_T_G + 
                     chr10_96774459_C_A + chr10_96774562_T_G + chr10_96789820_T_C))$coefficient[,1],
  pc4 = summary(lm(data = lcor_pca, formula = pc4 ~ chr10_96748564_G_A + 
                     chr10_96749834_G_A + chr10_96754057_T_C + chr10_96773208_T_G + 
                     chr10_96774459_C_A + chr10_96774562_T_G + chr10_96789820_T_C))$coefficient[,1],
  pc5 = summary(lm(data = lcor_pca, formula = pc5 ~ chr10_96748564_G_A + 
                     chr10_96749834_G_A + chr10_96754057_T_C + chr10_96773208_T_G + 
                     chr10_96774459_C_A + chr10_96774562_T_G + chr10_96789820_T_C))$coefficient[,1]
) %>%
  # apply(MARGIN = 2, function(x){
  #   (x-mean(x))/sd(x)
  # }) %>%
  as.data.frame()
effect_by_pc[["variant"]] = rownames(effect_by_pc)

write_xlsx(effect_by_pc, path = "lcor_effect_pca.xlsx")
