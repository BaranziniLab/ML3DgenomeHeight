library(dplyr)
library(readxl)
library(purrr)
library(ggplot2)
library(ggpubr)

d = read_excel("lcor_chip.xlsx")

as.data.frame(table(unlist(strsplit( paste(filter(d, Targets == "CTCF")$Organ, collapse = ", "),split = ", " )))) %>%
  arrange(desc(Freq))

# Summarize data: Calculate mean and standard deviation for each target
summarized_data = d %>%
  group_by(Targets) %>%
  summarise(
    mean_value = mean(Value),
    sd_value = sd(Value, na.rm = TRUE),
    n = length(Value)
  ) %>%
  filter(n >= 2) %>%
  arrange(desc(n))

summarized_data$Targets = factor(summarized_data$Targets,
                                 levels = summarized_data$Targets)

# Create the bar plot with error bars
plt = ggplot(summarized_data, aes(x = Targets, y = n)) +
  geom_bar(stat = "identity", fill = "lightblue", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Binding Targets", 
       y = "Frequency of Detection",
       title = "chr10:96789820 Locus Binding Pattern") + 
  theme(text = element_text(size = 15)); plt

ggsave(filename = "lcor_binding_pattern.png", dpi = 1200,
       device = "png", width = 5, height = 4.5)

