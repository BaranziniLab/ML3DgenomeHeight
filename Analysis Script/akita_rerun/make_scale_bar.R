library(ggplot2)

# Create a data frame for the scale bar
scale_data <- data.frame(
  x = seq(-0.9922, 1.5287, length.out = 100),
  y = 1
)

# Plot the scale bar
plt = ggplot(scale_data, aes(x = x, y = y, fill = x)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("#4590c1", "white", "#dc715b"),
    limits = c(-0.9922, 1.5287),
    name = "",
    breaks = c(-0.9922, (1.5287 - 0.9922) / 2, 1.5287),
    labels = round(c(-0.9922, (1.5287 - 0.9922) / 2, 1.5287), 2)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  coord_fixed(ratio = 10)  # Adjust aspect ratio for better visibility


ggsave("HiC_scale_bar.png", plot = plt, width = 6, height = 5, dpi = 1200)
