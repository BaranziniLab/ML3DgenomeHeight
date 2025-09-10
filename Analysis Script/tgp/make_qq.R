library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

file_list = list.files(path = "divergence_mat", pattern = ".csv")

divergence_list = vector()
for(f in file_list){
  
  cat(paste("Loading ", f, "...\n"))
  d = fread(file.path("divergence_mat", f), header = TRUE)[,-1]
  max_divergence = max(as.matrix(d))
  
  divergence_list = rbind.data.frame(
    divergence_list,
    tibble(
      filename = f,
      max_divergence = max_divergence
    )
  )
  
}


# Assuming max_divergence is your vector of numbers
max_divergence = divergence_list$max_divergence
# max_divergence = -log(max_divergence, base = 10)
max_divergence = max_divergence[is.finite(max_divergence)]

pdf(file = "qqplot.pdf", width = 4.5, height = 5)
# Create a QQ plot

# qqnorm(max_divergence, main = expression("Q-Q Plot for "-log((1-rho)[max])),
qqnorm(max_divergence, main = expression("Q-Q Plot for "(1-rho)[max]),
       xlab = "Theoretical Quantiles",
       ylab = "Maximum Divergence")
qqline(max_divergence, col = "#F8766D", lwd = 2, lty = 2)  # Add a reference line
legend("topleft", legend = "Reference Line", col = "#F8766D", lwd = 2, lty = 2)
dev.off()
