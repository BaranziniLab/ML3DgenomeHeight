library(dplyr)
library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

make_interaction_map = function(prediction_file, plot_title = "REF"){
  d = unlist(fread(prediction_file))
  names(d) = NULL
  
  d_mat = matrix(NA, nrow = 447, ncol = 447)
  d_mat[lower.tri(d_mat, diag = FALSE)] = d
  
  melted_d_mat = melt(d_mat)
  
  plt = ggplot(melted_d_mat, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(na.rm = TRUE) +
    # geom_text(aes(label = round(value, 2)), vjust = 1) +
    scale_fill_viridis_c(option = "A", na.value = "white") + 
    # scale_fill_gradient(low = "white", high = "blue", na.value = "white") +
    # geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
    labs(x = "Variable", y = "Variable", fill = "Correlation") +
    ggtitle(label = plot_title) + 
    coord_fixed(ratio = 1) +
    theme_void() + 
    theme(legend.position = "none")
  return(plt)
}

data_dir = "top_results"

dirs = list.files(path = data_dir)
subdirs = file.path(data_dir, dirs)

subdir = "top_results/EUR_10_96789820_T_C"
for(subdir in subdirs){
  files = list.files(path = subdir, pattern = "akita_")
  prediction_mat = as.list(files) %>%
    map(function(x){
      seq_name = unlist(strsplit(x, split = "_"))[2]
      d = fread(file.path(subdir, x))
      names(d) = seq_name
      return(d)
    }) %>%
    reduce(cbind.data.frame)
  divergence_matrix = 1 - cor(prediction_mat, method = "spearman")
  
  ref_divergence = divergence_matrix[,which(colnames(divergence_matrix) == "REF")]
  max_ref_divergence = max(ref_divergence)
  seq_diff_from_ref = names(ref_divergence)[which.max(ref_divergence)]
  seq_diff_from_each = rownames(which(divergence_matrix == max(divergence_matrix), arr.ind = TRUE))
  seq_to_plot = c(seq_diff_from_each, seq_diff_from_ref, "REF")
  
  prediction_files = paste0("akita_", seq_to_plot, "_", unlist(strsplit(subdir, split = "[/]"))[2], ".txt")
  prediction_files = file.path(subdir, prediction_files)
}



alt_v1 = make_interaction_map(prediction_files[1], plot_title = "    ALT1")
alt_v2 = make_interaction_map(prediction_files[2], plot_title = "    ALT2")
alt_vm = make_interaction_map(prediction_files[3], plot_title = "    ALTm")
ref_v1 = make_interaction_map(prediction_files[4], plot_title = "    REF")

combined_plt = grid.arrange(alt_v1, alt_v2, alt_vm, ref_v1)

title_text = "Combined Plot Title"
combined_plot = ggdraw(combined_plt) + 
  draw_label("Combined Plot Title",
             fontface = 'bold',
             size = 14, x = 0.5, y = 1, vjust = 0)

combined_plot



f = file_path[10]

snp_name = gsub(pattern = "akita_", replacement = "", file_names[1])
snp_name = gsub(pattern = ".txt", replacement = "", snp_name)

seq_name = strsplit(file_names, split = "_") %>%
  lapply(function(x){x[2]}) %>%
  unlist()

d = as.list(file_names) %>%
  map(function(f){
    seq_name = strsplit(f, split = "_") %>%
      lapply(function(x){x[2]}) %>%
      unlist()
    file_path = file.path(data_dir, f)
    d = unlist(fread(file_path))
    names(d) = NULL
    d = data.frame(pred = d)
    names(d) = seq_name
    return(d)
  }) %>%
  reduce(cbind.data.frame)

divergence_matrix = 1-cor(d, method = "spearman")

ref_divergence = divergence_matrix[,which(colnames(divergence_matrix) == "REF")]
max_divergence = max(ref_divergence)
names(ref_divergence)[which.max(ref_divergence)]

max_divergence = max(divergence_matrix)

which(divergence_matrix == max(divergence_matrix), arr.ind = TRUE)

# divergence_matrix[which.max(divergence_matrix, arr.ind = TRUE)]



pca_result = prcomp(d, scale. = TRUE)
pca_result = as.data.frame(pca_result$rotation[,1:2])
pca_result[["seq_name"]] = rownames(pca_result)

ggplot(data = pca_result, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = seq_name), size = 2, shape = 21)


for(f in file_names){
  snp_name = gsub(pattern = "akita_", replacement = "", f)
  snp_name = gsub(pattern = ".txt", replacement = "", snp_name)
  snp_name = unlist(strsplit(snp_name, split = "_"))[-1] %>%
    paste0(collapse = "_")
  seq_name = strsplit(f, split = "_") %>%
    lapply(function(x){x[2]}) %>%
    unlist()
  file_path = file.path(data_dir, f)
  d = unlist(fread(file_path))
  names(d) = NULL
  
  
  
}



f = file_path[4]
d = unlist(fread(f))
d1 = d
names(d) = NULL

d_mat = matrix(NA, nrow = 447, ncol = 447)
d_mat[lower.tri(d_mat, diag = FALSE)] = d

melted_d_mat = melt(d_mat)

plt1 = ggplot(melted_d_mat, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(na.rm = TRUE) +
  # geom_text(aes(label = round(value, 2)), vjust = 1) +
  scale_fill_viridis_c(option = "A", na.value = "white") + 
  # scale_fill_gradient(low = "white", high = "blue", na.value = "white") +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
  labs(x = "Variable", y = "Variable", fill = "Correlation") +
  coord_fixed(ratio = 1) +
  theme_void() + 
  theme(legend.position = "none"); plt1

grid.newpage()
print(plt1,
      vp = viewport(width = 0.5,
                    height = 0.5,
                    angle = -45))

f = file_path[10]
d = unlist(fread(f))
d2 = d
names(d) = NULL

d_mat = matrix(NA, nrow = 447, ncol = 447)
d_mat[lower.tri(d_mat, diag = FALSE)] = d

melted_d_mat = melt(d_mat)

plt2 = ggplot(melted_d_mat, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(na.rm = TRUE) +
  # geom_text(aes(label = round(value, 2)), vjust = 1) +
  scale_fill_viridis_c(option = "A", na.value = "white") + 
  # scale_fill_gradient(low = "white", high = "blue", na.value = "white") +
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Variable", y = "Variable", fill = "Correlation") +
  theme_void() + 
  theme(legend.position = "none"); plt2

grid.arrange(plt1, plt2, nrow = 1)


d3 = abs(d1 - d2)
d = d3

d_mat = matrix(NA, nrow = 447, ncol = 447)
d_mat[lower.tri(d_mat, diag = FALSE)] = d

melted_d_mat = melt(d_mat)

plt3 = ggplot(melted_d_mat, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(na.rm = TRUE) +
  # geom_text(aes(label = round(value, 2)), vjust = 1) +
  scale_fill_viridis_c(option = "A", na.value = "white") + 
  # scale_fill_gradient(low = "white", high = "blue", na.value = "white") +
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Variable", y = "Variable", fill = "Correlation") +
  theme_void() + 
  theme(legend.position = "top"); plt3
