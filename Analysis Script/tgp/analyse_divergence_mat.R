library(dplyr)
library(purrr)
library(data.table)
library(writexl)
library(stringdist)

load("haplotyping_results.rda")
file_list = list.files(path = "divergence_mat", pattern = ".csv")

divergence_list = vector()
for(f in file_list){
  
  file_name = unlist(strsplit(f, split = "[_]"))
  marker_name = paste(file_name[3:6], collapse = ":")
  marker_name = gsub(pattern = ".csv", replacement = "", marker_name)
  
  haplotype_index = ((lapply(haplotyping_results, function(x){
    x[2]
  }) %>%
    unlist()) == marker_name) %>% which()
  haplotype_index = haplotype_index[1]
  
  highlighted_haplotype = haplotyping_results[[haplotype_index]]
  
  remaining_hap_index = highlighted_haplotype[["haplotype_table"]][[tolower(highlighted_haplotype[["ancestry"]])]] >= 10
  highlighted_haplotype[["haplotype_table"]] = highlighted_haplotype[["haplotype_table"]][remaining_hap_index,]
  
  haplotype_count = dim(highlighted_haplotype[["haplotype_table"]])[1]
  
  haplotype_table = data.frame(
    haplotype = highlighted_haplotype[["haplotype_table"]][["haplotype"]],
    haplotype_name = c("REF", paste0("ALT", 1:(haplotype_count-1)))
  )
  
  ancestry = file_name[2]
  start_pos = paste0("chr", file_name[3], "_", as.numeric(file_name[4]) - (2^20)/2)
  
  cat(paste("Loading ", f, "...\n"))
  d = fread(file.path("divergence_mat", f), header = TRUE)[,-1]
  max_divergence = max(as.matrix(d))
  
  divergence_matrix = as.data.frame(d)
  rownames(divergence_matrix) = names(divergence_matrix)
  divergence_matrix = as.matrix(divergence_matrix)
  
  ref_divergence = divergence_matrix[,which(colnames(divergence_matrix) == "REF")]
  max_ref_divergence = max(ref_divergence)
  seq_diff_from_ref = names(ref_divergence)[which.max(ref_divergence)]
  seq_diff_from_each = rownames(which(divergence_matrix == max(divergence_matrix), arr.ind = TRUE))
  
  ref_sequence = filter(haplotype_table, haplotype_name == "REF")$haplotype
  max_ref_sequence = filter(haplotype_table, haplotype_name == seq_diff_from_ref)$haplotype
  max_ref_distance = as.vector(stringdistmatrix(max_ref_sequence, ref_sequence, method = "lv"))
  seq_diff_sequence1 = filter(haplotype_table, haplotype_name == seq_diff_from_each[1])$haplotype
  seq_diff_sequence2 = filter(haplotype_table, haplotype_name == seq_diff_from_each[2])$haplotype
  max_diff_distance = as.vector(stringdistmatrix(seq_diff_sequence1, seq_diff_sequence2, method = "lv"))
  
  divergence_list = rbind.data.frame(
    divergence_list,
    tibble(
      filename = f,
      ancestry = ancestry,
      start_pos = start_pos,
      max_divergence = max_divergence,
      seq_diff_from_ref = seq_diff_from_ref,
      seq_diff_from_each1 = seq_diff_from_each[1],
      seq_diff_from_each2 = seq_diff_from_each[2],
      ref_sequence, max_ref_sequence, max_ref_distance,
      seq_diff_sequence1, seq_diff_sequence2, max_diff_distance
    )
  )
  
}

filtered_list = divergence_list %>%
  filter(max_divergence > 0.001) %>%
  arrange(desc(max_divergence))

plot_command = divergence_list %>%
  filter(max_divergence > 0.01) %>%
  arrange(desc(max_divergence))

command_toref = paste0("map_comparison_no_delta_with_annotations('/Users/wanjun/Desktop/akita_visualize/akita_predictions/akita_", 
       plot_command$seq_diff_from_ref, "_",
       gsub(pattern = ".csv", replacement = ".txt", gsub(pattern = "dmat_", replacement = "", plot_command$filename)),
       "','/Users/wanjun/Desktop/akita_visualize/akita_predictions/akita_", "REF_",
       gsub(pattern = ".csv", replacement = ".txt", gsub(pattern = "dmat_", replacement = "", plot_command$filename)),
       "','", plot_command$start_pos, "','yes', 10, 10,'", paste0(plot_command$start_pos, "_TOREF.pdf')"))

command_toeach = paste0("map_comparison_no_delta_with_annotations('/Users/wanjun/Desktop/akita_visualize/akita_predictions/akita_", 
       plot_command$seq_diff_from_each1, "_",
       gsub(pattern = ".csv", replacement = ".txt", gsub(pattern = "dmat_", replacement = "", plot_command$filename)),
       "','/Users/wanjun/Desktop/akita_visualize/akita_predictions/akita_",  plot_command$seq_diff_from_each2, "_",
       gsub(pattern = ".csv", replacement = ".txt", gsub(pattern = "dmat_", replacement = "", plot_command$filename)),
       "','", plot_command$start_pos, "','yes', 10, 10,'", paste0(plot_command$start_pos, "_TOEACH.pdf')"))
       

write_xlsx(filtered_list, path = "top_divergent_regions_annotated.xlsx")
