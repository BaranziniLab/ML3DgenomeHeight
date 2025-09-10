library(dplyr)
library(data.table)
library(purrr)

data_dir = "haplotype_sequence"

if(!dir.exists("divergence_mat")){
  dir.create("divergence_mat")
}

data_subdirs = list.files(path = data_dir)

for(data_subdir in data_subdirs){
  
  cat(paste0("Calculating divergence in ", data_subdir, "...\n"))
  
  full_dir = file.path(data_dir, data_subdir)
  
  file_list = list.files(path = full_dir, pattern = "akita_")
  
  if(length(file_list)>0){
    prediction_mat = map(as.list(file_list), function(x){
      filename = unlist(strsplit(x, split = "_"))[2]
      d = fread(file.path(full_dir, x))
      names(d) = filename
      return(d)
    }) %>%
      reduce(cbind.data.frame)
    
    divergence_mat = 1 - cor(prediction_mat, method = "spearman")
    divergence_mat = as.data.frame(divergence_mat)
    fwrite(divergence_mat, file = file.path("divergence_mat", paste0("dmat_", data_subdir, ".csv")), 
           row.names = TRUE, col.names = TRUE)
  }
  

}


