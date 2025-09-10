source("make_sequence.R")

for(i in 1:5335){
  cat(paste0("Calculating hit #", i, "...\n"))
  hit = haplotyping_results[[i]]
  sequence_library = hit_to_sequence(hit)
  if(!is.null(sequence_library)){
    try({
      write_sequence(sequence_library = sequence_library,
                     file_dir = "haplotype_sequence")
    }, silent = TRUE)
  }
}