library(dplyr)
library(data.table)

afr = fread("afr_1000_tgp_clumping_result.clumps")
eur = fread("eur_1000_tgp_clumping_result.clumps")
amr = fread("amr_1000_tgp_clumping_result.clumps")
eas = fread("eas_893_tgp_clumping_result.clumps")
sas = fread("sas_103_tgp_clumping_result.clumps")


make_map = function(d, file_name){
  d = d[,c(1,2)]
  names(d) = c("chr", "pos")
  d[["start"]] = d[["pos"]] - 1000000
  d[["end"]] = d[["pos"]] + 1000000
  d = d %>%
    arrange(chr, start) %>%
    mutate(start = ifelse(start < 0, 0, start)) %>%
    mutate(end = ifelse(end < 0, 0, end))
  d = d[,-2]
  cat(paste0("Making ", file_name, "...\n"))
  fwrite(x = d, file = file_name, sep = " ", col.names = FALSE)
}

make_map(afr, file_name = "afr_height_loci.map")
make_map(amr, file_name = "amr_height_loci.map")
make_map(eas, file_name = "eas_height_loci.map")
make_map(eur, file_name = "eur_height_loci.map")
make_map(sas, file_name = "sas_height_loci.map")



