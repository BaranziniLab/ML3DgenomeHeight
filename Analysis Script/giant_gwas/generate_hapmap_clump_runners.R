library(dplyr)

hapmap_template = "~/software/plink2 --bfile /salemlab/REF/TopMed_Ref/topmed_ref_ANCESTRY --thin-indiv-count SAMPLE_SIZE --extract hapmap_snplist.txt --maf 0 --clump ANCESTRY_plink.txt --clump-p1 5e-08 --clump-r2 0.1 --clump-kb 500 --threads NCORES --out hapmap_result/ANCESTRY_SAMPLE_SIZE_hapmap_clumping_result"
tgp_template = "~/software/plink2 --bfile /salemlab/REF/TopMed_Ref/topmed_ref_ANCESTRY --thin-indiv-count SAMPLE_SIZE --extract tgp_snplist.txt --maf 0 --clump ANCESTRY_plink.txt --clump-p1 5e-08 --clump-r2 0.1 --clump-kb 500 --threads NCORES --out tgp_result/ANCESTRY_SAMPLE_SIZE_tgp_clumping_result"
topmed_template = "~/software/plink2 --bfile /salemlab/REF/TopMed_Ref/topmed_ref_ANCESTRY --thin-indiv-count SAMPLE_SIZE --maf 0 --clump ANCESTRY_plink.txt --clump-p1 5e-08 --clump-r2 0.1 --clump-kb 500 --threads NCORES --out topmed_result/ANCESTRY_SAMPLE_SIZE_topmed_clumping_result"

afr_n = c(seq(100, 1000, 100), seq(1000, 10000, 1000), 11000, 11720)
amr_n = c(seq(100, 1000, 100), seq(1000, 10000, 1000), 10381)
eas_n = c(seq(100, 800, 100), 893)
eur_n = c(seq(100, 1000, 100), seq(1000, 10000, 1000), 15000, 20000, 25000, 27322)
sas_n = c(100, 103)

command_snippet = ""
for(n in afr_n){
  command_snippet = paste(command_snippet,
  gsub(pattern = "ANCESTRY", replacement = "afr", hapmap_template) %>%
    gsub(pattern = "SAMPLE_SIZE", replacement = n) %>%
    gsub(pattern = "NCORES", replacement = 40),
  sep = "\n\n"
  )
}

command_snippet = paste("#! /bin/bash
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=45
#SBATCH --ntasks=1
#SBATCH -t 21-00:00
#SBATCH --job-name  clump_afr_hapmap
#SBATCH --output=clump_afr_hapmap.out
#SBATCH --partition=salem-compute

date", command_snippet, "date", sep = "\n\n")

afr_commands = command_snippet
rm(command_snippet)





command_snippet = ""
for(n in amr_n){
  command_snippet = paste(command_snippet,
                          gsub(pattern = "ANCESTRY", replacement = "amr", hapmap_template) %>%
                            gsub(pattern = "SAMPLE_SIZE", replacement = n) %>%
                            gsub(pattern = "NCORES", replacement = 40),
                          sep = "\n\n"
  )
}

command_snippet = paste("#! /bin/bash
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=45
#SBATCH --ntasks=1
#SBATCH -t 21-00:00
#SBATCH --job-name  clump_amr_hapmap
#SBATCH --output=clump_amr_hapmap.out
#SBATCH --partition=salem-compute

date", command_snippet, "date", sep = "\n\n")

amr_commands = command_snippet
rm(command_snippet)




command_snippet = ""
for(n in eur_n){
  command_snippet = paste(command_snippet,
                          gsub(pattern = "ANCESTRY", replacement = "eur", hapmap_template) %>%
                            gsub(pattern = "SAMPLE_SIZE", replacement = n) %>%
                            gsub(pattern = "NCORES", replacement = 70),
                          sep = "\n\n"
  )
}

command_snippet = paste("#! /bin/bash
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=75
#SBATCH --ntasks=1
#SBATCH -t 21-00:00
#SBATCH --job-name  clump_eur_hapmap
#SBATCH --output=clump_eur_hapmap.out
#SBATCH --partition=salem-compute

date", command_snippet, "date", sep = "\n\n")

eur_commands = command_snippet
rm(command_snippet)





command_snippet = ""
for(n in eas_n){
  command_snippet = paste(command_snippet,
                          gsub(pattern = "ANCESTRY", replacement = "eas", hapmap_template) %>%
                            gsub(pattern = "SAMPLE_SIZE", replacement = n) %>%
                            gsub(pattern = "NCORES", replacement = 20),
                          sep = "\n\n"
  )
}

command_snippet = paste("#! /bin/bash
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=25
#SBATCH --ntasks=1
#SBATCH -t 21-00:00
#SBATCH --job-name  clump_eas_hapmap
#SBATCH --output=clump_eas_hapmap.out
#SBATCH --partition=salem-compute

date", command_snippet, "date", sep = "\n\n")

eas_commands = command_snippet
rm(command_snippet)





command_snippet = ""
for(n in sas_n){
  command_snippet = paste(command_snippet,
                          gsub(pattern = "ANCESTRY", replacement = "sas", hapmap_template) %>%
                            gsub(pattern = "SAMPLE_SIZE", replacement = n) %>%
                            gsub(pattern = "NCORES", replacement = 20),
                          sep = "\n\n"
  )
}

command_snippet = paste("#! /bin/bash
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=25
#SBATCH --ntasks=1
#SBATCH -t 21-00:00
#SBATCH --job-name  clump_sas_hapmap
#SBATCH --output=clump_sas_hapmap.out
#SBATCH --partition=salem-compute

date", command_snippet, "date", sep = "\n\n")

sas_commands = command_snippet
rm(command_snippet)


writeLines(afr_commands, con = "afr_hapmap_clump_runner.sh")
writeLines(amr_commands, con = "amr_hapmap_clump_runner.sh")
writeLines(eur_commands, con = "eur_hapmap_clump_runner.sh")
writeLines(sas_commands, con = "sas_hapmap_clump_runner.sh")
writeLines(eas_commands, con = "eas_hapmap_clump_runner.sh")


