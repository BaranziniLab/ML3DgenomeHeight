library(dplyr)
library(Biostrings)

load("haplotyping_results.rda")
fasta_file = "ref_sequence/hg38_reference.fa"
index_file = "ref_sequence/hg38_reference.fa.fai"
fasta_sequences = readDNAStringSet(fasta_file, use.names = TRUE)

hit_to_sequence = function(hit, mhc = 10){
  
  sequence_library = NULL
  try({
    
    if(!is.list(hit)){
      return(NULL)
    }
    
    chromosome = paste0("chr", hit$fix$CHROM[1])
    hap_range = range(as.numeric(hit$fix$POS))
    hap_length = diff(hap_range)
    
    if(hap_length >= 2^20){
      return(NULL) 
    }
    
    if(dim(hit[["haplotype_table"]])[1] == 0){
      return(NULL)
    }
    
    remaining_hap_index = hit[["haplotype_table"]][[tolower(hit[["ancestry"]])]] >= mhc
    hit[["haplotype_table"]] = hit[["haplotype_table"]][remaining_hap_index,]
    
    if(dim(hit[["haplotype_table"]])[1] == 0){
      return(NULL)
    }
    
    n_base_edit = dim(hit$fix)[1]
    haplotype_table = hit$haplotype_table
    hap_from = hap_range[1]
    hap_to = hap_range[2]
    
    hap_sequence = fasta_sequences[which(fasta_sequences@ranges@NAMES == chromosome)] %>%
      subseq(start = hap_from, end = hap_to) %>%
      as.character() %>%
      strsplit(split = "") %>%
      unlist()
    names(hap_sequence) = NULL
    
    ref_sequence = paste(hap_sequence, collapse = "")
    
    for(i in 1:n_base_edit){
      snp_meta = hit$fix[i,]
      snp_pos = as.numeric(snp_meta$POS)
      snp_index = snp_pos - hap_from + 1
      hap_sequence[snp_index] = paste(c("x",i,"x"), collapse = "")
    }
    
    sequence_library = data.frame(
      ancestry = hit$ancestry,
      leading_marker = hit$leading_marker,
      haplotype = NA,
      sequence = ref_sequence
    )
    
    for(h in haplotype_table$haplotype){
      
      alt_sequence = hap_sequence
      if(nchar(h) == n_base_edit){
        haplotype = unlist(strsplit(h, split = ""))
        
        for(i in 1:n_base_edit){
          anchor = paste(c("x",i,"x"), collapse = "")
          alt_sequence[alt_sequence == anchor] = haplotype[i]
        }
        alt_sequence = paste(alt_sequence, collapse = "")
      }
      
      sequence_library = rbind.data.frame(
        sequence_library,
        data.frame(
          ancestry = hit$ancestry,
          leading_marker = hit$leading_marker,
          haplotype = h,
          sequence = alt_sequence
        )
      )
      
    }
    
    sequence_library[["sequence_type"]] = c("REF", paste0("ALT", 1:dim(haplotype_table)[1]))
    
    padding_length = (2^20 - length(hap_sequence)) / 2
    
    seq_from = hap_from - ceiling(padding_length)
    seq_to = hap_to + floor(padding_length)
    
    if(seq_from < 1){
      return(NULL)
    }
    
    upstream_padding = fasta_sequences[which(fasta_sequences@ranges@NAMES == chromosome)] %>%
      subseq(start = seq_from, end = hap_from - 1) %>%
      as.character() %>%
      strsplit(split = "") %>%
      unlist() %>%
      paste(collapse = "")
    
    downstream_padding = fasta_sequences[which(fasta_sequences@ranges@NAMES == chromosome)] %>%
      subseq(start = hap_to + 1, end = seq_to) %>%
      as.character() %>%
      strsplit(split = "") %>%
      unlist() %>%
      paste(collapse = "")
    
    sequence_library[["upstream_padding"]] = upstream_padding
    sequence_library[["downstream_padding"]] = downstream_padding
    
    sequence_library[["full_sequence"]] = paste(
      sequence_library[["upstream_padding"]],
      sequence_library[["sequence"]],
      sequence_library[["downstream_padding"]],
      sep = ""
    )
    
    sequence_library[["hap_from"]] = hap_from
    sequence_library[["hap_to"]] = hap_to
    sequence_library[["hap_length"]] = hap_length
    
    sequence_library[["file_name"]] = paste0(paste(
      sequence_library[["sequence_type"]],
      sequence_library[["ancestry"]],
      gsub(pattern = "[:]", replacement = "_", sequence_library[["leading_marker"]]),
      sep = "_"
    ), ".txt")
    
  }, silent = TRUE)
  
  return(sequence_library)
  
}

write_sequence = function(sequence_library, output_dir = "haplotype_sequence"){
  leader_marker = sequence_library$leading_marker[1]
  ancestry = sequence_library$ancestry[1]
  cat(paste0("Writting haplotypes tagged by ", leader_marker, " into file...\n"))
  dir_name = paste0(ancestry, "_", gsub(pattern = "[:]", replacement = "_", leader_marker))
  if(!dir.exists(file.path(output_dir, dir_name))){
    dir.create(file.path(output_dir, dir_name))
    for(i in 1:dim(sequence_library)[1]){
      full_sequence = sequence_library$full_sequence[i]
      file_name = sequence_library$file_name[i]
      writeLines(text = full_sequence, con = file.path(output_dir, dir_name, file_name),
                 sep = "")
    }
  }
}

insilico_mutagen = function(leading_marker, output_dir){
  
  gene_locus = haplotyping_results[[which(names(haplotyping_results) == leading_marker)]]
  locus_fix = gene_locus$fix[,c(4,5)]
  
  allele_list = list()
  for(i in 1:dim(locus_fix)[1]){
    allele_list[[i]] = unlist(locus_fix[i,])
  }
  
  haplotype = expand.grid(allele_list) %>%
    apply(MARGIN = 1, function(x){
      paste(x, collapse = "")
    })
  
  gene_locus$haplotype_table = tibble(
    haplotype = haplotype,
    afr = 100,
    amr = 100,
    eur = 100,
    eas = 100,
    sas = 100
  )
  
  gene_sequence = hit_to_sequence(gene_locus)
  
  if(!is.null(gene_sequence)){
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    
    write_sequence(sequence_library = gene_sequence,
                   output_dir = output_dir)
  }
}

insilico_mutagen(leading_marker = "10:96789820:T:C", output_dir = "lcor_insilico")



