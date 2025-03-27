library(Biostrings)
library(ggplot2)
library(dplyr)

fasta_file <- "./data/Lp02_protein_ncbi.faa"
protein_seqs <- readAAStringSet(fasta_file)

# function to perform the sliding window count of H
count_histidines <- function(prot_seq, window_size = 6) {
  # Convert the protein sequence to a character vector of single AAs
  seq_chars <- unlist(strsplit(as.character(prot_seq), split = ""))
  n <- length(seq_chars)
  
  # If the sequence is shorter than the window, return NA
  if(n < window_size) {
    return(NULL)
  }
  
  # Slide the window along the sequence and count "H" in each window.
  counts <- sapply(1:(n - window_size + 1), function(i) {
    window <- seq_chars[i:(i + window_size - 1)]
    sum(window == "H")
  })
  
  # Return a data frame with the window positions and the count of H.
  data.frame(window_start = 1:(n - window_size + 1),
             window_end   = window_size:(n),
             count_H      = counts)
}

# Apply the sliding window analysis to each protein CDS separately.
results_list <- lapply(names(protein_seqs), function(seq_name) {
  df <- count_histidines(protein_seqs[[seq_name]], window_size = 6)
  if (!is.null(df)) {
    df$CDS <- seq_name  # annotate which CDS this window belongs to
  }
  df
})

hist_freqs <- sliding_window_ORF$count_H






