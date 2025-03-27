library(Biostrings)

sequence = readAAStringSet("./Hist Frequency/lpg1596.fasta")
sequence <- as.character(sequence[[1]]) # convert to character vector

h_window <- function(seq, window_size = 6) {
  n <- nchar(seq)
  counts <- sapply(1:(n - window_size + 1), function(i) {
    substr(seq, i, i + window_size - 1) |> strsplit("") |> unlist() |> table()
  })
  
  histidine_counts <- sapply(counts, function(tbl) ifelse("H" %in% names(tbl), tbl["H"], 0))
  return(histidine_counts)
}

# Compute hist window
sliding_window_1596 <- h_window(sequence)

df <- data.frame(Position = 1:length(sliding_window_1596), Histidine_Count = sliding_window_1596)
write.csv(df, file = 'sliding_window_1596.csv', quote = F, row.names = F)
