# Load required packages
library(Biostrings)
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(plyranges)

# Read data
protein_seq <- readAAStringSet("./data/ncbi_Lp_protein.faa")[[1]]
gff_data <- import("./data/NCBI_Lp02_annotations.gff")
hist_freq <- read.csv("./out/NCBI_Lp_histFreq.csv") |> 
  distinct()

# Process CDS annotations
cds_annotation <- gff_data[gff_data$type == "CDS"]
cds_df <- as.data.frame(cds_annotation) |> 
  drop_na(protein_id) |> 
  distinct()

# Merge data
merged_df <- merge(cds_df, hist_freq, by = "protein_id", all.x = TRUE)

# Create genomic ranges object with proper chromosome and strand information
lp_granges <- makeGRangesFromDataFrame(
  merged_df, 
  keep.extra.columns = TRUE
)
strand(lp_granges) <- "*"

genome <- unique(as.character(seqnames(lp_granges)))
options(ucscChromosomeNames = FALSE)

axisTrack <- GenomeAxisTrack(
  fontsize = 12,
  lwd = 2
)

dataTrack <- DataTrack(
  range = lp_granges,
  name = "Hist Freq/12 Genes",
  data = "Freq",
  type = "histogram",
  window = -1,
  windowSize = 12,
  fill.histogram = "black",
  col.histogram = "black",
  ylim = c(0, max(merged_df$Freq, na.rm = TRUE) * 1.1),  # Add 10% padding to y-axis
  col.axis = "black",
  fontsize = 12,
  background.panel = "white",
  cex.axis = 0.8,
  col.grid = "lightgray",
  lwd = 1.5
)


geneTrack <- GeneRegionTrack(
  range = gff_data[gff_data$type == "gene"],
  name = "Genes",
  fill = "salmon",
  col = "black",
  fontsize = 12
)

displayPars(dataTrack) <- list(
  showTitle = TRUE,
  titleWidth = 1.2,
  cex.title = 1,
  col.title = "black",
  rotation.title = 90
)

pdf("./plots/genome_hist_freq_black.pdf", width = 10, height = 7, useDingbats = FALSE)

plotTracks(
  list(axisTrack, dataTrack),
  from = 1,
  to = max(end(lp_granges)),
  cex.main = 1.5,
  col.main = "black",
  background.title = "white",
  fontcolor.title = "black",
  sizes = c(0.2, 0.3),
  innerMargin = 0.01
)

dev.off()



