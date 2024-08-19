if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
library(GenomicRanges)

# Assuming your data frame is called df_bed
gr <- GRanges(seqnames = df_bed$chrom,
              ranges = IRanges(start = df_bed$start, end = df_bed$end),
              score = df_bed$score)

overlaps <- findOverlaps(gr, gr)
overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]
unique_overlaps <- unique(queryHits(overlaps))
num_overlapping_regions <- length(unique_overlaps)

cat("Number of overlapping regions:", num_overlapping_regions, "\n")
