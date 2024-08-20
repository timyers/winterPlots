# Description:  filter out and keep only one of the overlapping 
# intervals in bed file that are causing the error message, 
# "Data ranges cannot overlap.  Please check `start` and `end` 
# column ranges" from `plotSignal` #function in the 
# `plotgardener` R package.

ca_file_path <- "data/chip-atlas/filtered_chip-atlas_H3K27ac_HEK293T_hg19.bed"
ca_file_path <- "data/chip-atlas/filtered_chip-atlas_H3K27ac_ACHN_hg19.bed"


bed_file <- data.table::fread(ca_file_path, sep = "\t", header = TRUE)

# Rename columns for clarity
colnames(bed_file) <- c("chrom", "start", "end", "name")

# Convert the BED file to a GRanges object
gr <- GRanges(seqnames = bed_file$chrom,
              ranges = IRanges(start = bed_file$start, end = bed_file$end),
              name = bed_file$name)

# Find overlaps
overlaps <- findOverlaps(gr, gr)

# Identify and remove overlapping intervals
non_overlapping <- gr[-queryHits(overlaps[duplicated(queryHits(overlaps))])]

# Prepare the filtered data for output
filtered_bed <- data.frame(chrom = seqnames(non_overlapping),
                           start = start(non_overlapping),
                           end = end(non_overlapping),
                           score = mcols(non_overlapping)$name)

# Create filename
full_file_name <- "data/chip-atlas/merged_intervals_filtered_chip-atlas_H3K27ac_HEK293T_hg19.bed"
# Or
full_file_name <- "data/chip-atlas/merged_intervals_filtered_chip-atlas_H3K27ac_ACHN_hg19.bed"

# Write the filtered BED to a new file
write.table(filtered_bed, full_file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




