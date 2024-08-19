# Convert bigwig (signal p-value) to bed file

library(rtracklayer)

bigwig_file <- "/Users/myersta/Downloads/ENCFF110KNX.bigWig" # hg19
granges <- import(bigwig_file, format = "bigWig")

# bed_file <- "/Users/myersta/Downloads/encode_file.bed"

# filter by desired threshold (-log10(p-value))
threshold <- 3  
filtered_granges <- granges[score(granges) >= threshold]

# Convert the GRanges object to a data frame
result <- as.data.frame(filtered_granges)

# Clean up column names and return
# delete columns not needed
result$width <- NULL
result$strand <- NULL
colnames(result) <- c("chrom", "start", "end", "score")

# write bigwig to bed to file
# set new filename
bed_file <- "/Users/myersta/Downloads/encode_file.bed"

# write dataframe to file
write.table(result, file = bed_file, 
            sep = "\t", row.names = FALSE, 
            col.names = TRUE, quote = FALSE
)
