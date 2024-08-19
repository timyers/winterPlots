# Description: function to convert bed file with hg38 
# coordinates to hg19 and save to text file.

library(rtracklayer)
library(GenomicRanges)
library(R.utils)

# load function
convert_hg38_to_hg19 <- function(bed_data) {
  
  # Convert the data frame to a GRanges object
  gr <- GRanges(seqnames = bed_data$chrom,
                ranges = IRanges(start = bed_data$chromStart, end = bed_data$chromEnd),
                score = bed_data$score)
  
  # Chain file for hg38 to hg19 conversion
  # downloaded from UCSC Genome Browser
  # URL: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
  # File name for decompressed chain file
  chain_file <- "data/liftover_chain_file/hg38ToHg19.over.chain"
  
  # Import the decompressed chain file
   chain <- import.chain(chain_file)
  
  # Perform the liftOver
  gr_hg19 <- liftOver(gr, chain)
  gr_hg19 <- unlist(gr_hg19) # Unlist to get a GRanges object directly
  
  # Convert the GRanges object back to a data frame
  result <- as.data.frame(gr_hg19)
  
  # Clean up column names and return
  # result <- result[, c("seqnames", "start", "end", "score")]
  # delete columns not needed
  result$width <- NULL
  result$strand <- NULL
  colnames(result) <- c("chrom", "start", "end", "score")
  
  return(result)
}

# Example usage
bed_data <- data.frame(
  chrom = c("chr15", "chr15", "chr15"),
  chromStart = c(20100537, 20100540, 20100654),
  chromEnd = c(20101338, 20101612, 20101387),
  score = c(259, 273, 462)  # Numeric values
)

# read in file
file_path <- "data/chip-atlas/subset_chr15_hg38_parsed_Histone.Kidney.100.H3K27ac.AllCell.bed.txt"
df_bed <- data.table::fread(file_path, sep = "\t", header = TRUE)

# call function to convert coordinates
result_hg19 <- convert_hg38_to_hg19(df_bed)

# write file converted to hg19
# set new filename
file_path_hg19 <- gsub("hg38", "hg19", file_path)

# write dataframe to file
write.table(result_hg19, file = file_path_hg19, 
            sep = "\t", row.names = FALSE, 
            col.names = TRUE, quote = FALSE
           )


print(result_hg19)
