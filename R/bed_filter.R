# Description:  filters parsed ChIP-Atlas bed file (hg38) from
# my `UpSetToolkit` project by chromosome and cell type
# for use as input for signal track in project 
# `winterPlots`.  Output file may need to be converted
# to hg19 using script `bed_liftoverR.R`

# Load the necessary library
library(readr)
library(dplyr)


# define chip-atlas assay to be filtered
chip_atlas_assay <- "H3K27ac"

# set genome build
genome_build <- "hg38" # chip-atlas data was downloaded with hg38 coordinates

# define parsed ChIP-Atlas .bed file to filter
# this file includes data from multiple cell lines
bed_filename <- "/Users/myersta/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/Rprojects_OneDrive/winter/UpSetTookit/data/output/H3K27ac/parsed_Histone.Kidney.100.H3K27ac.AllCell.bed.txt"

# Read the BED file
bed_data <- read_tsv(bed_filename, col_names = TRUE)

# define filters
chrom_filter <- "chr15"
sample_filter <- "HEK293T"
# Orrr
sample_filter <- "ACHN"

# Filter rows where `chrom_filter` is in the "chrom" column and
# `sample_filter` appears in any column
filtered_data <- bed_data %>%
  filter(chrom == chrom_filter & apply(., 1, function(row) any(grepl(sample_filter, row))))


# Select and rename columns
final_data <- filtered_data %>%
  select(chrom, chromStart, chromEnd, score) %>%
  rename(start = chromStart, end = chromEnd)

# Save the final data to a new BED file
file_path <- "data/chip-atlas/"
full_file_name <- paste0(file_path,
                         "filtered_chip-atlas_",
                         chip_atlas_assay, "_",
                         sample_filter, "_",
                         genome_build,
                         ".bed"
                         )
write.table(final_data, full_file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
