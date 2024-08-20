# rs3087660
# hg19 used because Leandro's data from 2019
# was generated using that build

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("plotgardener")
BiocManager::install("plotgardenerData")

## Libraries
library(rtracklayer)
library(GenomicRanges)
library(plotgardener)
library(plotgardenerData)
library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# data("IMR90_HiC_10kb")
# data("IMR90_DNAloops_pairs")

################################
# Original location of Leandro's Capture Hi-C data files
# /Volumes/ifs/DCEG/Branches/LTG/Chanock/Leandro/MPRA_kidney/CaptureC/HEK806_2/data/MPRA.SNPs
################################

# read Capture Hi-C data
rs3087660_hek <- read.table("data/captureC/HEK/washU.chr15.rs3087660")
rs3087660_achn <- read.table("data/captureC/ACHN/washU.chr15.rs3087660")

# Create page
pageCreate(width = 6.5, height = 8.0,
           showGuides = TRUE,
           default.units = "inches")

# Set the coordinates
params <- pgParams(chrom = "chr15",
                   # chromstart = 66374958,
                   chromstart = 66395000,
                   # chromend = 67050000,
                   chromend = 67060000,
                   assembly = "hg19",
                   width = 5.25
                  )

# arch plots for rs3087660 HEK293T
my_arches_h <- plotPairsArches(
  data = rs3087660_hek,
  params = params,
  assembly = "hg19",
  archHeight = "V7", alpha = 1,
  x = 0.75, y = 1.5, width = 5.25, height = 1.5,
  just = c("left", "top"), default.units = "inches",
  fill = "#FFFFFF", linecolor = "blue", flip = FALSE
)

## Annotate genome label
annoGenomeLabel(plot = my_arches_h,
                x = 0.75, y =3.0, 
                scale = "Mb",
                # at = c(66400000, 66600000, 66800000, 67000000),
                at = c(66679251, 66783882), # add tick marks for MAP2K1
                tcl = 0.75
)

# rs3087660
annoYaxis(plot = my_arches_h, 
          at = c(0,5,10,15,20,25,30),
          fontsize = 10,
          axisLine = TRUE
)

# arch plots for rs3087660 ACHN
my_arches_a <- plotPairsArches(
  data = rs3087660_achn,
  params = params,
  assembly = "hg19",
  archHeight = "V7", alpha = 1,
  x = 0.75, y = 1.5, width = 5.25, height = 1.5,
  just = c("left", "top"), default.units = "inches",
  fill = "#FFFFFF", linecolor = "red", flip = FALSE
)

## Add label for Capture-C arch plot
plotText(label = "Capture-C Score",
         fontsize = 11,
         fontface = "bold",
         rot = 90,
         x = 0.25,
         y = 2.2
        ) 

# rs3087660
plotText(label = "rs3087660",
         fontsize = 14,
         fontface = "bold",
         x = 3.5,
         y = 0.15
)

# Add chromosome label
# chr15 - rs3087660
plotText(label = "chr15",
         fontsize = 8,
         x = 1.2,
         y = 0.75
        )

## Add ideogram
# chr15 - rs3087660
ideogramPlot <- plotIdeogram(
  chrom = "chr15", assembly = "hg19",
  orientation = "h",
  x = 1.0, y = 0.5, width = 4.75, height = 0.3, just = "left"
)

## Add highlight to ideogram
# chr15 - rs3087660
region <- pgParams(chrom = "chr15", chromstart = 66374958, chromend = 67050000 )
annoHighlight(
  plot = ideogramPlot, params = region,
  fill = "red",
  y = 0.25, height = 0.5, just = c("left", "top"), default.units = "inches"
)

## Add zoom lines
annoZoomLines(
  plot = ideogramPlot, params = region,
  y0 = 0.75, x1 = c(0.85, 6.0), y1 = 1.5,
  linecolor = "black",
  default.units = "inches"
)

# Add legend plot
# legendPlot <- plotLegend(legend = c("HEK293T", "ACHN"),
#                          fill = c("blue", "red"),
#                          fontsize = 8,
#                          border = FALSE,
#                          x = 4.75, y = 1.5,
#                          width = 1.3, height = 0.5,
#                          just = c("left", "top"),
#                          default.units = "inches"
#                         )

## Plot gene track
# rs3087660
orderGenes <- c("TIPIN", "MAP2K1", "SMAD6", "RPL4",
                "LCTL", "DIS3L", "SNORD16", "SCARNA14")

genePlot <- plotGenes(params = params,
                      x = 0.75, y = 3.5, height = 1,
                      geneOrder = orderGenes,
                      bg = "lightgray",
                      geneHighlights = data.frame(gene = "MAP2K1", color = "#E69F00"),
                      geneBackground = "black"
                     )     

# Add label for Genes track
plotText(label = "Genes",
         fontsize = 11,
         fontface = "bold",
         rot = 90,
         x = 0.25,
         y = 4
        ) 


## Plot signal track
# Read in chip-atlas file for chr15
########## hg38 -> need to convert to hg19 ##########
# file_path <- "data/chip-atlas/subset_chr15_hg19_parsed_Histone.Kidney.100.H3K27ac.AllCell.bed.txt"
# file above did not work. error msg about not allowing overlaps?

# original encode url https://www.encodeproject.org/experiments/ENCSR000DTU/
# H3K4me3, HEK293, filter by hg19 and bigwig files;
# downloaded bigwig file ENCFF110KNX.bigwig (signal p-value, 
# isogenic replicate #1, hg19) and converted to .bed file using 
# script `bitwig_to_bed.R`
file_path <- "data/encode/encode_H3K4me3_HEK293_hg19_ENCFF110KNX.bed" # this worked
df_bed <- data.table::fread(file_path, sep = "\t", header = TRUE)

# Plot signal track
plotSignal(
  data = df_bed,
  params = params,
  # chrom = "chr15",
  # chromstart = 66395000,
  # chromend = 67050000,
  # chromend = 67060000,
  # chrom = "chr21", chromstart = 28000000, chromend = 30300000,
  # assembly = "hg19",
  x = 0.75, y = 4.655, 
  width = 5.25, height = 0.5,
  just = c("left", "top"), default.units = "inches"
)

# Add label for Encode track
plotText(label = "Encode",
         fontsize = 11,
         fontface = "bold",
         rot = 90,
         x = 0.25,
         y = 4.9
        )

# Add sub-label for Encode track
plotText(label = "H3K4me3",
         fontsize = 8,
         fontface = "bold",
         rot = 90,
         x = 0.375,
         y = 4.9
) 

## Add SNP track
# SNP rs3087660 data
snp_data <- data.frame(chr = "chr1",
                       start = 66797492,
                       end = 66797492,
                       snp_id = "rs3087660"
)
# chr15 - rs3087660
snp_region <- pgParams(chrom = "chr15", 
                       chromstart = 66794000,
                       # chromstart = 66796492,
                       chromend = 66797492)
annoHighlight(
  plot = genePlot, params = snp_region,
  fill = "gray36",
  y = 3.25, height = 2.70, just = c("left", "top"), default.units = "inches"
)

# Add SNP rsID label
plotText(label = "rs3087660",
         fontsize = 7,
         x = 3.85,
         y = 6.02
)

# Add label for SNP track
plotText(label = "SNP",
         fontsize = 11,
         fontface = "bold",
         rot = 90,
         x = 0.25,
         y = 6.0
)

## Hide page guides
pageGuideHide()
pageGuideShow()

