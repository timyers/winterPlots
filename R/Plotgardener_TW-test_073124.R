if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("plotgardener")
BiocManager::install("plotgardenerData")

## Load plotgardener
library(plotgardener)

## Load example Hi-C data
library(plotgardenerData)


data("IMR90_HiC_10kb")
data("IMR90_DNAloops_pairs")

pageCreate(
     width = 3, height = 5, default.units = "inches",
     showGuides = FALSE, xgrid = 0, ygrid = 0
)
plotHicSquare(
  data = IMR90_HiC_10kb,
  chrom = "chr21", chromstart = 28000000, chromend = 30300000,
  assembly = "hg19",
  x = 0.5, y = 0.5, width = 2, height = 2,
  just = c("left", "top"), default.units = "inches"
)
plotPairsArches(
   data = IMR90_DNAloops_pairs,
   chrom = "chr21", chromstart = 28000000, chromend = 30300000,
   assembly = "hg19",
   x = 0.5, y = 2.5, width = 2, height = 0.25,
   just = c("left", "top"), default.units = "inches",
   fill = "black", linecolor = "black", flip = TRUE
)


rs8133 <- read.table("washU.chr1.rs8133")
pageCreate(
  width = 3, height = 5, default.units = "inches",
  showGuides = FALSE, xgrid = 0, ygrid = 0
)
plotPairsArches(
  data = rs8133,
  chrom = "chr1", chromstart = 165600000, chromend = 165800000,
  assembly = "hg19",
  x = 0.5, y = 2.5, width = 2, height = 0.25,
  just = c("left", "top"), default.units = "inches",
  fill = "black", linecolor = "black", flip = TRUE
)