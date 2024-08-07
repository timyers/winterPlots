if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("plotgardener")
BiocManager::install("plotgardenerData")

## Libraries
library(plotgardener)
library(plotgardenerData)
library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# data("IMR90_HiC_10kb")
# data("IMR90_DNAloops_pairs")

################################
# Leandro's Capture Hi-C data files original location
# HEK293T:
# /Volumes/ifs/DCEG/Branches/LTG/Chanock/Leandro/MPRA_kidney/CaptureC/HEK806_2/data/MPRA.SNPs
# ACHN:
# /Volumes/ifs/DCEG/Branches/LTG/Chanock/Leandro/MPRA_kidney/CaptureC/ACHN806_2/data/MPRA.SNPs
# rs8133 <- read.table("data/example/washU.chr1.rs8133")
rs3087660 <- read.table("data/captureC/HEK/washU.chr15.rs3087660")

pageCreate(width = 6.5, height = 7.5, default.units = "inches")
# pageCreate(
#  width = 3, height = 5, default.units = "inches",
#  showGuides = FALSE, xgrid = 0, ygrid = 0
#)

## Set the coordinates
# rs8133
# params <- pgParams(
#  chrom = "chr1",
#  chromstart = 165500000, chromend = 165900000,
#  assembly = "hg19",
#  width =10
# )

# rs3087660
c_start <- min(rs3087660$V2)
c_end <- max(rs3087660$V6)

params <- pgParams(chrom = "chr15",
                   chromstart = 66374958, 
                   chromend = 67050000,
                   assembly = "hg19",
                   width = 5
                  )

my_arches <- plotPairsArches(
  data = rs3087660,
  params = params,
  assembly = "hg19",
  archHeight = "V7", alpha = 1,
  x = 0.75, y = 1.5, width = 5, height = 1.5,
  just = c("left", "top"), default.units = "inches",
  fill = "#FFFFFF", linecolor = "black", flip = FALSE
)

## Annotate genome label
annoGenomeLabel(plot = my_arches,
                x = 0.75, y =3.0, 
                scale = "Mb",
                # at = c(66400000, 66600000, 66800000, 67000000),
                at = c(66679251, 66783882), # tick marks for MAP2K1
                tcl = 0.75
                )

## Add a vertical guide at x = 0.75 inches
# pageGuideVertical(x = 0.75,
#                   linecolor = "black",
#                   default.units = "inches")

## Add standard y-axis to Hi-C plot
# rs8133
# annoYaxis(plot = my_arches, 
#          at = c(0,2,4,6,8,10,12),
#          fontsize = 10
#         )

# rs3087660
 annoYaxis(plot = my_arches, 
           at = c(0,5,10,15,20,25,30),
           fontsize = 10,
           axisLine = TRUE
          )

## Add label for Capture-C arch plot
plotText(label = "Capture-C Score",
         fontsize = 11,
         fontface = "bold",
         rot = 90,
         x = 0.25,
         y = 2.2
        ) 

## Add title
# rs8133
# plotText(label = "rs8133",
#          fontsize = 14,
#          fontface = "bold",
#          x = 3.5,
#          y = 0.15
#         )

# rs3087660
plotText(label = "rs3087660",
         fontsize = 14,
         fontface = "bold",
         x = 3.5,
         y = 0.15
)

## Add chromosome number
# chr1- rs8133
# plotText(label = "chr1",
#          fontsize = 8,
#          x = 1.2,
#          y = 0.75
#         )

# chr15 - rs3087660
plotText(label = "chr15",
         fontsize = 8,
         x = 1.2,
         y = 0.75
        )

## Add ideogram
# chr1 - rs8133
# ideogramPlot <- plotIdeogram(chrom = "chr1", assembly = "hg19",
#                              orientation = "h",
#                              x = 1.0, y = 0.5, 
#                              width = 4.75, height = 0.3, 
#                              just = "left"
#                             )

# chr15 - rs3087660
ideogramPlot <- plotIdeogram(
  chrom = "chr15", assembly = "hg19",
  orientation = "h",
  x = 1.0, y = 0.5, width = 4.75, height = 0.3, just = "left"
)

## Add highlight to ideogram
# chr1 - rs8133
# region <- pgParams(chrom = "chr1", chromstart = 165400000, chromend = 166100000 )
# annoHighlight(
#   plot = ideogramPlot, params = region,
#   fill = "red",
#   y = 0.25, height = 0.5, just = c("left", "top"), default.units = "inches"
# )

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
  y0 = 0.75, x1 = c(0.85, 5.75), y1 = 1.5,
  linecolor = "black",
  default.units = "inches"
)


## Plot gene track
# rs3087660
orderGenes <- c("TIPIN", "MAP2K1", "SMAD6", "RPL4",
                "LCTL", "DIS3L", "SNORD16", "SCARNA14")

genePlot <- plotGenes(params = params,
                      x = 0.75, y = 3.5, height = 1,
                      geneOrder = orderGenes,
                      geneHighlights = data.frame(gene = "MAP2K1", color = "#E69F00"),
                      geneBackground = "blue"
                     )     

# Add label for Genes track
plotText(label = "Genes",
         fontsize = 11,
         fontface = "bold",
         rot = 90,
         x = 0.25,
         y = 4
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
  y = 3.25, height = 1.35, just = c("left", "top"), default.units = "inches"
)

# Add SNP rsID label
plotText(label = "rs3087660",
         fontsize = 6,
         x = 3.85,
         y = 4.655
)

# Add label for SNP track
plotText(label = "SNP",
         fontsize = 11,
         fontface = "bold",
         rot = 90,
         x = 0.25,
         y = 4.625
        )

## Hide page guides
pageGuideHide()
