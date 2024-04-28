## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

has_pkg = requireNamespace("dplyr", quietly = TRUE) &&
  requireNamespace("tidyr", quietly = TRUE) &&
  requireNamespace("gplots", quietly = TRUE)
knitr::opts_chunk$set(eval = has_pkg)

## -----------------------------------------------------------------------------
library(scPloidy)
library(dplyr)
library(tidyr)
library(gplots)

## ----eval = FALSE-------------------------------------------------------------
#  SU008_Tumor_Pre_windowcovariates =
#    read.table(
#      "multi_tissue_peaks.hg37.20MB.bed",
#      header = FALSE)
#  colnames(SU008_Tumor_Pre_windowcovariates) =
#    c("chr", "start", "end", "window", "peaks")

## ----eval = FALSE-------------------------------------------------------------
#  simpleRepeat = readr::read_tsv(
#    "~/human/publichuman/hg37_ucsc/simpleRepeat.chrom_chromStart_chromEnd.txt.gz",
#    col_names = c("chrom", "chromStart", "chromEnd"))
#  rmsk = readr::read_tsv(
#    "~/human/publichuman/hg37_ucsc/rmsk.Simple_repeat.genoName_genoStart_genoEnd.txt.gz",
#    col_names = c("chrom", "chromStart", "chromEnd"))
#  simpleRepeat = rbind(simpleRepeat, rmsk)
#  rm(rmsk)
#  # convert from 0-based position to 1-based
#  simpleRepeat[, 2] = simpleRepeat[, 2] + 1
#  simpleRepeat = GenomicRanges::makeGRangesFromDataFrame(
#    as.data.frame(simpleRepeat),
#    seqnames.field = "chrom",
#    start.field    = "chromStart",
#    end.field      = "chromEnd")
#  # remove duplicates
#  simpleRepeat = GenomicRanges::union(simpleRepeat, GenomicRanges::GRanges())

## ----eval = FALSE-------------------------------------------------------------
#  window = read.table("window.hg37.20MB.bed", header = FALSE)
#  colnames(window) = c("chr", "start", "end", "window")
#  at = GenomicRanges::makeGRangesFromDataFrame(window[, 1:3])
#  barcodesuffix = paste0(".", window$window)

## ----eval = FALSE-------------------------------------------------------------
#  sc = read.csv(
#    "GSE129785_scATAC-TME-All.cell_barcodes.txt",
#    header = TRUE,
#    sep = "\t")

## ----eval = FALSE-------------------------------------------------------------
#  sample = "GSM3722064"
#  tissue = "SU008_Tumor_Pre"
#  
#  bc = sc$Barcodes[sc$Group == tissue]
#  SU008_Tumor_Pre_fragmentoverlap =
#    fragmentoverlapcount(
#      paste0("SRX5679934/", sample, "_", tissue, "_fragments.tsv.gz"),
#      at,
#      excluderegions = simpleRepeat,
#      targetbarcodes = bc,
#      Tn5offset = c(1, 0),
#      barcodesuffix = barcodesuffix
#    )

## -----------------------------------------------------------------------------
data(GSE129785_SU008_Tumor_Pre)

## -----------------------------------------------------------------------------
levels = c(2, 4)
result = cnv(SU008_Tumor_Pre_fragmentoverlap,
             SU008_Tumor_Pre_windowcovariates,
             levels = levels,
             deltaBICthreshold = -600)

## -----------------------------------------------------------------------------
windowcovariates = SU008_Tumor_Pre_windowcovariates
windowcovariates$w =
  as.numeric(sub("window_", "", windowcovariates$window))

fragmentoverlap = SU008_Tumor_Pre_fragmentoverlap
fragmentoverlap$cell =
  sub(".window.*", "", fragmentoverlap$barcode)
fragmentoverlap$window =
  sub(".*window", "window", fragmentoverlap$barcode)
fragmentoverlap$w =
  as.numeric(sub("window_", "", fragmentoverlap$window))

x = match(fragmentoverlap$barcode,
        result$cellwindowCN$barcode)
fragmentoverlap$CN = result$cellwindowCN$CN[x]
fragmentoverlap$ploidy.moment.cell = result$cellwindowCN$ploidy.moment.cell[x]
fragmentoverlap = fragmentoverlap[!is.na(fragmentoverlap$CN), ]

# For better hierarchical clustering
fragmentoverlap$pwindownormalizedcleanedceiled =
  pmin(fragmentoverlap$CN, min(levels) * 2)

## -----------------------------------------------------------------------------
dataplot =
  fragmentoverlap %>%
  dplyr::select("w", "cell", "pwindownormalizedcleanedceiled") %>%
  tidyr::pivot_wider(names_from = "w", values_from = "pwindownormalizedcleanedceiled")
dataplot = as.data.frame(dataplot)
rownames(dataplot) = dataplot$cell
dataplot = dataplot[, colnames(dataplot) != "cell"]
dataplot = as.matrix(dataplot)
n = max(as.numeric(colnames(dataplot)))
dataplot = dataplot[, match(as.character(1:n), colnames(dataplot))]
colnames(dataplot) = as.character(1:n)

## ----eval = FALSE-------------------------------------------------------------
#  breaks = c(0, min(levels) - 1, min(levels) + 1, min(levels) * 2)
#  
#  x = windowcovariates
#  x$chr[duplicated(windowcovariates$chr)] = NA
#  x = x$chr[match(colnames(dataplot), x$w)]
#  
#  RowSideColors =
#    unlist(
#      lapply(
#        fragmentoverlap$ploidy.moment.cell[
#          match(rownames(dataplot), fragmentoverlap$cell)],
#        function (x) { which(sort(levels) == x)}))
#  RowSideColors = topo.colors(length(levels))[RowSideColors]
#  
#  gplots::heatmap.2(
#    dataplot,
#    Colv = FALSE,
#    dendrogram = "none",
#    breaks = breaks,
#    col = c("blue", "gray80", "red"),
#    trace = "none", labRow = FALSE, na.color = "white",
#    labCol = x,
#    RowSideColors= RowSideColors)

