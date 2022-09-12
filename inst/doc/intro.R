## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scPloidy)

## -----------------------------------------------------------------------------
?fragmentoverlapcount
?ploidy

## -----------------------------------------------------------------------------
targetregions =
  GenomicRanges::GRanges(
    c("chr19", "chr20"),
    IRanges::IRanges(c(1, 1), width = 500000000))

## -----------------------------------------------------------------------------
simpleRepeat = readr::read_tsv(
  system.file("extdata", "simpleRepeat.chr19_20.txt.gz", package = "scPloidy"),
  col_names = c("chrom", "chromStart", "chromEnd"))
simpleRepeat[, 2] = simpleRepeat[, 2] + 1 # convert from 0-based position to 1-based
simpleRepeat = GenomicRanges::makeGRangesFromDataFrame(
  as.data.frame(simpleRepeat),
  seqnames.field = "chrom",
  start.field    = "chromStart",
  end.field      = "chromEnd")

## -----------------------------------------------------------------------------
fragmentoverlap =
  fragmentoverlapcount(
    system.file("extdata", "SHR_m154211.10cells.chr19_20.fragments.txt.gz", package = "scPloidy"),
    targetregions,
    excluderegions = simpleRepeat,
    Tn5offset = c(4, -5))
fragmentoverlap
rm(fragmentoverlap)

## -----------------------------------------------------------------------------
data(SHR_m154211)
?SHR_m154211
fragmentoverlap = SHR_m154211$fragmentoverlap
fragmentoverlap

## -----------------------------------------------------------------------------
p = ploidy(fragmentoverlap,
           c(2, 4, 8))
head(p)

## -----------------------------------------------------------------------------
cells = SHR_m154211$cells
table(cells$celltype,
      p$ploidy.moment[match(cells$barcode, p$barcode)])

## ---- eval = FALSE------------------------------------------------------------
#  system.file("perl", "samtofragmentbed.pl", package = "scPloidy")

