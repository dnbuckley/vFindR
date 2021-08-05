require(Rsamtools)
require(rtracklayer)
require(parallel)
require(dplyr)
source("src/split_read_table.R")
source("src/chimeric_reads_table.R")
source("~/misc_R_scripts/combineInRange.R")

output.dir <- "vFindR_VIcaller-example-output/"

vFindRAmalgamate <- function(output.dir, combine.limit = 300) {
  sreads <- .splitReadTable(output.dir)
  sregions <- combineInRange(GRanges(sreads$read.location, strand = "*"), combine.limit)
  vchrs <- as.character(unique(seqnames(GRanges(sreads$virus.read.location))))
  hchrs <- as.character(unique(seqnames(sregions)))
  creads <- .chimericReadTable(output.dir)
  creads$strand <- "*"
  # chiral reads with split evidence
  cReads.SE <- subsetByOverlaps(GRanges(creads), sregions)$read.name
  creads.filter1 <- creads[creads$read.name %in% cReads.SE, ]
  creads.filter1ByRead <- split(creads.filter1, creads.filter1$read.name)
  creads.filter1ByRead <- creads.filter1ByRead[sapply(creads.filter1ByRead,
                                                      function(x, vc, hc) {
                                                      any(x$seqnames %in% vc) & any(x$seqnames %in% hc)
                                                        }, vchrs, hchrs)]
  creads.filter2 <- do.call(rbind.data.frame, creads.filter1ByRead)
  res = as.data.frame(sregions)
  res$num_split_reads
}