require(Rsamtools)
require(rtracklayer)
require(parallel)
require(dplyr)
source("src/split_read_table.R")
source("src/chimeric_reads_table.R")
source("src/vFindR_utils.R")
source("~/misc_R_scripts/combineInRange.R")


vFindRAmalgamate <- function(output.dir, combine.limit = 300) {
  message("checking split reads")
  sreads <- .splitReadTable(output.dir)
  sregions <- .bamDF2GR(sreads)
  sregions.expanded <- GRanges(seqnames = seqnames(sregions), 
                               ranges = IRanges(start(sregions)-combine.limit,
                                                end = end(sregions)+combine.limit))
  sregions <- combineInRange(sregions, combine.limit)
  vchrs <- as.character(unique(sreads$virus.rname))
  hchrs <- as.character(unique(seqnames(sregions)))
  
  message("checking chimeric reads")
  creads <- .chimericReadTable(output.dir)
  creads.gr <- .bamDF2GR(creads)
  creads.gr$qname <- creads$qname
  
  # chiral reads with split evidence
  cReads.SE <- subsetByOverlaps(creads.gr, sregions.expanded)$qname
  creads.filter1 <- creads[creads$qname %in% cReads.SE, ]
  creads.filter1ByRead <- split(creads.filter1, creads.filter1$qname)
  
  creads.filter1ByRead <- creads.filter1ByRead[unlist(mclapply(creads.filter1ByRead,function(x, vc, hc) {
    any(x$rname %in% vc) & any(x$rname %in% hc)
    }, vchrs, hchrs, mc.cores = detectCores()-2))]
  
  creads.filter2ByRead <- creads.filter1ByRead[unlist(mclapply(creads.filter1ByRead, function(x, g) {
    y <- .bamDF2GR(x)
    any(from(findOverlaps(GRanges(y) , g)))
  }, sregions.expanded, mc.cores = detectCores()-2))]
  
  sl <- split(sreads, from(findOverlaps(.bamDF2GR(sreads), sregions.expanded)))
  
  # per region, get split and chimeric reads
  message("Amalgamating...")
  res <- mcmapply(function(r, re, sl, cl, hc, vc) {
    cl2 <- cl[sapply(cl, function(x, g) {
      y <- .bamDF2GR(x)
      any(from(findOverlaps(GRanges(y) , g)))
    }, GRanges(re))]
    sl2 <- sl[sapply(sl, function(x, g) {
      y <- .bamDF2GR(x)
      any(from(findOverlaps(GRanges(y) , g)))
    }, GRanges(re))]
    sre <- do.call(rbind.data.frame, sl2)
    cre <- do.call(rbind.data.frame, cl2)
    sdf <- data.frame(region = r,
                      split_read_evidence = nrow(sre),
                      chiral_read_evidence = length(unique(cre$qname)),
                      human_chr_split = paste0(unique(sre$rname), collapse = "|"),
                      viral_chr_split = paste0(unique(sre$virus.rname), collapse = "|"),
                      human_chr_chiral = paste0(unique(cre$rname[cre$rname %in% hc]), collapse = "|"),
                      viral_chr_chiral = paste0(unique(cre$rname[cre$rname %in% vc]), collapse = "|"))
    list(summary = sdf,
         region = r, 
         expanded.region = re,
         split.read.evidence = sre,
         chimeric.read.evidence = cre)
    }, paste0(sregions), paste0(sregions.expanded), 
    MoreArgs = list(cl = creads.filter2ByRead, sl = sl,
                  hc = hchrs, vc = vchrs), 
    mc.cores = detectCores()-2, SIMPLIFY = F)
  return(res)
}

output.dir <- "vFindR_VIcaller-example-output/"
res <- vFindRAmalgamate(output.dir)









