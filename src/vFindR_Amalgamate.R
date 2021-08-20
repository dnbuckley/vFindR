setwd("~/Desktop/salhia_lab/vFindR/")
require(Rsamtools)
require(rtracklayer)
require(parallel)
require(dplyr)
source("src/split_read_table.R")
source("src/chimeric_reads_table.R")
source("src/vFindR_utils.R")
source("~/misc_R_scripts/combineInRange.R")


vFindRAmalgamate <- function(output.dir, 
                             combine.limit = 300,
                             min.split.reads = 5,
                             viruses = NULL,
                             repeatMask = T) {
  
  # limit to viruses of interest
  if (!is.null(viruses)) {
    if(!grepl("^NC_", viruses[1])) {
      viruses <- .species2vChr(viruses)
    }
  } 
  message("reading split reads")
  sreads <- .splitReadTable(output.dir, viruses = viruses)
  if (!is.null(viruses)) {
    vchrs <- as.character(unique(sreads$virus.rname))
  } else {
    vchrs <- as.character(unique(sreads$virus.rname))
  }
  sreads <- sreads[sreads$virus.rname %in% vchrs, ]
  sregions <- .bamDF2GR(sreads)
  sregions <- combineInRange(sregions, combine.limit)
  sregions.expanded <- GRanges(seqnames = seqnames(sregions), 
                               ranges = IRanges(start(sregions)-combine.limit,
                                                end = end(sregions)+combine.limit))
  
  hchrs <- as.character(unique(seqnames(sregions)))
  sl <- split(sreads, to(findOverlaps(.bamDF2GR(sreads), sregions.expanded)))
  idx <- sapply(sl, nrow) >= min.split.reads
  if (!any(idx)) {
    stop("no regions found with ", min.split.reads, " split read evidence.")
  }
  sl <- sl[idx]
  sregions <- sregions[idx]
  sregions.expanded <- sregions.expanded[idx]
  
  message("reading chimeric reads")
  creads <- .chimericReadTable(output.dir, do.filter = F)
  creads <- creads[!is.na(creads$pos), ]
  creads.gr <- .bamDF2GR(creads)
  creads.gr$qname <- creads$qname
  
  # chimeric reads with split evidence
  ovlp <- as.data.frame(findOverlaps(creads.gr, sregions.expanded))
  temp.l <- split(ovlp, ovlp$subjectHits)
  oll <- as.list(rep(NA, length(sregions.expanded)))
  oll[as.numeric(names(temp.l))] <- temp.l
  cl2 <- suppressWarnings(lapply(oll, function(o, x){
    if (is.na(o)) {
      return(NULL)
    }
    rnames <- x$qname[o$queryHits]
    return(return(x[x$qname %in% rnames, ]))
  }, creads))
  
  sreads.gr <- .bamDF2GR(sreads)
  ovlp <- findOverlaps(sreads.gr, sregions.expanded)
  oll <- split(ovlp, to(ovlp))
  sl2 <- mclapply(oll, function(o, x){
    rnames <- x$qname[from(o)]
    return(return(x[x$qname %in% rnames, ]))
  }, sreads, mc.cores = detectCores()-2)
  
  message("identified ", length(sregions), " regions, amalgamating split/chimeric read information...")
  res <- mapply(function(r, re, cre, sre, hc, vc) {
    sdf <- data.frame(region = r,
                      split_read_evidence = nrow(sre),
                      chimeric_read_evidence = length(unique(cre$qname)),
                      human_chr_split = paste0(unique(sre$rname), collapse = "|"),
                      viral_chr_split = paste0(unique(sre$virus.rname), collapse = "|"),
                      human_chr_chimeric = paste0(unique(cre$rname[cre$rname %in% hc]), collapse = "|"),
                      viral_chr_chimeric = paste0(unique(cre$rname[cre$rname %in% vc]), collapse = "|"))
    list(summary = sdf,
         region = r, 
         expanded.region = re,
         split.read.evidence = sre,
         chimeric.read.evidence = cre)
    }, paste0(sregions), paste0(sregions.expanded), cl2, sl2,
    MoreArgs = list(hc = hchrs, vc = vchrs), SIMPLIFY = F)
  # add info on whether the region is known as repetitive
  if (repeatMask) {
    rmask <- readRDS("~/Desktop/salhia_lab/vFindR/data/hg38_repeatMask.rds")
    ovlp <- findOverlaps(GRanges(unlist(lapply(res, function(x){x$region}))), rmask)
    for (i in 1:length(res)) {
      lst <- res[[i]]
      if (i %in% from(ovlp)) {
        lst$repeatMask <- rmask[to(ovlp)[from(ovlp) %in% i]]
        lst$summary$repeatRegion <- T
      } else {
        lst$repeatMask <- GRanges()
        lst$summary$repeatRegion <- F
      }
      res[[i]] <- lst
    }
  }
  return(res)
}





