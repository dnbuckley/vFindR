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
                             viruses = NULL) {
  
  # limit to viruses of interest
  if (!is.null(viruses)) {
    if(!grepl("^NC_", viruses[1])) {
      viruses <- .species2vChr(viruses)
    }
  } else message("WARNING: Not specifying viruses of interest will greatly increase runtime and memory usage.")
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
  creads <- .chimericReadTable(output.dir, do.filter = T)
  creads <- creads[!is.na(creads$pos), ]
  creads.gr <- .bamDF2GR(creads)
  creads.gr$qname <- creads$qname
  
  # chimeric reads with split evidence
  ovlp <- as.data.frame(findOverlaps(creads.gr, sregions.expanded))
  temp.l <- split(ovlp, ovlp$subjectHits)
  oll <- list()
  for (i in 1:length(sregions.expanded)) {
    if (any(i %in% as.numeric(names(temp.l)))) {
      oll[[i]] <- temp.l[[which(as.numeric(names(temp.l)) %in% i)]]
    } else {
      oll[[i]] <- NULL
    }
  }
  cl2 <- lapply(oll, function(o, x){
    if (is.null(o)) {
      return(NULL)
    }
    rnames <- x$qname[o$queryHits]
    return(return(x[x$qname %in% rnames, ]))
  }, creads)
  
  sreads.gr <- .bamDF2GR(sreads)
  ovlp <- findOverlaps(sreads.gr, sregions.expanded)
  oll <- split(ovlp, to(ovlp))
  sl2 <- lapply(oll, function(o, x){
    rnames <- x$qname[from(o)]
    return(return(x[x$qname %in% rnames, ]))
  }, sreads)
  
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
  return(res)
}

# output.dir <- "vFindR_VIcaller-example-output//"
output.dir <- "sample_data/vFindR_output/"
# viruses <- "Human_gammaherpesvirus_4"
res <- vFindRAmalgamate(output.dir, min.split.reads = 10)
res.summary <- do.call(rbind.data.frame, lapply(res, function(x){x$summary}))






