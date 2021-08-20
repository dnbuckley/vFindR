setwd("~/Desktop/salhia_lab/vFindR/")
source("src/vFindR_Amalgamate.R")
library(BSgenome.Hsapiens.UCSC.hg38)
library(Gviz)
library(Biostrings)

.readFromCigar <- function(seq, cigar) {
  spl <- strsplit(cigar, "[A-Z]")
  spl <- mapply(function(s, c) {
    s <- as.numeric(s)
    x <- strsplit(c, "[0-9]{1,10}")[[1]][-1]
    l <- unlist(mapply(function(s, x){rep(x, s)}, s, x))
    return(l)
  }, spl, cigar, SIMPLIFY = F)
  seq <- mapply(function(x, y) {
    x <- strsplit(x, "")[[1]]
    x[y != "M"] <- tolower(x[y != "M"])
    return(paste0(x, collapse = ""))
  }, seq, spl, SIMPLIFY = F)
  return(seq)
}

.formatChimericReads <- function(df) {
  byRead <- split(df, df$qname)
  # byRead <- byRead[sapply(byRead, nrow) == 2]
  seqs <- lapply(byRead, function(x) {
    unlist(.readFromCigar(x$seq, x$cigar))
  })
  return(seqs)
}

.getReadLoc <- function(gr, cigar) {
  spl <- strsplit(cigar, "[A-Z]")
  spl <- mapply(function(s, c) {
    s <- as.numeric(s)
    names(s) <- strsplit(c, "[0-9]{1,10}")[[1]][-1]
    return(s)
  }, spl, cigar, SIMPLIFY = F)
  starts <- mapply(function(c, s) {
    ifelse(names(c)[1] == "M", return(s), return(s - c[1]))
  }, spl, start(gr))
  ends <- starts+nchar(gr$seq)
  reads <- GRanges(seqnames = seqnames(gr), ranges = IRanges(start = starts,
                                                             end = ends-1))
  return(reads)
}

.inferBreakPoint <- function (gr, cigar) {
  spl <- strsplit(cigar, "[A-Z]")
  spl <- mapply(function(s, c) {
    s <- as.numeric(s)
    names(s) <- strsplit(c, "[0-9]{1,10}")[[1]][-1]
    return(s)
  }, spl, cigar, SIMPLIFY = F)
  starts <- mapply(function(c, s) {
    ifelse(names(c)[1] == "S", return(s), return(s + c[1]))
  }, spl, start(gr))
  bpoints <- GRanges(seqnames = seqnames(gr), ranges = IRanges(start = starts,
                                                               end = starts))
  return(bpoints)
}

.toDNAStringSet <- function(gr, seq) {
  x <- DNAStringSet(seq, start = start(gr), end = end(gr), width = width(gr), use.names = F)
  x <- DNAStringSet(seq, use.names = F)
}

.splitReads <- function(gr, cigar) {
  reads <- .getReadLoc(gr, cigar)
  spl <- strsplit(cigar, "[A-Z]")
  
  # # complicated cigars confuse Gviz
  # idx <- sapply(spl, function(x){length(x) > 3})
  # spl <- spl[!idx]
  # gr <- gr[!idx]
  # reads <- reads[!idx]
  # cigar <- cigar[!idx]
  
  spl <- mapply(function(s, c) {
    s <- as.numeric(s)
    names(s) <- strsplit(c, "[0-9]{1,10}")[[1]][-1]
    return(s)
  }, spl, cigar, SIMPLIFY = F)
  spl.reads <- mapply(function(c, s, e, chr, q, st) {
    s2 <- s
    stl <- vector()
    edl <- vector()
    # break up read, 2bp gap to account for GViz color fill weirdness
    for (i in 1:length(c)) {
      e2 <- s2+c[i]
      stl[i] <- s2
      edl[i] <- e2
      s2 <- e2 + 2
    }
    gspl <- GRanges(seqnames = chr, 
                    ranges = IRanges(start = stl, end = edl),
                    strand = st)
    gspl$mapping <- names(c)
    gspl$qname <- q
    return(gspl)
  }, spl, start(reads), end(reads), 
     as.character(seqnames(reads)), gr$qname, 
     as.character(strand(gr)), SIMPLIFY = F)
  spl.reads <- do.call("c", spl.reads)
  return(spl.reads)
}

.foo <- function(res, 
                 showSplitReads = F, 
                 expand = 0,
                 pdf = NULL) {
  reg <- GRanges(res$region)
  df.c <- res$chimeric.read.evidence
  df.s <- res$split.read.evidence
  options(ucscChromosomeNames=FALSE)
  df.c <- df.c[order(df.c$qname), ]
  # df.c$seq[df.c$strand == "-"] <- as.character(reverseComplement(DNAStringSet(df.c$seq[df.c$strand == "-"])))
  # df.c$strand <- "+"
  gr <- .bamDF2GR(df.c)
  strand(gr) <- df.c$strand
  # end(gr) <- df.c$pos+nchar(df.c$seq)
  seqs <- unlist(.formatChimericReads(df.c))
  gr$seq <- unlist(.formatChimericReads(df.c))
  gr$qname <- df.c$qname
  gr$cigar <- df.c$cigar
  gr.hg <- gr[seqnames(gr) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38)]
  gr.hg.reads <- .getReadLoc(gr.hg, gr.hg$cigar)
  gr.hg.reads$qname <- gr.hg$qname
  reads <- .splitReads(gr.hg, gr.hg$cigar)
  readsByGroup <- split(reads, reads$mapping)
  gr.s <- .bamDF2GR(df.s)
  gr.s$cigar <- df.s$cigar
  gr.s$seq <- df.s$seq
  gr.s$qname <- df.s$qname
  reads.s <- .splitReads(gr.s, gr.s$cigar)
  axisTrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(chromosome = unique(as.character(seqnames(gr.hg))), genome = "hg38")
  strack <- SequenceTrack(BSgenome.Hsapiens.UCSC.hg38)
  btrack <- AnnotationTrack(reduce(.inferBreakPoint(gr.hg, gr.hg$cigar)), fill = "red")
  # t1 <- AnnotationTrack(readsByGroup$M, group = readsByGroup$M$qname, genome = "hg38", fill = "blue", stacking = "full")
  # t2 <- AnnotationTrack(readsByGroup$S, group = readsByGroup$S$qname, genome = "hg38", fill = "red", stacking = "full")
  t3 <- AnnotationTrack(reads.s, group, genome = "hg38")
  # atrack <- OverlayTrack(c(t1, t2))
  atrack <- AnnotationTrack(reads, 
                            fill = ifelse(reads$mapping == "M", "blue",
                                          ifelse(reads$mapping == "S", "red", "gray25")),
                            # feature = reads$mapping, 
                            # group = reads$qname,
                            strand = "*")
  if (showSplitReads) {
    tracks <- c(itrack, axisTrack, strack, btrack, atrack, t3)
  } else {
    tracks <- c(itrack, axisTrack, strack, btrack, atrack)
  }
  if (is.null(pdf)) {
    pdf <- paste0(reg, ".pdf")
  }
  pdf(pdf, height = 14, width = 11)
  plotTracks(tracks, showOverplotting = TRUE, 
             chromosome = seqnames(reg), 
             from = min(start(reads))-expand, 
             to = max(end(reads))+expand)
  dev.off()
}

res.l <- vFindRAmalgamate("sample_data/CNSMets_RNAseq/J00255_C3TRYACXX/vFindR_output/", 
                          combine.limit = 20, repeatMask = T)
summary <- do.call(rbind.data.frame, lapply(res.l, function(x){x$summary}))
# summary <- summary[grepl("007605", summary$viral_chr_split) | grepl("007605", summary$viral_chr_chimeric), ]
summary <- summary %>% 
  dplyr::filter(split_read_evidence > 10 &
                  chimeric_read_evidence > 10)
summary <- summary[order(summary$split_read_evidence+summary$chimeric_read_evidence, decreasing = T), ]
res <- res.l$`chr4:54727404-54727475`

.foo(res, expand = 0, showSplitReads = F)











