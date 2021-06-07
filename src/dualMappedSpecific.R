require(Rsamtools)
require(rtracklayer)

dualMappedSpecific <- function(path.to.vFindR.outdir, 
                      refseq.virus, 
                      samtools.e = NULL,
                      output.dir = NULL,
                      output.file = NULL,
                      make.depth.track = T) {
  if (is.null(samtools.e)) {
    samtools.e <- system("which samtools", intern = T)
  }
  # get reads mapped to specific viral chromosomes
  dual.bam.file <- list.files(path.to.vFindR.outdir, "dual.*bam$", full.names = T)
  first.bam.file <- list.files(path.to.vFindR.outdir, "first.*virus.bam$", full.names = T)
  second.bam.file <- list.files(path.to.vFindR.outdir, "second.*virus.bam$", full.names = T)
  stopifnot(length(dual.bam.file) == 1 &
              length(first.bam.file) == 1 &
              length(second.bam.file) == 1)
  dual.bam <- rtracklayer::import(dual.bam.file, 
                                  paired = F,
                                  use.names = T)
  first.bam <- rtracklayer::import(first.bam.file, 
                                   paired = F,
                                   use.names = T)
  second.bam <- rtracklayer::import(second.bam.file, 
                                    paired = F,
                                    use.names = T)
  reads.mapped.specific <- c(names(first.bam)[as.logical(seqnames(first.bam) %in% refseq.virus)],
                             names(second.bam)[as.logical(seqnames(second.bam) %in% refseq.virus)])
  if (length(reads.mapped.specific) == 0) {
    warning("No reads found.")
    return("Done.")
  }
  dual.bam.specific <- dual.bam[as.logical(names(dual.bam) %in% reads.mapped.specific)]
  output.dir <- ifelse(is.null(output.dir), ".", output.dir)
  output.file <- ifelse(is.null(output.file), gsub("dual", "specific", dual.bam.file), output.file)
  output.file <- gsub(".*\\/", "", output.file)
  output.file <- paste0(output.dir, "/", output.file)
  message("writing output to ", output.file)
  export(dual.bam.specific, output.file, "bam")
  system(paste("samtools index", output.file))
  if (make.depth.track) {
    depth.txt <- gsub("\\.bam$", ".depth.txt", output.file)
    depth.bg <- gsub("\\.bam$", ".depth.bedGraph", output.file)
    system(paste(samtools.e, "depth", output.file, ">", depth.txt))
    depth <- read.table(depth.txt)
    names(depth) <- c("seqnames", "start", "score")
    depth$end <- depth$start
    depth <- GenomicRanges::GRanges(depth)
    rtracklayer::export.bedGraph(depth, depth.bg)
    system(paste0("rm -v ", depth.txt))
  }
  return("Done.")
}








