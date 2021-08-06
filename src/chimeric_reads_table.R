require(Rsamtools)
require(rtracklayer)
require(parallel)
require(dplyr)
source("src/vFindR_utils.R")

# shoud make this default in the pipleine
.filterLocalBamFile <- function(bam.file) {
  filtered.file <- gsub("\\.bam$", ".filtered.bam", bam.file)
  cmd <- paste0("samtools view -F 4 -O SAM -h ", bam.file , 
                " | /Users/dbuckley/Desktop/salhia_lab/vFindR/src/filterChimericBam | ",
                "samtools sort -O BAM > ", filtered.file)
  message(cmd)
  system(cmd)
  cmd <- paste("samtools index", filtered.file)
  message(cmd)
  system(cmd)
}

# need to rewrite pipeline to preprocess for easy reading
.chimericReadTable <- function(output.dir,
                               ref.species = "homo_sapiens",
                               virus.table = NULL,
                               do.filter = T) {
  bam.files <- list.files(output.dir, "localalign.bam$", full.names = T)
  if (do.filter) {
    # don't redo if already done
    allBams <- list.files(output.dir, ".bam$", full.names = T)
    if (!all(gsub("\\.bam$", ".filtered.bam", bam.files) %in% allBams)) {
      mclapply(bam.files, .filterLocalBamFile, mc.cores = detectCores()-2)
    }
    bam.files <- list.files(output.dir, "localalign.filtered.bam$", full.names = T)
  }
  bams <- mclapply(bam.files, function(f){
    d <- .readBAM(f)
    d <- d[d$qname %in% names(table(d$qname))[table(d$qname)>1], ]
    return(d)
  }, mc.cores = detectCores()-2)
  read.tab <- do.call(rbind.data.frame, bams)
  return(read.tab)
}