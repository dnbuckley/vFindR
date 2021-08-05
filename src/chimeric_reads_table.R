require(Rsamtools)
require(rtracklayer)
require(parallel)
require(dplyr)
source("src/vFindR_utils.R")

# need to rewrite pipeline to preprocess for easy reading
.chimericReadTable <- function(output.dir,
                               ref.species = "homo_sapiens",
                               virus.table = NULL) {
  bam.files <- list.files(output.dir, "localalign.bam$", full.names = T)
  bams <- mclapply(bam.files, function(f){
    d <- .readBAM(f)
    d <- d[d$qname %in% names(table(d$qname))[table(d$qname)>1], ]
    return(d)
  }, mc.cores = detectCores()-2)
  read.tab <- do.call(rbind.data.frame, bams)
  return(read.tab)
}