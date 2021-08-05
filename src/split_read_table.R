require(Rsamtools)
require(rtracklayer)
require(parallel)
require(dplyr)
source("src/vFindR_utils.R")

.splitReadTable <- function(output.dir,
                            ref.species = "homo_sapiens",
                            virus.table = NULL) {
  virus.table <- openxlsx::read.xlsx("viral_seqnames.xlsx")
  ref.file <- list.files(output.dir, paste0(ref.species, "_dual-mapped.bam$"), full.names = T)
  virus.files <- list.files(paste0(output.dir, "/perVirus/"), "bam$", full.names = T)
  ref.df <- .readBAM(ref.file)
  ref.df <- ref.df[!is.na(ref.df$pos), ]
  virus.df <- do.call(rbind.data.frame, mclapply(virus.files, function(f) {
    d <- .readBAM(f)
    d <- d[!is.na(d$pos), ]
    names(d) <- paste0("virus.", names(d))
    d$refseq_virus <- gsub("(.*\\.REF_|\\.bed$)", "", f)
    return(d)
  }, mc.cores = detectCores()-2))
  merge.df <- ref.df %>% 
    dplyr::inner_join(virus.df, by = c("qname" = "virus.qname"))
  merge.df <- merge.df %>% 
    dplyr::left_join(virus.table, by = c("refseq_virus" = "seqnames"))
  merge.df$mapq.average <- (merge.df$mapq+merge.df$virus.mapq)/2
  return(merge.df)
}











