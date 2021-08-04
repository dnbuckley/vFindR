require(Rsamtools)
require(rtracklayer)
require(parallel)
require(dplyr)

.splitReadTable <- function(output.dir,
                            ref.species = "homo_sapiens",
                            virus.table = NULL) {
  virus.table <- openxlsx::read.xlsx("viral_seqnames.xlsx")
  ref.file <- list.files(output.dir, paste0(ref.species, "_dual-mapped.bed$"), full.names = T)
  virus.files <- list.files(paste0(output.dir, "/perVirus/"), "bed$", full.names = T)
  ref.gr <- import.bed(ref.file)
  ref.gr <- ref.gr[width(ref.gr) > 0]
  ref.df <- cbind.data.frame(read.location = paste0(granges(ref.gr)),
                             read.length = width(ref.gr),
                             read.name = ref.gr$name,
                             map.quality = ref.gr$score)
  virus.grs <- mclapply(virus.files, import.bed, mc.cores = detectCores()-2)
  virus.df <- do.call(rbind.data.frame, mapply(function(g, f) {
    g <- g[width(g) > 0]
    d <- cbind.data.frame(virus.read.location = paste0(granges(g)),
                          virus.read.length = width(g),
                          virus.read.name = g$name,
                          virus.map.quality = g$score)
    d$refseq_virus <- gsub("(.*\\.REF_|\\.bed$)", "", f)
    return(d)
  }, virus.grs, virus.files, SIMPLIFY = F))
  merge.df <- ref.df %>% 
    dplyr::inner_join(virus.df, by = c("read.name" = "virus.read.name"))
  merge.df <- merge.df %>% 
    dplyr::left_join(virus.table, by = c("refseq_virus" = "seqnames"))
  return(merge.df)
}