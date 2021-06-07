setwd("~/Desktop/salhia_lab/genomes/refseq/viral/")
library(ballgown)
library(tidyverse)
library(seqinr)
library(ShortRead)

getFileDF <- function(files) {
  file.df <- data.frame(ID = gsub("^\\.\\/", "", files),
                        file = files)
  file.df$ID <- gsub("\\/.*", "", file.df$ID)
  return(file.df)
}

fasta.files <- list.files(".", "fna", recursive = T, full.names = T)
fasta.files <- fasta.files[!grepl("(cds|rna)_from_genomic", fasta.files)]
species.files <- list.files(".", "_organism", recursive = T, full.names = T)
gtf.files <- list.files(".", "gff", recursive = T, full.names = T)

fasta.files <- getFileDF(fasta.files)
names(fasta.files)[2] <- "fasta.file"
species.files <- getFileDF(species.files)
gtf.files <- getFileDF(gtf.files)
names(gtf.files)[2] <- "gtf.file"

species.files$species <- mclapply(species.files$file, function(x) {
  x <- read.table(x, comment.char = "")
  x <- paste0(x[, 4:(ncol(x)-1)], collapse = "_")
  return(x)
}, mc.cores = 10)

df <- fasta.files %>% 
  left_join(gtf.files, by = "ID") %>% 
  left_join(select(species.files, ID, species), by = "ID")

df <- df[rowSums(is.na(df)) == 0, ]

gtfs <- mclapply(df$gtf.file, gffRead, mc.cores = 10)
chrs <- unlist(sapply(gtfs, function(x){unique(x[, 1])}))
length(unique(chrs)) == length(chrs)

fastas <- mclapply(df$fasta.file, readFasta, mc.cores = 6)
fastas <- lapply(fastas, function(x) {
  x@id <- BStringSet(gsub(" .*", "", x@id))
  return(x)
})
ids <- unlist(sapply(fastas, function(x) {as.character(x@id)}))
length(unique(ids)) == length(ids)

df$species <- unlist(df$species)
write.table(df, "../file_table.txt", sep = "\t", 
            row.names = F, col.names = F, quote = F)

df.chr <- do.call(rbind, mapply(function(f, s) {
  data.frame(seqnames = unique(as.character(f@id)),
             organism = s)
}, fastas, df$species,  SIMPLIFY = F))

write.table(df.chr, "../../../vFindeR/viral_seqnames.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)
openxlsx::write.xlsx(df.chr, "../../../vFindR/viral_seqnames.xlsx")

fa <- fastas[[1]]
for (i in 2:length(fastas)) {
  message(i)
  f <- fastas[[i]]
  fa <- append(fa, f)
}
gtf <- do.call(rbind, gtfs)
writeFasta(fa, "../../refseq_all_virus/refseq_all_virus.fa")
write.table(gtf, "../../refseq_all_virus/refseq_all_virus.gtf", 
            row.names = F, col.names = F, quote = F, sep = "\t")











