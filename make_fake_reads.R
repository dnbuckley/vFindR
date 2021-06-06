setwd("~/Desktop/salhia_lab/vFindeR/")
library(seqinr)
library(ShortRead)

sampleString = function(s, l) {
  nStart = sample(1:(nchar(s)-l), 1)
  substr(s, nStart, nStart+l-1)
}


subRead.single <- function(r, f, n = 1000) {
  rlen <- median(nchar(sread(r)))
  false.reads <- replicate(n, sampleString(as.character(sread(f)), rlen))
  stopifnot(all(nchar(false.reads) == rlen))
  idx <- sample(1:length(r), n, replace = F)
  false.reads <- DNAStringSet(false.reads)
  r@sread[idx] <- false.reads
  swapped.reads <- r@id[idx]
  return(list(read = r, swapped.reads = swapped.reads))
}

shrink <- T

r1 <- readFastq("sample_data/HG00099.EWS-ctrl.hg19.sorted.mkdup.chr21_R1.fastq.gz")
r2 <- readFastq("sample_data/HG00099.EWS-ctrl.hg19.sorted.mkdup.chr21_R2.fastq.gz")
fa <- readFasta("~/Desktop/salhia_lab/genomes/ebv/ebv.fa")

if (shrink) {
  idx <- sample(1:length(r1), 40000, replace = F)
  r1 <- r1[idx]
  r2 <- r2[idx]
}

s1 <- subRead.single(r1, fa)
s2 <- subRead.single(r2, fa)
system("rm -v fake_data/*")
writeFastq(s1$read, "fake_data/fake_R1.fastq", compress = F)
writeFastq(s2$read, "fake_data/fake_R2.fastq", compress = F)
writeFastq(s1$read, "fake_data/fake_R1.fastq.gz", compress = T)
writeFastq(s2$read, "fake_data/fake_R2.fastq.gz", compress = T)
write.table(as.character(s1$swapped.reads), "fake_data/R1_swapped.txt",
            row.names = F, col.names = F, quote = F)
write.table(as.character(s2$swapped.reads), "fake_data/R2_swapped.txt",
            row.names = F, col.names = F, quote = F)
write.table(as.character(r1@id), "fake_data/R1_names.txt",
            row.names = F, col.names = F, quote = F)
write.table(as.character(r1@id), "fake_data/R2_names.txt",
            row.names = F, col.names = F, quote = F)





