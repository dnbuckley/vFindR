setwd("~/Desktop/salhia_lab/vFindeR/")
library(Rsamtools)
library(ggplot2)
library(reshape2)
library(caret)
source("~/misc_R_scripts/samFlagTranslation.R")

bam.files <- list.files("fake_data/vFindR_output/", "bam$", full.names = T)
bai.files <- list.files("fake_data/vFindR_output/", "bai$", full.names = T)
param <- ScanBamParam(what = c("qname", "rname", "pos", "qwidth"))
bams <- mapply(scanBam, bam.files, bai.files, SIMPLIFY = T)
names(bams) <- gsub(".*\\/", "", bam.files)
x <- bams[[1]]
# bams <- bams[-grep("test_human.bam", names(bams))]

r1.swapped <- read.table("fake_data/R1_swapped.txt", col.names = "read.name")
r2.swapped <- read.table("fake_data/R2_swapped.txt", col.names = "read.name")
r1.swapped$read.name <- gsub("\\/.*", "", r1.swapped$read.name)
r2.swapped$read.name <- gsub("\\/.*", "", r2.swapped$read.name)
df.swap <- rbind.data.frame(data.frame(read.name = r1.swapped$read.name,
                                       read = "R1"),
                            data.frame(read.name = r2.swapped$read.name,
                                       read = "R2"))
df.swap <- cbind.data.frame(df.swap, do.call(cbind.data.frame,
                                             lapply(bams, function(b, n) {
                                               return(n %in% b$qname)
                                             }, df.swap$read.name)))
dfm <- melt(df.swap, id.vars = c("read.name", "read"))
dfm <- dfm[!dfm$variable %in% "test_human.bam", ]
dfm$value <- as.numeric(dfm$value)
ggplot(dfm, aes(x = variable, y = value, fill = read)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90))

df <- data.frame(read.name = bams$test_human.bam$qname,
                 read = sapply(bams$test_human.bam$flag, samFlagTranslation, check = "first_in_pair"))
df$R1.fake <- df$read.name %in% r1.swapped$read.name
df$R2.fake <- df$read.name %in% r2.swapped$read.name
df$both.fake <- df$R1.fake & df$R2.fake
df$first.unmapped.hg <- df$read.name %in% bams$test_umapped.first.bam$qname
df$second.unmapped.hg <- df$read.name %in% 
df$both.unmapped.hg <- df$read.name %in% bams$test_umapped.both.bam$qname
df$first.unmapped.vg <- df$read.name %in% bams$`test_unmapped-first-human_virus.bam`$qname
df$second.unmapped.vg <- df$read.name %in% bams$`test_unmapped-second-human_virus.bam`$qname
df$both.unmapped.vg <- df$read.name %in% bams$`test_unmapped-both-human_virus.bam`$qname


confusionMatrix(factor(df$first.unmapped.hg), factor(df$R1.fake))
confusionMatrix(factor(df$second.unmapped.hg), factor(df$R2.fake))
confusionMatrix(factor(df$both.unmapped.hg), factor(df$both.fake))




