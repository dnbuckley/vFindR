setwd("~/Desktop/salhia_lab/vFindeR/")

if (Sys.info()["user"] == "dbuckley") {
  Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/Users/dbuckley/Desktop/bin", sep=":"))
  hg.idx.base <- "/Users/dbuckley/Desktop/salhia_lab/genomes/GRCh38_21/GRCh38_21.fa"
  R1 <- "/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/fake_R1.fastq.gz"
  R2 <- "/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/fake_R2.fastq.gz"
  output.dir <- "/Users/dbuckley/Desktop/salhia_lab/vFindeR/test_output/"
  output.stub <- paste0(output.dir, "test")
  threads <- "8"
  # first alignment
  # aln.hg.1.temp <- paste0(output.stub, "_human.bam")
  aln.hg.1.out <- paste0(output.stub, "_human.sam")
  aln.hg.1.sam.header <- paste0(output.stub, "_samheader.txt")
  # get first/second/unaligned pairs
  unmapped.first.hg.r1.out <- paste0(output.stub, "_human.unmapped.first.R1.sam")
  unmapped.first.hg.r2.out <- paste0(output.stub, "_human.unmapped.first.R2.sam")
  unmapped.second.hg.r1.out <- paste0(output.stub, "_human.unmapped.second.R1.sam")
  unmapped.second.hg.r2.out <- paste0(output.stub, "_human.unmapped.second.R2.sam")
  unmapped.both.hg.r1.out <- paste0(output.stub, "_human.unmapped.both.R1.sam")
  unmapped.both.hg.r2.out <- paste0(output.stub, "_human.unmapped.both.R2.sam")
  mapped.both.hg.out <- paste0(output.stub, "_human.mapped.sam")
  # convert to fastq for realignment
  unmapped.first.hg.r1.fastq <- paste0(output.stub, "_human.unmapped.first_R1.fastq.gz")
  unmapped.first.hg.r2.fastq <- paste0(output.stub, "_human.unmapped.first_R2.fastq.gz")
  unmapped.second.hg.r1.fastq <- paste0(output.stub, "_human.unmapped.second_R1.fastq.gz")
  unmapped.second.hg.r2.fastq <- paste0(output.stub, "_human.unmapped.second_R2.fastq.gz")
  unmapped.both.hg.r1.fastq <- paste0(output.stub, "_human.unmapped.both_R1.fastq.gz")
  unmapped.both.hg.r2.fastq <- paste0(output.stub, "_human.unmapped.both_R2.fastq.gz")
  
}

# first alignment
cmds <- vector()
cmds['clear'] <- paste0("rm -v /Users/dbuckley/Desktop/salhia_lab/vFindeR/test_output/*")
bt2.e <- system("which bowtie2", intern = T)
samtools.e <- system("which samtools", intern = T)
java.e <- system("which java", intern = T)
picard.jar <- "/Users/dbuckley/Desktop/bin/picard.jar"
aln.hg.1 <- paste(bt2.e, "-p", threads, 
                  "-x", hg.idx.base, 
                  "-1", R1, 
                  "-2", R2, " >", aln.hg.1.out)
cmds['update1'] <- paste0("echo Running first pass human alignment...")
cmds['aln.hg.1'] <- aln.hg.1
cmds['get.first.reads.r1'] <- paste(samtools.e, "view -f 69 -G 9 -o", unmapped.first.hg.r1.out, aln.hg.1.out)
cmds['get.first.reads.r2'] <- paste(samtools.e, "view -f 137 -G 5 -o", unmapped.first.hg.r2.out, aln.hg.1.out)
cmds['get.second.reads.r1'] <- paste(samtools.e, "view -f 133 -G 9 -o", unmapped.second.hg.r1.out, aln.hg.1.out)
cmds['get.second.reads.r2'] <- paste(samtools.e, "view -f 73 -G 5 -o", unmapped.second.hg.r2.out, aln.hg.1.out)
cmds['get.both.reads.r1'] <- paste(samtools.e, "view -f 77 -o", unmapped.both.hg.r1.out, aln.hg.1.out)
cmds['get.both.reads.r2'] <- paste(samtools.e, "view -f 141 -o", unmapped.both.hg.r2.out, aln.hg.1.out)
cmds['get.mapped'] <- paste(samtools.e, "view -F 4 -o", mapped.both.hg.out, aln.hg.1.out)
cmds['get.header'] <- paste("samtools view -H ", aln.hg.1.out, ">", aln.hg.1.sam.header)

# convert bam to fastq and realign to viral genome
cmds['temp'] <- paste("cat", aln.hg.1.sam.header, unmapped.first.hg.r1.out, unmapped.first.hg.r2.out, "> temp.sam")
cmds['picard.samtofastq.unmapped.first.R1'] <-
  paste0(java.e, " -jar ", picard.jar, " SamToFastq I=", "temp.sam", " F=", unmapped.first.hg.r1.fastq, " F2=", unmapped.first.hg.r2.fastq)


write.table(cmds, "commands.sh", row.names = F, col.names = F, quote = F)















