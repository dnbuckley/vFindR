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
  # get first/second/unaligned pairs
  unmapped.first.hg.r1.out <- paste0(output.stub, "_human.unmapped.first.R1.sam")
  unmapped.first.hg.r2.out <- paste0(output.stub, "_human.unmapped.first.R2.sam")
  unmapped.second.hg.r1.out <- paste0(output.stub, "_human.unmapped.second.R1.sam")
  unmapped.second.hg.r2.out <- paste0(output.stub, "_human.unmapped.second.R2.sam")
  unmapped.both.hg.r1.out <- paste0(output.stub, "_human.unmapped.both.R1.sam")
  unmapped.both.hg.r2.out <- paste0(output.stub, "_human.unmapped.both.R2.sam")
  mapped.both.hg.out <- paste0(output.stub, "_human.mapped.sam")
}

cmds <- vector()
cmds['clear'] <- paste0("rm -v /Users/dbuckley/Desktop/salhia_lab/vFindeR/test_output/*")
bt2.e <- system("which bowtie2", intern = T)
samtools.e <- system("which samtools", intern = T)
# ALIGN="$BOWTIE2 -p $THREADS -x $GENOMEDIR/$GENOME.fa --met 60 --met-file $STUB.bt2.log ${RGFLAG} -1 $R1 -2 $R2"
aln.hg.1 <- paste(bt2.e, "-p", threads, 
                  "-x", hg.idx.base, 
                  "-1", R1, 
                  "-2", R2, " >", aln.hg.1.out)
cmds['update1'] <- paste0("echo Running first pass human alignment...")
cmds['aln.hg.1'] <- aln.hg.1
# cmds['aln.hg.1.sort'] <- paste(samtools.e, "sort -n -@", threads, "-o", aln.hg.1.out, aln.hg.1.temp)
# cmds['aln.hg.1.rmtemp'] <- paste("rm -v", aln.hg.1.temp)
# get left/right/unaligned pairs
# cmds['get.first.reads.r1'] <- paste(samtools.e, "view -f 69 -o", unmapped.first.hg.r1.out, aln.hg.1.out)
# cmds['get.first.reads.r2'] <- paste(samtools.e, "view -f 137 -o", unmapped.first.hg.r2.out, aln.hg.1.out)
# cmds['get.second.reads.r1'] <- paste(samtools.e, "view -f 133 -o", unmapped.second.hg.r1.out, aln.hg.1.out)
# cmds['get.second.reads.r2'] <- paste(samtools.e, "view -f 73 -o", unmapped.second.hg.r2.out, aln.hg.1.out)
# cmds['get.both.reads.r1'] <- paste(samtools.e, "view -f 77 -o", unmapped.both.hg.r1.out, aln.hg.1.out)
# cmds['get.both.reads.r2'] <- paste(samtools.e, "view -f 141 -o", unmapped.both.hg.r2.out, aln.hg.1.out)
cmds['get.first.reads.r1'] <- paste(samtools.e, "view -f 69 -G 9 -o", unmapped.first.hg.r1.out, aln.hg.1.out)
cmds['get.first.reads.r2'] <- paste(samtools.e, "view -f 137 -G 5 -o", unmapped.first.hg.r2.out, aln.hg.1.out)
cmds['get.second.reads.r1'] <- paste(samtools.e, "view -f 133 -G 9 -o", unmapped.second.hg.r1.out, aln.hg.1.out)
cmds['get.second.reads.r2'] <- paste(samtools.e, "view -f 73 -G 5 -o", unmapped.second.hg.r2.out, aln.hg.1.out)
cmds['get.both.reads.r1'] <- paste(samtools.e, "view -f 77 -o", unmapped.both.hg.r1.out, aln.hg.1.out)
cmds['get.both.reads.r2'] <- paste(samtools.e, "view -f 141 -o", unmapped.both.hg.r2.out, aln.hg.1.out)
cmds['get.mapped'] <- paste(samtools.e, "view -f 141 -o", unmapped.both.hg.r2.out, aln.hg.1.out)
cmds['get.mapped'] <- paste(samtools.e, "view -F 4 -o", mapped.both.hg.out, aln.hg.1.out)


write.table(cmds, "commands.sh", row.names = F, col.names = F, quote = F)















