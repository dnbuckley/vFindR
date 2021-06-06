# main function

vFindR <- function(sample.dir = NULL,
                   ref.genome.idx, vir.genome.idx, 
                   path.to.picard.jar = NULL,
                   R1 = NULL,
                   R2 = NULL,
                   output.name = NULL,
                   bt2.e = NULL,
                   samtools.e = NULL,
                   java.e = NULL,
                   ref.species = "homo_sapiens",
                   threads = 8,
                   output.dir = "vFindR_output",
                   mode = "sh",
                   slurm.header.args = c("-t 7-00:00:00", 
                                         "-p dtg", 
                                         "-A davidwcr_263", 
                                         "-c 16", 
                                         "--mem=55G",
                                         "-J vFindR", 
                                         "-o %x\\_%j.out")) {
  # some idiot proofing
  if ((!is.null(sample.dir) & !is.null(R1)) |
      (!is.null(sample.dir) & !is.null(R2)) |
      (is.null(sample.dir) & (is.null(R1) | is.null(R2)))) {
    stop("Specify directory with reads (gzipped) OR R1 & R2")
  }
  # get paths and files
  if (Sys.info()["user"] == "dbuckley") {
    Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/Users/dbuckley/Desktop/bin", sep=":"))
  }
  if (is.null(bt2.e)) {
    bt2.e <- system("which bowtie2", intern = T)
  } 
  if (is.null(samtools.e)) {
    samtools.e <- system("which samtools", intern = T)
  }
  if (is.null(java.e)) {
    java.e <- system("which java", intern = T)
  } 
  if (is.null(path.to.picard.jar)) {
    path.to.picard.jar <- system("which picard.jar", intern = T)
  }
  if (is.null(R1) & is.null(R2)) {
    R1 <- normalizePath(list.files(sample.dir, "R1.*fastq.gz", full.names = T))
    R2 <- normalizePath(list.files(sample.dir, "R2.*fastq.gz", full.names = T))
    output.dir <- paste0(normalizePath(sample.dir), "/", output.dir)
  }
  # some defaults
  if(is.null(output.name)) {
    output.name <- "vFindR"
  }
  output.stub <- paste0(output.dir, "/", output.name)
  
  # set up all the files
  aln.hg.1.out <- paste0(output.stub, "_", ref.species, ".bam")
  aln.hg.1.sam.header <- paste0(output.stub, "_samheader.txt")
  # get first/second/unaligned pairs
  unmapped.first.bam <- paste0(output.stub, "_", ref.species, ".umapped.first.bam")
  unmapped.second.bam <- paste0(output.stub, "_", ref.species, ".umapped.second.bam")
  unmapped.both.bam <- paste0(output.stub, "_", ref.species, ".umapped.both.bam")
  # convert to fastq for realignment
  unmapped.first.hg.r1.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.first_R1.fastq.gz")
  unmapped.first.hg.r2.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.first_R2.fastq.gz")
  unmapped.second.hg.r1.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.second_R1.fastq.gz")
  unmapped.second.hg.r2.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.second_R2.fastq.gz")
  unmapped.both.hg.r1.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.both_R1.fastq.gz")
  unmapped.both.hg.r2.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.both_R2.fastq.gz")
  aln.vir.first.bam <- paste0(output.stub, "_unmapped-first-", ref.species, "_virus.bam")
  aln.vir.second.bam <- paste0(output.stub, "_unmapped-second-", ref.species, "_virus.bam")
  aln.vir.both.bam <- paste0(output.stub, "_unmapped-both-", ref.species, "_virus.bam")
  
  # wirite all the commands
  cmds <- vector()
  cmds['make.output.dir'] <- paste0("mkdir ", output.dir)
  cmds['change.dir'] <- paste0("cd ", output.dir)
  cmds['aln.hg.1'] <- paste(bt2.e, "-p", threads, 
                            "-x", ref.genome.idx, 
                            "-1", R1, 
                            "-2", R2, " | samtools sort -O BAM -@", threads, ">", aln.hg.1.out)
  
  cmds['get.unmapped.first'] <- paste0("{ ", 
                                       paste(samtools.e, "view -H ", aln.hg.1.out), " && ",
                                       paste(samtools.e, "view -f 69 -G 9", aln.hg.1.out), " && ",
                                       paste(samtools.e, "view -f 137 -G 5", aln.hg.1.out), 
                                       "; } | samtools sort -O BAM -@ ", threads, " > ", unmapped.first.bam)
  cmds['get.unmapped.second'] <- paste0("{ ", 
                                        paste(samtools.e, "view -H ", aln.hg.1.out), " && ",
                                        paste(samtools.e, "view -f 133 -G 9", aln.hg.1.out), " && ",
                                        paste(samtools.e, "view -f 73 -G 5", aln.hg.1.out), 
                                        "; } | samtools sort -O BAM -@ ", threads, " > ", unmapped.second.bam)
  cmds['get.unmapped.both'] <- paste0("{ ", 
                                      paste(samtools.e, "view -H ", aln.hg.1.out), " && ",
                                      paste(samtools.e, "view -f 77", aln.hg.1.out), " && ",
                                      paste(samtools.e, "view -f 141", aln.hg.1.out), 
                                      "; } | samtools sort -O BAM -@ ", threads, " > ", unmapped.both.bam)
  
  # # convert bam to fastq and realign to viral genome
  cmds['picard.samtofastq.unmapped.first.R1'] <-
    paste0(java.e, " -jar ", path.to.picard.jar, " SamToFastq I=", unmapped.first.bam, 
           " F=", unmapped.first.hg.r1.fastq, 
           " F2=", unmapped.first.hg.r2.fastq)
  cmds['picard.samtofastq.unmapped.second.R1'] <-
    paste0(java.e, " -jar ", path.to.picard.jar, " SamToFastq I=", unmapped.second.bam, 
           " F=", unmapped.second.hg.r1.fastq, 
           " F2=", unmapped.second.hg.r2.fastq)
  cmds['picard.samtofastq.unmapped.both.R1'] <-
    paste0(java.e, " -jar ", path.to.picard.jar, " SamToFastq I=", unmapped.both.bam, 
           " F=", unmapped.both.hg.r1.fastq, 
           " F2=", unmapped.both.hg.r2.fastq)
  
  cmds['aln.vir.first'] <- paste(bt2.e, "--very-sensitive",
                                 "-p", threads, 
                                 "-x", vir.genome.idx, 
                                 "-1", unmapped.first.hg.r1.fastq, 
                                 "-2", unmapped.first.hg.r2.fastq, 
                                 " | samtools sort -O BAM -@", threads, ">", aln.vir.first.bam)
  cmds['aln.vir.second'] <- paste(bt2.e, "--very-sensitive",
                                  "-p", threads, 
                                  "-x", vir.genome.idx, 
                                  "-1", unmapped.second.hg.r1.fastq, 
                                  "-2", unmapped.second.hg.r2.fastq, 
                                  " | samtools sort -O BAM -@", threads, ">", aln.vir.second.bam)
  cmds['aln.vir.both'] <- paste(bt2.e, "--very-sensitive",
                                "-p", threads, 
                                "-x", vir.genome.idx, 
                                "-1", unmapped.both.hg.r1.fastq, 
                                "-2", unmapped.both.hg.r2.fastq, 
                                " | samtools sort -O BAM -@", threads, ">", aln.vir.both.bam)
  cmds['index'] <- paste0("for i in ", output.dir, "/*bam; do ", samtools.e, " index $i; done")
  cmds['idxstats'] <- paste0("for i in ", output.dir, "/*bam; do ", samtools.e, " idxstats $i > $i.idxstats; done")
  cmds <- .addEcho(cmds)
  if (mode == "local") {
      lapply(cmds, system, ignore.stdout = F, ignore.stderr = F)
  } else if (mode == "sh") {
    write.table(cmds, paste0(output.name, "_commands.sh"),
                row.names = F, col.names = F, quote = F)
  } else if (mode == "slurm") {
    slurm.header.args <- paste0("#SBATCH ", slurm.header.args)
    cmds <- c("#!/bin/bash", slurm.header.args, "\n\n\n", cmds)
    write.table(cmds, paste0(output.name, "_commands.slurm"),
                row.names = F, col.names = F, quote = F)
  } else {
    stop("Mode must be one of: 'local', 'sh', or 'slurm'")
  }
  return("Done.")
}

.addEcho <- function(cmds) {
  c <- vector()
  n <- 1
  for (i in 1:length(cmds)) {
    c[n] <- "printf \"\\n######## COMMAND: ########\\n\""
    c[n+1] <- paste0("echo ", "\"", cmds[i], "\"")
    c[n+2] <- "printf \"##########################\\n\""
    c[n+3] <- cmds[i]
    n <- n+4
  }
  return(c)
}










