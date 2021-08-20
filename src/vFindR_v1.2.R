# main function

vFindR <- function(sample.dir = NULL,
                   ref.genome.idx, vir.genome.idx, ref_vir.genome.idx,
                   path.to.picard.jar = NULL,
                   R1 = NULL,
                   R2 = NULL,
                   output.name = NULL,
                   bt2.e = NULL,
                   samtools.e = NULL,
                   bamtools.e = NULL,
                   java.e = NULL,
                   path.to.extract.py = NULL,
                   python.e = NULL,
                   ref.species = "homo_sapiens",
                   threads = 8,
                   output.dir = "vFindR_output",
                   mode = "sh",
                   vFindR.dir = NULL,
                   cleanup = T,
                   slurm.header.args = c("-t 7-00:00:00", 
                                         "-p dtg", 
                                         "-A davidwcr_263", 
                                         "-c 16", 
                                         "--mem=55G",
                                         "-J vFindR", 
                                         "-o %x\\_%j.out")) {
  # some checks
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
  if (is.null(bamtools.e)) {
    bamtools.e <- system("which bamtools", intern = T)
  }
  if (is.null(java.e)) {
    java.e <- system("which java", intern = T)
  } 
  if (is.null(path.to.picard.jar)) {
    path.to.picard.jar <- system("which picard.jar", intern = T)
  }
  if (is.null(python.e)){
    python.e <- system("which python3", intern = T)
  }
  if (is.null(path.to.extract.py)) {
    path.to.extract.py <- "~/bin/extract_reads.py"
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
  output.stub.perVirus <- paste0(output.dir, "/", "perVirus/", output.name)
  filterChimericBam <- paste0(vFindR.dir, "/src/filterChimericBam")
  
  # set up all the files
  aln.hg.1.bam <- paste0(output.stub, "_", ref.species, ".bam")
  aln.hg.1.sam.header <- paste0(output.stub, "_samheader.txt")
  # files for first/second/unaligned pairs
  unmapped.first.bam <- paste0(output.stub, "_", ref.species, ".unmapped.first.bam")
  unmapped.second.bam <- paste0(output.stub, "_", ref.species, ".unmapped.second.bam")
  unmapped.both.bam <- paste0(output.stub, "_", ref.species, ".unmapped.both.bam")
  # files for convert to fastq for realignment
  unmapped.first.hg.r1.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.first_R1.fastq.gz")
  unmapped.first.hg.r2.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.first_R2.fastq.gz")
  unmapped.second.hg.r1.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.second_R1.fastq.gz")
  unmapped.second.hg.r2.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.second_R2.fastq.gz")
  unmapped.both.hg.r1.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.both_R1.fastq.gz")
  unmapped.both.hg.r2.fastq <- paste0(output.stub, "_", ref.species, ".unmapped.both_R2.fastq.gz")
  # files for align to viral genome
  aln.vir.first.bam <- paste0(output.stub, "_unmapped-first-", ref.species, "_virus.bam")
  aln.vir.second.bam <- paste0(output.stub, "_unmapped-second-", ref.species, "_virus.bam")
  aln.vir.both.bam <- paste0(output.stub, "_unmapped-both-", ref.species, "_virus.bam")
  # aln.vir.first.mapped.readnames <- paste0(output.stub, "_", ref.species, "_first-remapped-to-virus_readnames.txt")
  # aln.vir.second.mapped.readnames <- paste0(output.stub, "_", ref.species, "_second-remapped-to-virus_readnames.txt")
  # realn.mapped.first.bam <- paste0(output.stub, "_", ref.species, "_first-remapped-to-virus_", ref.species, ".bam")
  # realn.mapped.second.bam <- paste0(output.stub, "_", ref.species, "_second-remapped-to-virus_", ref.species, ".bam")  
  dual.mapped.readnames <- paste0(output.stub, "_", ref.species, "_dual-mapped_readnames.txt")
  dual.mapped.temp <- paste0(output.stub, "_", ref.species, "_dual-mapped_temp")
  dual.mapped.bam <- paste0(output.stub, "_", ref.species, "_dual-mapped.bam")
  aln.vir.first.perVirus.stub <- paste0(output.stub.perVirus, "_unmapped-first-", ref.species, "_virus")
  aln.vir.second.perVirus.stub <- paste0(output.stub.perVirus, "_unmapped-second-", ref.species, "_virus")
  
  # potential.chimeric.reads.header <- paste0(output.stub, "_potential_chimeric.samHeader")
  # potential.chimeric.reads.sam <- paste0(output.stub, "_potential_chimeric.sam")
  # potential.chimeric.reads.bam <- paste0(output.stub, "_potential_chimeric.bam")
  potential.chimeric.reads.first.fastq <- paste0(output.stub, "_potential_chimeric_first.fastq.gz")
  potential.chimeric.reads.second.fastq <- paste0(output.stub, "_potential_chimeric_second.fastq.gz")
  potential.chimeric.reads.both.fastq <- paste0(output.stub, "_potential_chimeric_both.fastq.gz")
  local.first.bam <- paste0(output.stub, "_first_localalign.bam")
  local.second.bam <- paste0(output.stub, "_second_localalign.bam")
  local.both.bam <- paste0(output.stub, "_both_localalign.bam")
  
  
  
  
  
  # wirite all the commands
  cmds <- vector()
  cmds['make.output.dir'] <- paste0("mkdir ", output.dir)
  cmds['make.output.dir.perVirus'] <- paste0("mkdir ", output.dir, "/", "perVirus")
  cmds['change.dir'] <- paste0("cd ", output.dir)
  cmds['aln.hg.1'] <- paste(bt2.e, "-p", threads, 
                            "-x", ref.genome.idx, 
                            "-1", R1, 
                            "-2", R2, " | samtools sort -n -O BAM >", aln.hg.1.bam)
  
  cmds['get.unmapped.first'] <- paste0("{ ", 
                                       paste(samtools.e, "view -H ", aln.hg.1.bam), " && ",
                                       paste(samtools.e, "view -f 69 -G 9", aln.hg.1.bam), " && ",
                                       paste(samtools.e, "view -f 137 -G 5", aln.hg.1.bam), 
                                       "; } | samtools sort -O BAM > ", unmapped.first.bam)
  cmds['get.unmapped.second'] <- paste0("{ ", 
                                        paste(samtools.e, "view -H ", aln.hg.1.bam), " && ",
                                        paste(samtools.e, "view -f 133 -G 9", aln.hg.1.bam), " && ",
                                        paste(samtools.e, "view -f 73 -G 5", aln.hg.1.bam), 
                                        "; } | samtools sort -O BAM > ", unmapped.second.bam)
  cmds['get.unmapped.both'] <- paste0("{ ", 
                                      paste(samtools.e, "view -H ", aln.hg.1.bam), " && ",
                                      paste(samtools.e, "view -f 77", aln.hg.1.bam), " && ",
                                      paste(samtools.e, "view -f 141", aln.hg.1.bam), 
                                      "; } | samtools sort -O BAM > ", unmapped.both.bam)
  
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
  # align to viral genome
  cmds['aln.vir.first'] <- paste(bt2.e, "--very-sensitive",
                                 "-p", threads, 
                                 "-x", vir.genome.idx, 
                                 "-1", unmapped.first.hg.r1.fastq, 
                                 "-2", unmapped.first.hg.r2.fastq, 
                                 " | samtools sort -O BAM >", aln.vir.first.bam)
  cmds['aln.vir.second'] <- paste(bt2.e, "--very-sensitive",
                                  "-p", threads, 
                                  "-x", vir.genome.idx, 
                                  "-1", unmapped.second.hg.r1.fastq, 
                                  "-2", unmapped.second.hg.r2.fastq, 
                                  " | samtools sort -O BAM >", aln.vir.second.bam)
  cmds['aln.vir.both'] <- paste(bt2.e, "--very-sensitive",
                                "-p", threads, 
                                "-x", vir.genome.idx, 
                                "-1", unmapped.both.hg.r1.fastq, 
                                "-2", unmapped.both.hg.r2.fastq, 
                                " | samtools sort -O BAM >", aln.vir.both.bam)
  
  # # extract reads that mapped discordantly in human and did map to virus
  cmds['index.first.bam'] <- paste(samtools.e, "index", aln.hg.1.bam)
  cmds['get.dual.mapped.readnames'] <- paste0("{ ", 
                                              paste(samtools.e, "view -F 4", aln.vir.first.bam), " && ",
                                              paste(samtools.e, "view -F 4", aln.vir.second.bam),
                                              "; } | cut -f 1 | sort > ", dual.mapped.readnames)
  # cmds['extract.dual.mapped.reads'] <- paste(python.e, path.to.extract.py, "-b", aln.hg.1.bam, "-n", dual.mapped.readnames,
  #                                            "-o /dev/stdout |", samtools.e,"sort -O BAM -@", threads, ">", dual.mapped.bam)
  # cmds['extract.dual.mapped.reads'] <- paste(python.e, path.to.extract.py, "-b", aln.hg.1.bam, "-n", dual.mapped.readnames, 
  #                                            "-o /dev/stdout | samtools sort -O BAM -@", threads, ">", dual.mapped.bam)
  cmds['extract.dual.mapped.reads'] <- paste(python.e, path.to.extract.py, "-b", aln.hg.1.bam, "-n", dual.mapped.readnames,
                                             "-o", dual.mapped.temp)
  cmds['convert.dual.mapped.to.bam'] <- paste(samtools.e,"sort -O BAM -@", threads, dual.mapped.temp, ">", dual.mapped.bam)
  cmds['remove.dual.mapped.temp'] <- paste("rm -v", dual.mapped.temp)
  cmds['split.perVirus.first'] <- paste(bamtools.e, "split -reference -stub", aln.vir.first.perVirus.stub, "-in", aln.vir.first.bam)
  cmds['split.perVirus.second'] <- paste(bamtools.e, "split -reference -stub", aln.vir.second.perVirus.stub, "-in", aln.vir.second.bam)
  cmds['remove.unmapped.bams'] <- paste0("rm -vf ", output.dir, "/perVirus/*_unmapped.bam")
  # cmds['dual.bam.to.bed'] <- paste(bamtools.e, "convert -format bed -in", dual.mapped.bam, "-out", gsub("bam$", "bed", dual.mapped.bam))
  # cmds['viral.bams.to.bed'] <- paste0("bash ", vFindR.dir, "/src/", "convert_all_virus_bams.sh ", 
  #                                     paste0(output.dir, "/", "perVirus/"), " ", bamtools.e)
  
  # take everything that remains unmapped; potential chimeric reads
  # get unmapped R1 reads | remove bowtie information | awk trick picard | sam2fastq
  # the nasty awk command is because picard is so fucking pedantic
  cmds['picard.samtofastq.potential.chimeric.first'] <-
    paste0("samtools view -f 68 ", aln.vir.first.bam, " | cut -f 1-11 | ",
           "awk -F'\\t' 'BEGIN {OFS = FS} {$2=4; $3=$7=\"*\"; $4=$5=$8=0; print}' | ",
           java.e, " -jar ", path.to.picard.jar, " SamToFastq I=/dev/stdin", 
           " F=", potential.chimeric.reads.first.fastq,
           " INCLUDE_NON_PF_READS=true")
  cmds['picard.samtofastq.potential.chimeric.second'] <-
    paste0("samtools view -f 133 ", aln.vir.second.bam, " | cut -f 1-11 | ",
           "awk -F'\\t' 'BEGIN {OFS = FS} {$2=4; $3=$7=\"*\"; $4=$5=$8=0; print}' | ",
           java.e, " -jar ", path.to.picard.jar, " SamToFastq I=/dev/stdin", 
           " F=", potential.chimeric.reads.second.fastq,
           " INCLUDE_NON_PF_READS=true")
  cmds['picard.samtofastq.potential.chimeric.both'] <-
    paste0("samtools view -f 5 ", aln.vir.both.bam, " | cut -f 1-11 | ",
           "awk -F'\\t' 'BEGIN {OFS = FS} {$2=4; $3=$7=\"*\"; $4=$5=$8=0; print}' | ",
           java.e, " -jar ", path.to.picard.jar, " SamToFastq I=/dev/stdin", 
           " F=", potential.chimeric.reads.both.fastq,
           " INCLUDE_NON_PF_READS=true")
  
  # run local alignment on combined human/viral genome
  cmds['aln.ref_vir.1'] <- paste(bt2.e, "-p", threads, 
                                 "-x", ref_vir.genome.idx, 
                                 "-U", potential.chimeric.reads.first.fastq,
                                 "--very-sensitive-local",
                                 "-k 2 -R 2", "|", 
                                 samtools.e, "view -F 4 -O SAM -h |",
                                 samtools.e, "sort -n -O SAM |", 
                                 filterChimericBam,"|", samtools.e, "sort -O BAM >", 
                                 local.first.bam)
  cmds['aln.ref_vir.2'] <- paste(bt2.e, "-p", threads, 
                                 "-x", ref_vir.genome.idx, 
                                 "-U", potential.chimeric.reads.second.fastq,
                                 "--very-sensitive-local",
                                 "-k 2 -R 2", "|", 
                                 samtools.e, "view -F 4 -O SAM -h |",
                                 samtools.e, "sort -n -O SAM |", 
                                 filterChimericBam,"|", samtools.e, "sort -O BAM >",
                                 local.second.bam)
  cmds['aln.ref_vir.3'] <- paste(bt2.e, "-p", threads, 
                                 "-x", ref_vir.genome.idx, 
                                 "-U", potential.chimeric.reads.both.fastq,
                                 "--very-sensitive-local",
                                 "-k 2 -R 2", "|", 
                                 samtools.e, "view -F 4 -O SAM -h |",
                                 samtools.e, "sort -n -O SAM |", 
                                 filterChimericBam,"|", samtools.e, "sort -O BAM >",
                                 local.both.bam)
  if (cleanup) {
    cmds['mktemp'] <- "mkdir temp"
    cmds['mvtemp'] <- "mv -v *localalign.bam *dual-mapped.bam temp/"
    cmds['remove_bullshit'] <- "rm -v *bam* *fastq*"
    cmds['mvback'] <- "mv -v temp/* ."
    cmds['rmtemp'] <- "rmdir temp/"
  }
  
  
  # index everything
  # cmds['index'] <- paste0("for i in ", output.dir, "/*bam; do ", samtools.e, " index $i; done")
  # cmds['idxstats'] <- paste0("for i in ", output.dir, "/*bam; do ", samtools.e, " idxstats $i > $i.idxstats; done")
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
  c <- c(c, "echo Done.")
  return(c)
}










