setwd("~/Desktop/hpc_salhia_lab/dbuckley/vFindR_test/VIcaller_example_input/")
source("~/Desktop/salhia_lab/vFindR/src/vFindR_v1.1.R")

vFindR(sample.dir = "rundir/",
       ref.genome.idx = "/project/salhia_618/dbuckley/genomes/hg38/hg38.fa",
       vir.genome.idx = "/project/salhia_618/dbuckley/genomes/refseq_all_virus/refseq_all_virus.fa", 
       path.to.picard.jar = "/project/salhia_618/bin/picard.jar", 
       bt2.e = "/project/salhia_618/bin/bowtie2", 
       path.to.extract.py = "/project/salhia_618/scripts/other_NGS_tools/vFindR/src/extract_reads.py", 
       python.e = "~/.conda/envs/def/bin/python",
       samtools.e = "/project/salhia_618/bin/samtools", 
       bamtools.e = "/spack/apps/linux-centos7-x86_64/intel-19.0.4/bamtools-2.5.1-sve2uhoew43gadknqjuaedcvwx2scgc2/bin/bamtools",
       java.e = "/project/salhia_618/bin/java", 
       threads = 20, 
       vFindR.dir = "/project/salhia_618/scripts/other_NGS_tools/vFindR",
       mode = "sh")


vFindR(sample.dir = "rundir/",
       ref.genome.idx = "/project/salhia_618/dbuckley/genomes/hg38/hg38.fa",
       vir.genome.idx = "/project/salhia_618/dbuckley/genomes/refseq_all_virus/refseq_all_virus.fa", 
       path.to.picard.jar = "/project/salhia_618/bin/picard.jar", 
       bt2.e = "/project/salhia_618/bin/bowtie2", 
       path.to.extract.py = "/project/salhia_618/scripts/other_NGS_tools/vFindR/src/extract_reads.py", 
       python.e = "~/.conda/envs/def/bin/python",
       samtools.e = "/project/salhia_618/bin/samtools", 
       bamtools.e = "/spack/apps/linux-centos7-x86_64/intel-19.0.4/bamtools-2.5.1-sve2uhoew43gadknqjuaedcvwx2scgc2/bin/bamtools",
       java.e = "/project/salhia_618/bin/java", 
       threads = 20, 
       vFindR.dir = "/project/salhia_618/scripts/other_NGS_tools/vFindR",
       mode = "slurm")

dirs <- list.dirs(".", recursive = F)
dirs <- gsub("^\\.\\/", "", dirs)
for (dir in dirs) {
  vFindR(sample.dir = dir,
         output.name = paste0("vFindR_", dir),
         ref.genome.idx = "/project/salhia_618/dbuckley/genomes/hg38/hg38.fa",
         vir.genome.idx = "/project/salhia_618/dbuckley/genomes/refseq_all_virus/refseq_all_virus.fa", 
         ref_vir.genome.idx = "/project/salhia_618/dbuckley/genomes/hg38-refseq_all_virus/hg38-refseq_all_virus.fa",
         path.to.picard.jar = "/project/salhia_618/bin/picard.jar", 
         bt2.e = "/project/salhia_618/bin/bowtie2", 
         path.to.extract.py = "/project/salhia_618/scripts/other_NGS_tools/vFindR/src/extract_reads.py", 
         python.e = "~/.conda/envs/def/bin/python",
         samtools.e = "/project/salhia_618/bin/samtools", 
         bamtools.e = "/spack/apps/linux-centos7-x86_64/intel-19.0.4/bamtools-2.5.1-sve2uhoew43gadknqjuaedcvwx2scgc2/bin/bamtools",
         java.e = "/project/salhia_618/bin/java", 
         threads = 20, 
         vFindR.dir = "/project/salhia_618/scripts/other_NGS_tools/vFindR",
         mode = "slurm",
         slurm.header.args = c("-t 2-00:00:00", "-c 16", "--mem=55G"))
}


dirs <- list.dirs("rundir/", recursive = F, full.names = T)
dirs <- gsub("^\\.\\/", "", dirs)
for (dir in dirs) {
  vFindR(sample.dir = dir,
         output.name = paste0("vFindR_", dir),
         ref.genome.idx = "/project/salhia_618/dbuckley/genomes/hg38/hg38.fa",
         vir.genome.idx = "/project/salhia_618/dbuckley/genomes/refseq_all_virus/refseq_all_virus.fa", 
         ref_vir.genome.idx = "/project/salhia_618/dbuckley/genomes/hg38-refseq_all_virus/hg38-refseq_all_virus.fa",
         path.to.picard.jar = "/project/salhia_618/bin/picard.jar", 
         bt2.e = "/project/salhia_618/bin/bowtie2", 
         path.to.extract.py = "/project/salhia_618/scripts/other_NGS_tools/vFindR/src/extract_reads.py", 
         python.e = "~/.conda/envs/def/bin/python",
         samtools.e = "/project/salhia_618/bin/samtools", 
         bamtools.e = "/spack/apps/linux-centos7-x86_64/intel-19.0.4/bamtools-2.5.1-sve2uhoew43gadknqjuaedcvwx2scgc2/bin/bamtools",
         java.e = "/project/salhia_618/bin/java", 
         threads = 20, 
         vFindR.dir = "/project/salhia_618/scripts/other_NGS_tools/vFindR",
         mode = "slurm")
}










