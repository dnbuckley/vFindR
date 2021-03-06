printf "\n######## COMMAND: ########\n"
echo "mkdir /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output"
printf "##########################\n"
mkdir /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output
printf "\n######## COMMAND: ########\n"
echo "cd /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output"
printf "##########################\n"
cd /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output
printf "\n######## COMMAND: ########\n"
echo "/Users/dbuckley/Desktop/bin/bowtie2 -p 8 -x ~/Desktop/salhia_lab/genomes/GRCh38_21/GRCh38_21.fa -1 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/fake_R1.fastq.gz -2 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/fake_R2.fastq.gz  | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam"
printf "##########################\n"
/Users/dbuckley/Desktop/bin/bowtie2 -p 8 -x ~/Desktop/salhia_lab/genomes/GRCh38_21/GRCh38_21.fa -1 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/fake_R1.fastq.gz -2 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/fake_R2.fastq.gz  | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam
printf "\n######## COMMAND: ########\n"
echo "{ /usr/local/bin/samtools view -H  /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 69 -G 9 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 137 -G 5 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam; } | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.first.bam"
printf "##########################\n"
{ /usr/local/bin/samtools view -H  /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 69 -G 9 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 137 -G 5 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam; } | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.first.bam
printf "\n######## COMMAND: ########\n"
echo "{ /usr/local/bin/samtools view -H  /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 133 -G 9 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 73 -G 5 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam; } | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.second.bam"
printf "##########################\n"
{ /usr/local/bin/samtools view -H  /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 133 -G 9 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 73 -G 5 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam; } | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.second.bam
printf "\n######## COMMAND: ########\n"
echo "{ /usr/local/bin/samtools view -H  /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 77 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 141 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam; } | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.both.bam"
printf "##########################\n"
{ /usr/local/bin/samtools view -H  /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 77 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam && /usr/local/bin/samtools view -f 141 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam; } | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.both.bam
printf "\n######## COMMAND: ########\n"
echo "/usr/bin/java -jar ~/Desktop/bin/picard.jar SamToFastq I=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.first.bam F=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.first_R1.fastq.gz F2=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.first_R2.fastq.gz"
printf "##########################\n"
/usr/bin/java -jar ~/Desktop/bin/picard.jar SamToFastq I=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.first.bam F=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.first_R1.fastq.gz F2=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.first_R2.fastq.gz
printf "\n######## COMMAND: ########\n"
echo "/usr/bin/java -jar ~/Desktop/bin/picard.jar SamToFastq I=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.second.bam F=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.second_R1.fastq.gz F2=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.second_R2.fastq.gz"
printf "##########################\n"
/usr/bin/java -jar ~/Desktop/bin/picard.jar SamToFastq I=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.second.bam F=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.second_R1.fastq.gz F2=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.second_R2.fastq.gz
printf "\n######## COMMAND: ########\n"
echo "/usr/bin/java -jar ~/Desktop/bin/picard.jar SamToFastq I=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.both.bam F=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.both_R1.fastq.gz F2=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.both_R2.fastq.gz"
printf "##########################\n"
/usr/bin/java -jar ~/Desktop/bin/picard.jar SamToFastq I=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.umapped.both.bam F=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.both_R1.fastq.gz F2=/Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.both_R2.fastq.gz
printf "\n######## COMMAND: ########\n"
echo "/Users/dbuckley/Desktop/bin/bowtie2 --very-sensitive -p 8 -x ~/Desktop/salhia_lab/genomes/refseq_all_virus/refseq_all_virus.fa -1 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.first_R1.fastq.gz -2 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.first_R2.fastq.gz  | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_unmapped-first-human_virus.bam"
printf "##########################\n"
/Users/dbuckley/Desktop/bin/bowtie2 --very-sensitive -p 8 -x ~/Desktop/salhia_lab/genomes/refseq_all_virus/refseq_all_virus.fa -1 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.first_R1.fastq.gz -2 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.first_R2.fastq.gz  | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_unmapped-first-human_virus.bam
printf "\n######## COMMAND: ########\n"
echo "/Users/dbuckley/Desktop/bin/bowtie2 --very-sensitive -p 8 -x ~/Desktop/salhia_lab/genomes/refseq_all_virus/refseq_all_virus.fa -1 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.second_R1.fastq.gz -2 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.second_R2.fastq.gz  | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_unmapped-second-human_virus.bam"
printf "##########################\n"
/Users/dbuckley/Desktop/bin/bowtie2 --very-sensitive -p 8 -x ~/Desktop/salhia_lab/genomes/refseq_all_virus/refseq_all_virus.fa -1 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.second_R1.fastq.gz -2 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.second_R2.fastq.gz  | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_unmapped-second-human_virus.bam
printf "\n######## COMMAND: ########\n"
echo "/Users/dbuckley/Desktop/bin/bowtie2 --very-sensitive -p 8 -x ~/Desktop/salhia_lab/genomes/refseq_all_virus/refseq_all_virus.fa -1 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.both_R1.fastq.gz -2 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.both_R2.fastq.gz  | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_unmapped-both-human_virus.bam"
printf "##########################\n"
/Users/dbuckley/Desktop/bin/bowtie2 --very-sensitive -p 8 -x ~/Desktop/salhia_lab/genomes/refseq_all_virus/refseq_all_virus.fa -1 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.both_R1.fastq.gz -2 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.unmapped.both_R2.fastq.gz  | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_unmapped-both-human_virus.bam
printf "\n######## COMMAND: ########\n"
echo "/usr/local/bin/samtools index /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam"
printf "##########################\n"
/usr/local/bin/samtools index /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam
printf "\n######## COMMAND: ########\n"
echo "{ /usr/local/bin/samtools view -F 4 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_unmapped-first-human_virus.bam && /usr/local/bin/samtools view -F 4 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_unmapped-second-human_virus.bam; } | cut -f 1 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human_dual-mapped_readnames.txt"
printf "##########################\n"
{ /usr/local/bin/samtools view -F 4 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_unmapped-first-human_virus.bam && /usr/local/bin/samtools view -F 4 /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_unmapped-second-human_virus.bam; } | cut -f 1 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human_dual-mapped_readnames.txt
printf "\n######## COMMAND: ########\n"
echo "/Users/dbuckley/anaconda3/bin/python3 ~/bin/extract_reads.py -b /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam -n /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human_dual-mapped_readnames.txt -o /dev/stdout | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human_dual-mapped.bam"
printf "##########################\n"
/Users/dbuckley/anaconda3/bin/python3 ~/bin/extract_reads.py -b /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human.bam -n /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human_dual-mapped_readnames.txt -o /dev/stdout | samtools sort -O BAM -@ 8 > /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/vFindR_human_dual-mapped.bam
printf "\n######## COMMAND: ########\n"
echo "for i in /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/*bam; do /usr/local/bin/samtools index $i; done"
printf "##########################\n"
for i in /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/*bam; do /usr/local/bin/samtools index $i; done
printf "\n######## COMMAND: ########\n"
echo "for i in /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/*bam; do /usr/local/bin/samtools idxstats $i > $i.idxstats; done"
printf "##########################\n"
for i in /Users/dbuckley/Desktop/salhia_lab/vFindeR/fake_data/vFindR_output/*bam; do /usr/local/bin/samtools idxstats $i > $i.idxstats; done
