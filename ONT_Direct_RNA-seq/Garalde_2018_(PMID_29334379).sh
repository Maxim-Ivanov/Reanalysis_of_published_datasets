# Garalde et al. 2018 (PMID 29334379)
# "Highly parallel direct RNA sequencing on an array of nanopores"
# Accession PRJNA408327
# Yeast strain S288C
# Lib prep protocol: Direct RNA sequencing

# Download FASTQ files:
echo -e "SRR6059707\trep1
SRR6059706\trep2
SRR6059712\trep3
SRR6059711\trep4
SRR6059710\trep5" > Garalde2018_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Garalde2018_"$2".fastq.gz "; \
  system(cmd)}' Garalde2018_acc.txt

# Count raw reads:
for file in Garalde2018*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Align to SacCer3 genome with Minimap2:
for file in Garalde2018*fastq.gz; do 
  echo $file && 
  minimap2 -t 8 -ax splice -k 14 -L --cs --secondary=no \
    -G 10000 sacCer3.fa $file > ${file/fastq.gz/sam}; 
done

# Convert SAM to sorted BAM:
for file in Garalde2018*sam; do 
  echo $file && 
  samtools view -hu $file | 
  samtools sort - -o ${file/sam/bam}; 
done

# Filter out unmapped reads and reads with low MAPQ:
for file in Garalde2018*bam; do 
  echo $file && 
  samtools view -hb -q 10 -F 4 $file > ${file/.bam/_mapq.bam}; 
done

# Count aligned reads:
for file in Garalde2018*mapq.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | 
    awk '{print $1}'); 
done

# Post-process filtered BAM files in R session using Filter_Direct_RNA-seq_BAM_files.R;
