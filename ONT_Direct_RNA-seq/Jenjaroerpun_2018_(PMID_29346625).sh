# Jenjaroerpun et al., 2018 (PMID 29346625)
# "Complete genomic and transcriptional landscape analysis using third-generation sequencing: a case study of Saccharomyces cerevisiae CEN.PK113-7D"
# Accession for Nanopore data: SRP116559
# Yeast strain: CEN.PK113-7D
# Lib prep protocol: Direct RNA sequencing (ONT SQK-RNA001 kit);

# Download FASTQ files:
echo -e "SRR5989374\tGlu_rep1
SRR6352887\tGlu_rep2
SRR6352888\tGlu_rep3
SRR6352889\tGlu_rep4
SRR5989373\tEth_rep1
SRR6352890\tEth_rep2
SRR6352891\tEth_rep3
SRR6352892\tEth_rep4" > Jenjar2018_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Jenjar2018_"$2".fastq.gz "; \
  system(cmd)}' Jenjar2018_acc.txt


# Count raw reads:
for file in Jenjar*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Download CENPK113-7D genome from Supplementary data to Boerlin et al., 2019 (PMID_30590648);

# Align to CEN.PK113-7D genome with Minimap2:
for file in Jenjar2018*fastq.gz; do 
  echo $file && 
  minimap2 -t 8 -ax splice -k 14 -L --cs --secondary=no \
    -G 10000 CENPK113-7D_Sequence.fasta $file > ${file/fastq.gz/sam}; 
done

# Convert SAM to sorted BAM:
for file in Jenjar2018*sam; do 
  echo $file && 
  samtools view -hu $file | 
  samtools sort - -o ${file/sam/bam}; 
done

# Filter out unmapped reads and reads with low MAPQ:
for file in Jenjar2018*bam; do 
  echo $file && 
  samtools view -hb -q 10 -F 4 $file > ${file/.bam/_mapq.bam}; 
done

# Count aligned reads:
for file in Jenjar2018*mapq.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | 
    awk '{print $1}'); 
done

# Post-process filtered BAM files in R session using Filter_Direct_RNA-seq_BAM_files.R;
