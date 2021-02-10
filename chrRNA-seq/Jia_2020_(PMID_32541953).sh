# Jia et al., 2020 (PMID 32541953)
# "Post-transcriptional splicing of nascent RNA contributes to widespread intron retention in plants"
# Isolation of chromosome-bound RNA followed by full-length cDNA sequencing on the Oxford Nanopore (ONT) platform;
# The chrosome-bound RNA is expected to be enriched for nascent RNA molecules;
# The cDNA synthesis protocol was based on Takara (Clontech) SMARTer kit, however the SMART CDS Primer II A was replaced by NEB Universal miRNA Cloning Linker. This modification allows to recover strand orientation of the original RNA molecule (see Clontech_SMARTer_oligo_system.pdf);

# Download FASTQ files:

acc="SRR10538401"; fastq-dump --gzip $acc && \
  mv ${acc}.fastq.gz Jia2020_rep1.fastq.gz
acc="SRR10538409"; fastq-dump --gzip $acc && \
  mv ${acc}.fastq.gz Jia2020_rep1.fastq.gz

# Count raw reads:
for file in Jia2020*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Align to TAIR10 genome with Minimap2:
for file in Jia2020*fastq.gz; do 
  echo $file && 
  minimap2 -t 8 -ax splice -k 14 -L --cs --secondary=no \
    -G 10000 TAIR10.fa $file > ${file/fastq.gz/sam}; 
done

# Convert SAM to sorted BAM:
for file in Jia2020*sam; do 
  echo $file && 
  samtools view -hu $file | 
  samtools sort - -o ${file/sam/bam}; 
done

# Filter out unmapped reads and reads with low MAPQ:
for file in Jia2020*bam; do 
  echo $file && 
  samtools view -hb -q 10 -F 4 $file > ${file/.bam/_mapq.bam}; 
done

# Count filtered reads:
for file in Jia2020*mapq.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# To recover strand info of the reads, continue in R session (Jia_2020_(PMID_32541953).R);