# Dutta et al., 2017 (PMID 28249159)
# "Composition and Function of Mutant Swi/Snf Complexes"
# Accession GSE81722;
# PE sequencing;

# Victorino et al., 2020 (PMID 32187185)
# "RNA Polymerase II CTD phosphatase Rtr1 fine-tunes transcription termination"
# Accession GSE135056;
# PE sequencing;


# Download WT samples from Dutta 2017:
echo -e "SRR3568018\trep1
SRR3568019\trep2" > Dutta2017_acc.txt

awk '{cmd="fastq-dump --gzip --split-files "$1" && \
  mv "$1"_1.fastq.gz RNAseq_Dutta2017_"$2"_R1.fq.gz && \
  mv "$1"_2.fastq.gz RNAseq_Dutta2017_"$2"_R2.fq.gz"; \
  system(cmd)}' Dutta2017_acc.txt

# Download WT and rrp6 samples from Victorino 2020:

echo -e "SRR9856226\tWT_rep1_part1
SRR9856227\tWT_rep1_part2
SRR9856228\tWT_rep1_part3
SRR9856229\tWT_rep2_part1
SRR9856230\tWT_rep2_part2
SRR9856231\tWT_rep2_part3
SRR9856232\tWT_rep3_part1
SRR9856233\tWT_rep3_part2
SRR9856234\tWT_rep3_part3
SRR9856235\tWT_rep4_part1
SRR9856236\tWT_rep4_part2
SRR9856237\tWT_rep4_part3
SRR9856238\tRRP6_rep1_part1
SRR9856239\tRRP6_rep1_part2
SRR9856240\tRRP6_rep1_part3
SRR9856241\tRRP6_rep2_part1
SRR9856242\tRRP6_rep2_part2
SRR9856243\tRRP6_rep2_part3
SRR9856244\tRRP6_rep3_part1
SRR9856245\tRRP6_rep3_part2
SRR9856246\tRRP6_rep3_part3
SRR9856247\tRRP6_rep4_part1
SRR9856248\tRRP6_rep4_part2
SRR9856249\tRRP6_rep4_part3" > Victorino2020_acc.txt


awk '{cmd="fastq-dump --gzip --split-files "$1" && \
  mv "$1"_1.fastq.gz RNAseq_Victorino2020_"$2"_R1.fq.gz && \
  mv "$1"_2.fastq.gz RNAseq_Victorino2020_"$2"_R2.fq.gz"; \
  system(cmd)}' Victorino2020_acc.txt

# Victorino 2020 data: merge all runs of the same replicate:
for f1 in RNAseq_Victorino2020*part1_R?.fq.gz; do 
  f2=${f1/part1/part2} && 
  f3=${f1/part1/part3} && 
  sample=${f1/_part1/} && 
  echo $sample && zcat $f1 $f2 $f3 | gzip > $sample; 
done

rm RNAseq_Victorino2020*part*fq.gz

# Count raw read pairs:
for file in RNAseq*R1.fq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done


# Align to the SacCer3 genome:
for f1 in RNAseq*R1.fq.gz; do 
  f2=${f1/_R1/_R2} && 
  echo $f1 $f2 && 
  STAR --genomeDir saccer3_star --readFilesIn $f1 $f2 \
    --runThreadN 4 --outFileNamePrefix ${file/R1.fq.gz/} \
    --clip3pAdapterSeq AGATCGGAAGAGC --outSAMmultNmax 1 \
    --alignEndsType Local --readFilesCommand zcat \
    --outSAMtype BAM Unsorted; 
done

rm -r *out *tab *STARtmp
for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done

# Count aligned read pairs:
for file in RNAseq*bam; do 
  echo $file $(samtools flagstat $file | \
    sed -n '9p' | awk '{sprintf("%i", $1/2)}'); 
done

# Filter for reads with high MAPQ in proper pairs:
for file in RNAseq*bam; do 
  echo $file && 
  samtools view -hb -f 2 -q 10 $file -o ${file/.bam/_mapq.bam}; 
done

# Run samtools fixmate (required by samtools markdup):
for file in RNAseq*mapq.bam; do 
  echo $file && 
  samtools fixmate -m $file ${file/.bam/_fixmate.bam}; 
done

# Sort by coordinates:
for file in RNAseq*fixmate.bam; do 
  echo $file && 
  samtools sort $file -o ${file/.bam/_sorted.bam}; 
done

# Deduplicate on start coordinates of paired reads:
for file in RNAseq*sorted.bam; do 
  echo $file && 
  samtools markdup -r $file ${file/.bam/_dedup.bam}; 
done

# Load deduplicated BAM files into R session and produce stranded normalized Bedgraph files using convert_RNAseq_BAM_to_normalized_bedGraph.R;



