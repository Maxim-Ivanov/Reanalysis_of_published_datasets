# Inagaki et al., 2021 (PMID 33649596)
# Illumina 2x150 sequencing
# KAPA Stranded mRNA-seq kit


# Download the data:
echo -e "DRR235371\tfld4_rep2
DRR235370\tWT_rep2
DRR235369\tfld4_rep1
DRR235368\tWT_rep1" | \
awk '{cmd="fastq-dump --gzip --split-files "$1" && \
mv "$1"_1.fastq.gz Inagaki2021_"$2"_R1.fq.gz && \
mv "$1"_2.fastq.gz Inagaki2021_"$2"_R2.fq.gz"; system(cmd)}'

# Count raw reads:
for file in Inagaki2021*R1.fq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Align to TAIR10 genome:
for f1 Inagaki2021*R1.fq.gz; do 
  f2=${f1/_R1/_R2} && 
  echo $f1 $f2 &&
  STAR --genomeDir tair10_star --readFilesIn $f1 $f2 \
    --readFilesCommand zcat --runThreadN 4 \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix ${f1/R1.fq.gz/} --outSAMmultNmax 1 \
    --alignEndsType Local --clip3pAdapterSeq AGATCGGAAGAGC;
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Count aligned read pairs:
for file in Inagaki2021*bam; do 
  echo $file $(samtools flagstat $file | \
    sed -n '9p' | awk '{printf("%i", $1/2)}'); 
done

# Filter for reads with high MAPQ in proper pairs:
for file in Inagaki2021*bam; do 
  echo $file && 
  samtools view -hb -f 2 -q 10 $file -o ${file/.bam/_mapq.bam}; 
done

# Run samtools fixmate (required by samtools markdup):
for file in Inagaki2021*mapq.bam; do 
  echo $file && 
  samtools fixmate -m $file ${file/.bam/_fixmate.bam}; 
done

# Sort by coordinates:
for file in Inagaki2021*fixmate.bam; do 
  echo $file && 
  samtools sort $file -o ${file/.bam/_sorted.bam}; 
done

# Deduplicate on start coordinates of paired reads:
for file in Inagaki2021*sorted.bam; do 
  echo $file && 
  samtools markdup -r $file ${file/.bam/_dedup.bam}; 
done

# Remove rRNA, tRNA and snoRNA reads:
# (for details on Arabidopsis short ncRNA annotation, see https://github.com/Maxim-Ivanov/Reanalysis_of_published_datasets/blob/main/pNET-seq/pNET-seq_pipeline.sh)

for file in Inagaki2021*dedup.bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b Araport11_short_ncRNA_ext100bp.bed \
    > ${file/.bam/_clean.bam}; 
done

# Load deduplicated BAM files into R session and produce stranded normalized Bedgraph files (see https://github.com/Maxim-Ivanov/Reanalysis_of_published_datasets/blob/main/RNA-seq/convert_RNAseq_BAM_to_normalized_bedGraph.R)