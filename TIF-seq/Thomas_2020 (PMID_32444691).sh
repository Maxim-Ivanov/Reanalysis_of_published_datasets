# Thomas et al., 2020 (PMID 32444691)
# "Transcript isoform sequencing reveals widespread promoter-proximal transcriptional termination in Arabidopsis"
# PE sequencing (Illumina NextSeq)
# UMI (8 nt) on R2
# R1 can start with poly(T) stretches of variable length;
# Due to unknown reason, these poly(T) stretches can start with a small offset (1-5 non-T bases) from 5' end of R1;
# Also R1 reads can contain internal poly(T) stretches;


# Download FASTQ files:
echo -e "SRR8870029\twt_rep1
SRR8870030\twt_rep2
SRR8870031\then2_rep1
SRR8870032\then2_rep2" > Thomas2020_acc.txt

awk '{cmd="fastq-dump --gzip --split-files "$1" && \
  mv "$1"_1.fastq.gz "$2"_R1.fastq.gz && \
  mv "$1"_2.fastq.gz "$2"_R2.fastq.gz"; \
  system(cmd)}' Thomas2020_acc.txt

# Count raw reads:
for file in Thomas2020*R1.fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | 
    awk '{print $1}') / 4 )); 
done

# Trim adapters from 3' ends of R1 and R2, trim UMI from 5' end of R2, trim polyT stretches from 5' ends of R1 (up to 4 non-T + at least 8 T), discard R1 reads with internal polyT (at least 12 T):
 
for f1 in Thomas2020*R1.fastq.gz; do 
  f2=${f1/R1/R2} && 
  name=${f1/R1.fastq.gz/} && 
  echo $f1 "+" $f2 && 
  cutadapt -j 4 -a "AGGTGACCGG" -A "AGGTGACCGG" \
    -a "AGATCGGAAG" -A "AGATCGGAAG" --nextseq-trim=20 \
    --match-read-wildcards --minimum-length 28 -o cut1 \
    -p cut2 <(zcat $f1) <(zcat $f2) && 
  umi_tools extract -I cut2 --extract-method=string \
    --bc-pattern=NNNNNNNN --read2-in=cut1 --stdout=umi2 \
    --read2-out=umi1 && 
  rm cut1 cut2 && 
  awk -v name=$name 'BEGIN{OFS="\n"}\
    {h=$0; getline seq; getline qh; getline qseq; \
    if (NR==FNR) {match(seq, /^[^T]{0,4}T{8,}/); \
      end=RSTART+RLENGTH; \
      if (RSTART!=0) \
        {if (length(seq)-end > 20) \
          {seq=substr(seq, end, length(seq)); \
          qseq=substr(qseq, end, length(qseq))} \
        else {arr[FNR]=1; next}}; \
      if (seq~/T{12,}/) {arr[FNR]=1; next}; \
      print h, seq, qh, qseq > name"R1_trimmed.fq"} \
    else {if (arr[FNR]==1) next; \
      print h, seq, qh, qseq > name"R2_trimmed.fq"}}' \
  umi1 umi2 && 
  rm umi1 umi2; 
done

# Count trimmed reads:
for file in Thomas2020*R1.trimmed.fq; do 
  echo $file $(( $(wc -l $file | \
    awk '{print $1}') / 4 )); 
done

# Gzip trimmed files:
for file in Thomas2020*trimmed.fq; do 
  echo $file && gzip $file; 
done

# Align to TAIR10:
for f1 in Thomas2020*R1_trimmed.fq.gz; do 
  f2=${f1/R1/R2} && 
  echo $f1 $f2 && 
  STAR --genomeDir tair10_star --readFilesIn $f1 $f2 \
    --runThreadN 4 --outSAMmultNmax 1 --alignEndsType EndToEnd \
    --outFileNamePrefix ${f1/R1_trimmed.fq.gz/} \
    --readFilesCommand zcat --outSAMtype BAM Unsorted \
    --alignIntronMax 10000 --alignMatesGapMax 25000; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Count properly paired reads:
for file in Thomas*bam; do 
  echo $file $(( $(samtools flagstat $file | sed -n '9p' | \
    awk '{print $1}') / 2 )); 
done

# Remove low MAPQ reads and reads not in proper pair. Sort and index BAM files:
for file in Thomas2020*bam; do 
  echo $file && 
  samtools view -hq 10 -f 2 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Deduplicate by UMIs:
for file in Thomas2020*mapq.bam; do 
  echo $file && 
  samtools index $file && 
  umi_tools dedup -I $file -S ${file/.bam/_dedup.bam} \
    --method cluster --paired; 
done

