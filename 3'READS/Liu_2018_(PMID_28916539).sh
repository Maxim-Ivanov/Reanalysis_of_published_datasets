# Liu et al. 2018 (PMID 28916539);
# "Comparative analysis of alternative polyadenylation in S. cerevisiae and S. pombe";
# S.cerevisiae, 3'READS lib prep protocol, single-end reads;
# 5' end: UMI (4 nt) + few non-template T's;
# 3' end: Illumina Small RNA adapter;

# Download FASTQ files from SRA:
echo -e "SRR5276076\tLiu2018_s01_wt_MM_R1
SRR5276078\tLiu2018_s02_wt_MM_R2
SRR5276077\tLiu2018_s03_wt_RM_R1
SRR5276079\tLiu2018_s04_wt_RM_R2" > Liu2018_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz "$2".fastq.gz"; system(cmd)}' Liu2018_acc.txt

# Count raw reads:
for file in Liu2018*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Trim 3' Illumina Small RNA adapters:
for file in Liu2018*fastq.gz; do 
  echo $file && 
  trim_galore --small_rna --three_prime_clip_R1 4 --gzip $file; 
done

# Process UMI on 5' end (SE mode):
for file in Liu2018*trimmed.fq.gz; do 
  echo $file && 
  umi_tools extract --stdin=${file} --bc-pattern=NNNN \
    --stdout=${file/.fq.gz/_UMI.fq.gz}; 
done

# Align to SacCer3 in Local mode:
for file in Liu2018*fq.gz; do 
  echo $file && 
  STAR --genomeDir saccer3_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Local \
    --readFilesCommand zcat --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Count aligned reads:
for file in Liu2018*bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' \
    | awk '{print $1}'); 
done

# Sort BAM files and remove low MAPQ reads:
for file in Liu2018*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Deduplicate on UMI:
for file in Liu2018*mapq.bam; do 
  echo $file && 
  samtools index $file && 
  umi_tools dedup --stdin=${file} \
    --stdout=${file/.bam/_dedup.bam}; 
done

# Postprocess BAM files in R: see Liu_2018_(PMID_28916539).R


