# Zhu et al. 2018 (PMID 30374093);
# "RNA polymerase II activity revealed by GRO-seq and pNET-seq in Arabidopsis"
# Accession GSE109974;
# Bioo Scientific Small RNA-seq kit v3;
# UMI (4 nt);
# PE sequencing;
# First non-UMI base of R2 = position of RNAPII;


# Download FASTQ files for WT samples:
echo -e "SRR6661081\ttotal_rep1
SRR6661082\ttotal_rep2
SRR6661083\tunph_rep1
SRR6661084\tunph_rep2
SRR6661085\tSer2P_rep1
SRR6661086\tSer2P_rep2
SRR6661087\tSer5P_rep1
SRR6661088\tSer5P_rep2" > pNET_Zhu2018_acc.txt

awk '{cmd="fastq-dump --gzip --split-files "$1" && \
  mv "$1"_1.fastq.gz pNET_Zhu2018_"$2"_R1.fastq.gz && \
  mv "$1"_2.fastq.gz pNET_Zhu2018_"$2"_R2.fastq.gz"; \
  system(cmd)}' pNET_Zhu2018_acc.txt

# Count raw reads:
for file in pNET_Zhu2018*R2.fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Process UMIs (PE mode):
for f1 in pNET_Zhu2018*R1.fastq.gz; do 
  f2=${f1/_R1/_R2} && 
  echo $f1 $f2 && 
  umi_tools extract --stdin=${f1} --read2-in=${f2} \
    --bc-pattern=NNNN --bc-pattern2=NNNN \
    --stdout=${f1/.fastq.gz/_UMI.fq.gz} \
    --read2-out=${f2/.fastq.gz/_UMI.fq.gz}; 
done

# Align to TAIR10 (SE mode, only R2):
for f2 in pNET_Zhu2018*R2_UMI.fq.gz; do 
  echo $f2 && 
  STAR --genomeDir tair10_star --readFilesIn $f2 \
    --runThreadN 4 --outFileNamePrefix ${f2/R2_UMI.fq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq GATCGTCGGACT \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARgenome *STARtmp *out *tab

# Sort BAM files:
for file in pNET_Zhu2018*bam; do 
  echo $file && 
  samtools sort $file -o ${file/.bam/_sorted.bam}; 
done

# Count aligned reads:
for file in pNET_Zhu2018*sorted.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Deduplicate (UMI-Tools):
for file in pNET_Zhu2018*sorted.bam; do 
  echo $file && 
  samtools index $file && 
  umi_tools dedup --stdin=${file} \
    --stdout=${file/.bam/_dedup.bam}; 
done

# Download Araport11 annotation:
ann="Araport11_GFF3_genes_transposons.201606.gff.gz"
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/${ann}

# Extract rRNA, tRNA, snRNA and snoRNA genes:
zcat $ann | sed '/^#/d' | 
awk 'BEGIN{OFS="\t"}\
  {if ($3=="rRNA" || $3=="tRNA" || $3=="snRNA" || $3=="snoRNA")\
  print $1, $4-100, $5+100, $3, ".", $7}' | 
sed 's/^ChrC/Pt/;s/^ChrM/Mt/;s/^Chr//' > \
Araport11_short_ncRNA_ext100bp.bed

# Add the whole Mt and Pt chromosomes to the BED file:
echo -e "Pt\t1\t154478\tPt\t.\t+
Pt\t1\t154478\tPt\t.\t-
Mt\t1\t366924\tMt\t.\t+
Mt\t1\t366924\tMt\t.\t-" >> Araport11_short_ncRNA_ext100bp.bed

# Remove intersections with non-RNAPII transcribed genes:
for file in pNET_Zhu2018*dedup.bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b Araport11_short_ncRNA_ext100bp.bed > \
    ${file/.bam/_clean.bam}; 
done

# Remove low MAPQ reads:
for file in pNET_Zhu2018*clean.bam; do 
  echo $file && 
  samtools view -hb -q 10 $file > ${file/.bam/_mapq.bam}; 
done

# Load filtered BAM files into R session (pNET-seq_pipeline.R).
