# Mayer et al. 2015 (PMID 25910208)
# "Native elongating transcript sequencing reveals human transcriptional activity at nucleotide resolution"
# Accession: GSE61332;
# SE sequencing;
# Adapter: ATCTCGTATGCCG;
# UMIs (6 nt);

echo -e "SRR1575919\tMayer2015_HeLa_rep1_p1
SRR1575920\tMayer2015_HeLa_rep1_p2
SRR1575914\tMayer2015_HeLa_rep2_p1
SRR1575915\tMayer2015_HeLa_rep2_p2
SRR1575916\tMayer2015_HeLa_rep2_p3
SRR1575917\tMayer2015_HeLa_rep2_p4
SRR1575918\tMayer2015_HeLa_rep2_p5
SRR1575921\tMayer2015_HEK293T_rep1_p1
SRR1575922\tMayer2015_HEK293T_rep1_p2
SRR1575923\tMayer2015_HEK293T_rep2_p1
SRR1575924\tMayer2015_HEK293T_rep2_p2
SRR1575925\tMayer2015_HEK293T_rep2_p3
SRR1575926\tMayer2015_HEK293T_rep2_p4
SRR1575927\tMayer2015_HEK293T_rep2_p5
SRR1575928\tMayer2015_HEK293T_rep2_p6
SRR1928223\tMayer2015_HEK293T_rep2_p7" > Mayer2015_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv ${acc}.fastq.gz "$2".fastq.gz"; \
  system(cmd)}' Mayer2015_acc.txt

for file in Mayer2015*p1.fastq.gz; do 
  echo ${file/p1/p?} && zcat ${file/p1/p?} | 
  gzip > ${file/_p1/} && rm ${file/p1/p?}; 
done


# Nojima et al., 2015 (PMID 25910207)
# "Mammalian NET-Seq Reveals Genome-wide Nascent Transcription Coupled to RNA Processing"
# Illumina TruSeq Small RNA kit, no UMIs;
# Paired-end reads. First base of R2 corresponds to the RNAPII position;
# Adapter sequence = GATCGTCGGACT (3' adapter for R2 in TruSeq Small RNA kit);

# Nojima et al., 2018 (PMID 30340024)
# "RNA Polymerase II Phosphorylated on CTD Serine 5 Interacts with the Spliceosome during Co-transcriptional Splicing"
# Human (HeLa) or mouse (TAP) cells;
# The same lib prep and sequencing protocol as in Nojima 2015;

echo -e "SRR1544602\tNojima2015_8WG16_rep1
SRR1544603\tNojima2015_8WG16_rep2
SRR1544606\tNojima2015_CMA601
SRR1544604\tNojima2015_CMA602
SRR1544605\tNojima2015_CMA603
SRR6290113\tNojima2018_HeLa_S5P_rep1
SRR6290115\tNojima2018_HeLa_S5P_rep2
SRR6290117\tNojima2018_HeLa_S2P_rep1
SRR6290119\tNojima2018_HeLa_S2P_rep2
SRR6290127\tNojima2018_TAP_S5P_rep1
SRR6290129\tNojima2018_TAP_S5P_rep2" > Nojima_acc.txt

awk '{cmd="fastq-dump --gzip --split-files "$1" && \
  mv "$1"_1.fastq.gz "$2"_R1.fastq.gz && \
  mv "$1"_.fastq.gz "$2"_R2.fastq.gz"; \
  system(cmd)}' Nojima_acc.txt

# Delete R1 files:
rm Nojima*R1.fastq.gz
for file in Nojima*R2.fastq.gz; do mv $file ${file/_R2/}; done


##### Prepare human genomic files --------------------------

# Generate STAR index for the human genome (hg38):
wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta

STAR --runMode genomeGenerate \
  --genomeFastaFiles Homo_sapiens_assembly38.fasta \
  --runThreadN 4 --genomeDir hg38_star

# Download hg38 annotation from Ensembl:
ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz

# Fix the chromosome names:
zcat Homo_sapiens.GRCh38.84.gtf.gz | 
sed '/^[^#]/s/^/chr/' > Homo_sapiens.GRCh38.84.chr.gtf

# Extract rRNA, snRNA and snoRNA genes from GFF3
cat Homo_sapiens.GRCh38.84.chr.gtf | 
awk 'BEGIN{OFS="\t"}{match($9, /Name=[^;][^;]*/); \
  if ($3=="rRNA" || $3=="snRNA" || $3=="snoRNA") \
  print "chr"$1,$4,$5,substr($9,RSTART+5,RLENGTH-5),$3,$7}' > \
  hg38_rRNA_snRNA_snoRNA.bed

# Make BED file with tRNA annotation from GtRNAdb:
http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/Hsapi38-gene-list.html

sed 's/[^\t.][^\t.]*\.//;s/(//;s/)//' hg38_tRNA_GtRNAdb.txt > hg38_tRNA_GtRNAdb.bed

# Add tRNA genes:
cat hg38_rRNA_snRNA_snoRNA.bed hg38_tRNA_GtRNAdb.bed | 
sort -k1,1 -k2,2n > hg38_rRNA_tRNA_snRNA_snoRNA.bed


##### Prepare mouse genomic files ---------------------------

# Download mm10 genome from UCSC:
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit

twoBitToFa mm10.2bit mm10.fa

# Generate STAR genome index from mm10:
STAR --runMode genomeGenerate --genomeFastaFiles mm10.fa \
  --runThreadN 4 --genomeDir mm10_star

# Find coordinates of rRNA, snRNA and snoRNA in mm10:
wget ftp://ftp.ensembl.org/pub/release-94/gff3/mus_musculus/Mus_musculus.GRCm38.94.chr.gff3.gz
zcat Mus_musculus.GRCm38.94.chr.gff3.gz | 
sed '/^#/d' | 
awk 'BEGIN{OFS="\t"}\
  {if ($3=="rRNA" || $3=="snRNA" || $3=="snoRNA") \
  print "chr"$1,$4,$5,$3,100,$7}' > mm10_rRNA_snRNA_snoRNA.bed

# Get coordinates of tRNA from GtRNAdb:
http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/Mmusc10-gene-list.html

cat mm10_tRNA_GtRNAdb.txt | 
sed 's/[^\t.][^\t.]*\.//;s/(//;s/)//' > mm10_tRNA_GtRNAdb.bed

# Combine all unwanted ncRNA:
cat mm10_tRNA_GtRNAdb.bed | 
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"tRNA",$5,$6}' | 
cat - mm10_rRNA_snRNA_snoRNA.bed | 
sort -k1,1 -k2,2n > mm10_rRNA_tRNA_snRNA_snoRNA.bed


##### Process Mayer 2015 data (human mNET-seq) ----------------

# Count input reads:
for file in Mayer2015*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Trim UMIs:
for file in Mayer2015*fastq.gz; do 
  echo $file && 
  umi_tools extract --stdin=${file} --bc-pattern=NNNNNN \
    --stdout=${file/.fastq.gz/_UMI.fastq.gz}; 
done

# Trim the adapter and align to Hg38:
for file in Mayer2015*UMI.fastq.gz; do 
  echo $file && 
  STAR --genomeDir hg38_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq ATCTCGTATGCCGT \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rmdir *STARtmp; rm *out *tab

# Sort by coordinates and filter by MAPQ values:
for file in Mayer2015*bam; do 
  echo $file && samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Filter out reads from unwanted genes:
for file in Mayer2015*mapq.bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b hg38_rRNA_tRNA_snRNA_snoRNA.bed > \
    ${file/.bam/_clean.bam}; 
done

# Deduplicate (UMI-Tools):
for file in Mayer2015*clean.bam; do 
  echo $file && 
  samtools index $file && 
  umi_tools dedup --stdin=${file} \
    --stdout=${file/.bam/_dedup.bam}; 
done

# Count output reads:
for file in Mayer2015*dedup.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Load deduplicated BAM files into R session for post-processing (Mammalian_NET-seq_pipeline.R);


##### Process human mNET-seq data from Nojima 2015/2018 -------

# Count input reads:
for file in Nojima2015*fastq.gz Nojima2018_HeLa*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Trim Illumina Small RNA adapter R2 (GATCGTCGGACT) and align human data to Hg38:
for file in Nojima2015*fastq.gz Nojima2016_HeLa*fastq.gz; do 
  echo $file && 
  STAR --genomeDir hg38_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq GATCGTCGGACT \
    --outSAMtype BAM Unsorted; 
done 

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *out *tab *STARtmp

# Filter out reads from unwanted human genes:
for file in Nojima2015*bam Nojima2016_HeLa*bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b hg38_rRNA_tRNA_snRNA_snoRNA.bed > \
    ${file/.bam/_clean.bam}; 
done

# Sort by coordinates and filter by MAPQ values:
for file in Nojima2015*clean.bam Nojima2016_HeLa*clean.bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Count output reads:
for file in Nojima2015*mapq.bam Nojima2016_HeLa*mapq.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Load filtered BAM files into R session for post-processing (Mammalian_NET-seq_pipeline.R);


##### Process mouse mNET-seq data from Nojima 2018 ------------

# Trim Illumina Small RNA adapter R2 (GATCGTCGGACT) and align mouse data to mm10:
for file in Nojima2018_TAP*fastq.gz; do 
  echo $file && 
  STAR --genomeDir mm10_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq GATCGTCGGACT \
    --outSAMtype BAM Unsorted; 
done

# Filter out reads from unwanted genes:
for file in Nojima2018_TAP*bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b mm10_rRNA_tRNA_snRNA_snoRNA.bed > \
    ${file/.bam/_clean.bam}; 
done

# Sort by coordinates and filter by MAPQ values:
for file in Nojima2018_TAP*clean.bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Load filtered BAM files into R session for post-processing (Mammalian_NET-seq_pipeline.R);
