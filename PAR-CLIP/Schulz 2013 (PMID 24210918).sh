# Schulz et al. 2013 (PMID 24210918);
# "Transcriptome surveillance by selective termination of noncoding RNA synthesis"
# PAR-CLIP protocol: crosslinking sites (detected as T-to-C conversion events) denote genomic positions where the protein of interest binds to its cognate RNA;
# In this study, PAR-CLIP was done for transcription termination complex Nrd1/Nab3;
# Finding T-to-C conversion events require to use an SNP calling pipeline;

# Download PAR-CLIP FASTQ files:
echo -e "Nrd1_rep1\tERR318750/ERR318750
Nrd1_rep2\tERR318753/ERR318753
Nab3_rep1\tERR318748/ERR318748
Nab3_rep2\tERR318746/ERR318746" > Schulz2013_acc.txt

awk '{cmd="wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR318/"$2".fastq.gz -O Schulz2013_PAR_CLIP_"$1".fastq.gz"; system(cmd)}' Schulz2013_acc.txt

# Count raw reads:
for file in Schulz2013_PAR_CLIP*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Trim Illumina Small RNA adapters:
for file in Schulz2013_PAR_CLIP*fastq.gz; do 
  echo $file && 
  trim_galore --small_rna --no_report_file $file; 
done

# Align to SacCer3 genome using STAR:
for file in Schulz2013_PAR_CLIP*fq.gz; do 
  echo $file && 
  STAR --genomeDir saccer3_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/fq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Local \
    --readFilesCommand zcat; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Convert unsorted SAM to sorted BAM:
for file in Schulz2013_PAR_CLIP*sam; do 
  echo $file && 
  samtools view -hu $file | 
  samtools sort - -o ${file/.sam/_sorted.bam}; 
done

# Count aligned reads:
for file in Schulz2013_PAR_CLIP*sorted.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Skip unwanted reads from non-RNAPII transcribed genes:
for file in Schulz2013_PAR_CLIP*sorted.bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b tRNA_rRNA_snRNA_snoRNA_genes_sacCer3.bed > \
    ${file/.bam/_clean.bam}; 
done

# MAPQ filtering:
for file in Schulz2013_PAR_CLIP*clean.bam; do 
  echo $file && 
  samtools view -hbq 10 $file > ${file/.bam/_mapq.bam}; 
done

# Split BAM files by strand orientation into "forward" and "reverse" counterparts: 
for file in Schulz2013_PAR_CLIP*mapq.bam; do 
  echo $file && 
  samtools view -hb -F 16 $file > \
    ${file/sorted_clean_mapq.bam/fwReads.bam} && 
  samtools view -hb -f 16 $file > \
    ${file/sorted_clean_mapq.bam/revReads.bam}; 
done

# Produce mpileup files for VarScan (Nrd1 and Nab3 samples are pooled together):
samtools mpileup -C 50 -f sacCer3.fa \
  Schulz2013_PAR_CLIP*fwReads.bam \
  -o Schulz2013_PAR_CLIP_all_samples_fwReads.mpileup

samtools mpileup -C 50 -f sacCer3.fa \
  Schulz2013_PAR_CLIP*revReads.bam \
  -o Schulz2013_PAR_CLIP_all_samples_revReads.mpileup

# Call "SNPs":
for file in Schulz2013_PAR_CLIP*mpileup; do 
  echo $file && 
  java -Xmx8G -jar VarScan.v2.4.3.jar mpileup2snp $file \
    --min-coverage 6 --min-reads2 2 --min-var-freq 0.01 \
    --min-avg-qual 5 --p-value 0.05 --strand-filter 0 \
    --output-vcf 1 > ${file/mpileup/vcf}; 
done

# Extract conversion events:
cat Schulz2013_PAR_CLIP_all_samples_fwReads.vcf | 
sed '/^#/d' | 
awk '{if ($4=="T" && $5=="C") print}' > \
  Schulz2013_PAR_CLIP_all_samples_fwReads_convTC.vcf

cat Schulz2013_PAR_CLIP_all_samples_revReads.vcf | 
sed '/^#/d' | 
awk '{if ($4=="A" && $5=="G") print}' > \
  Schulz2013_PAR_CLIP_all_samples_revReads_convAG.vcf

# Extract counts of supporting reads to separate Bedgraph files:
python3 Extract_PAR-CLIP_signal_from_VCF.py \
  Schulz2013_PAR_CLIP_all_samples_fwReads_convTC.vcf \
  "Schulz2013_PAR_CLIP" "fw"

python3 Schulz_2013_Extract_PAR-CLIP_signal_from_VCF.py \
  Schulz2013_PAR_CLIP_all_samples_revReads_convAG.vcf \
  "Schulz2013_PAR_CLIP" "rev"

# Combine "fw" and "rev" Bedgraph files belonging to the same samples:
f_str="fw"; 
r_str="rev"; 
ext=".bedgraph"; 
for file1 in Schulz2013_PAR_CLIP*${f_str}${ext}; do 
  file2=${file1/${f_str}/${r_str}} && 
  outfile=${file1/${f_str}/fw_rev} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' > $outfile; 
done

