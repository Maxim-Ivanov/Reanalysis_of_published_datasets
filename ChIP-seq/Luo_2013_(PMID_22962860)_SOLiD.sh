# Luo et al., 2013 (PMID 22962860)
# "Integrative analysis of chromatin states in Arabidopsis identified potential regulatory mechanisms for natural antisense transcript production"
# Color space data from SOLiD!
# Conversion from csFASTQ to Sanger FASTQ resulted in a very low fraction of aligned reads. Most probably, this happened due to high frequency of sequencing errors (which propagate when csfastq files are converted to FASTQ);
# Thus, this dataset was aligned in color space mode using Bowtie1;

# Download the csFASTQ files:

echo -e "SRR037789\tH3K4me2_rep1
SRR037790\tH3K4me2_rep2
SRR037791\tH3K4me3_rep1
SRR037792\tH3K4me3_rep2
SRR037879\tH3K9Ac_rep1
SRR037879\tH3K9Ac_rep2
SRR037793\tH3K9me2_rep1
SRR037794\tH3K9me2_rep2
SRR037710\tH3K18Ac
SRR037788\tH3K27me1
SRR037881\tH3K27me3
SRR037785\tH3K36me2
SRR037786\tH3K36me3
SRR037787\tH3_rep1
SRR037882\tH3_rep2
SRR037883\tH3_rep3
SRR037795\tInput_rep1
SRR037796\tInput_rep2
SRR037884\tInput_rep3
SRR037885\tInput_rep4" > Luo2013_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".csfastq.gz Luo2013_"$2".csfastq.gz"; \
  system(cmd)}' Luo2013_acc.txt

# Build color space index of TAIR10 genome with Bowtie:
bowtie-build -C TAIR10.fa tair10_cs

# Run Bowtie (v1.2.2):
for file in Luo2013*csfastq.gz; do 
  echo $file && 
  bowtie -C -p 4 --best -S tair10_cs $file ${file/.csfastq.gz/.sam}; 
done

# Remove unmapped and low MAPQ reads, secondary and split alignments, convert SAM to sorted BAM:
for file in Luo2013*sam; do echo $file && samtools view -h -F 260 -q 5 $file | awk '$6~/N/ {next}{print $0}' | samtools sort - -o ${file/.sam/_sorted.bam}; done

# Count filtered reads:
for file in Luo2013*sorted.bam; do echo $file && samtools flagstat $file | sed -n '1p' | awk '{print $1}'; done

# Extend reads to d/2 and generate gzipped bedGraph files using MACS2:
for file in Luo2013*sorted.bam; do echo $file && macs2 -t $file -n ${file/_sorted.bam/} -g 1.35e+08 -m 3,50 --half-ext --bdg > ${file/_sorted.bam/_log.txt} 2>&1; done

rm *model.r *pvalue.bdg *qvalue.bdg *control_lambda.bdg *peaks.bed *summits.bed *peaks.xls *encodePeak *log.txt

for file in Luo2013*pileup.bdg; do 
  mv $file ${file/_treat_pileup.bdg/.bedgraph}; 
done

# Sort, round and gzip Bedgraph files:
for file in Luo2013*bedgraph; do 
  echo $file && 
  sort -k1,1 -k2,2n $file | 
  awk '{printf("%s\t%d\t%d\t%d\n", $1,$2,$3,$4)}' | 
  gzip > ${file}.gz; 
done



