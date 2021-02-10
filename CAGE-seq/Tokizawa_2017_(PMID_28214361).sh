# Tokizawa et al. 2017 (PMID 28214361)
# DRA004921
# "Identification of Arabidopsis genic and non-genic promoters by paired-end sequencing of TSS tags"
# The single CAGE-seq library was prepared from a mix of A.thaliana tissues and sequenced in PE mode;

# Download FASTQ files:
prefix="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004921/DRX060399"
wget ${prefix}/DRR066436_1.fastq.bz2 \
  -O Tokizawa2017_R1.fastq.bz2
wget ${prefix}/DRR066436_2.fastq.bz2 \
  -O Tokizawa2017_R1.fastq.bz2

# Align to the TAIR10 genome:
prefix="Tokizawa2017_" &&
STAR --genomeDir tair10_star \
  --readFilesIn ${prefix}R1.fastq.bz2 ${prefix}R2.fastq.bz2 \
  --runThreadN 4 --outFileNamePrefix ${prefix} \
  --outSAMmultNmax 1 --alignEndsType Local \
  --readFilesCommand bzcat --outSAMtype BAM Unsorted \
  --clip3pAdapterSeq AGATCGGAAGAGC

rm *out *tab; rmdir *STARtmp
for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done

# Extract R1 aligned reads:
samtools view -h -f 66 Tokizawa2017.bam > Tokizawa2017_R1.bam

# Sort R1 BAM file:
samtools view -hu Tokizawa2017_R1.bam | 
samtools sort - -o Tokizawa2017_R1_sorted.bam

# Filter by MAPQ values:
samtools view -h -q 10 Tokizawa2017_R1_sorted.bam \
  -o Tokizawa2017_R1_sorted_mapq.bam; done

# Generate Fw and Rev Bedgraph files from first base of R1:
file="Tokizawa2017_R1_sorted_mapq.bam" && 
bedtools genomecov -ibam $file -bg -5 -strand + | 
sort -k1,1 -k2,2n > ${file/sorted_mapq.bam/fw.bg} && 
bedtools genomecov -ibam $file -bg -5 -strand - | 
sort -k1,1 -k2,2n > ${file/sorted_mapq.bam/rev.bg}

# Merge forward and reverse Bedgraph files:
file1="Tokizawa2017_R1_fw.bg" && 
file2=${file1/fw/rev} && 
outfile=${file1/fw.bg/fw_rev.bedgraph.gz} && 
echo $file1 "+" $file2 "=" $outfile && 
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
cat $file1 - | sort -k1,1 -k2,2n | 
sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
gzip > $outfile

