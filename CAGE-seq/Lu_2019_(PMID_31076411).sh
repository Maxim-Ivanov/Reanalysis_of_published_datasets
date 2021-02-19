# Lu a. Lin 2019 (PMID 31076411);
# "Pervasive and Dynamic Transcription Initiation in Saccharomyces cerevisiae";
# SRP155983;
# Single end sequencing;
# BY4741 strain;
# nAnT-iCAGE lib prep protocol;
# The first 9 bases must be (3 sample barcode + 6 UMI) but apparently they are not;

# Download FASTQ files:
echo -e "SRR7633070\tLu2019_wt_ypd_rep1
SRR7633069\tLu2019_wt_ypd_rep2" > Lu2019_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz "$2".fastq.gz"; system(cmd)}' Lu2019_acc.txt

# Align to SacCer3:
for file in Lu2019*fastq.gz; do 
  echo $file && 
  STAR --genomeDir saccer3_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --clip3pAdapterSeq AGATCGGAAGAGC --outSAMmultNmax 1 \
    --alignEndsType Extend5pOfRead1 --readFilesCommand zcat \
    --outSAMtype BAM Unsorted; 
done

rm *out *tab; rmdir *STARtmp
for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done

# Sort BAM files and remove low MAPQ reads:
for file in Lu2019*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Make stranded Bedgraph files (no strand switch, only 5' bases):
for str in "+" "-"; do 
  [ "$str" = "-" ] && n="rev" || n="fw"; 
  for file in Lu2019*mapq.bam; do 
    sample=${file/_mapq.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 -strand $str | 
    sort -k1,1 -k2,2n > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files for the same sample:
f_str="fw"; r_str="rev"; ext=".bg"; 
for file1 in Lu2019*fw.bg; do 
  file2=${file1/_fw/_rev} && 
  outfile=${file1/_fw.bg/.bedgraph.gz} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $outfile; 
done
