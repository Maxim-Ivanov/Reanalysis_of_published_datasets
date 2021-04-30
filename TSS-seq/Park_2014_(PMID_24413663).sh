# Park 2014 (PMID 24413663)
# "Simultaneous mapping of transcript ends at single-nucleotide resolution and identification of widespread promoter-associated non-coding RNA governed by TATA elements"
# Lib prep protocol: SMORE-seq
# 3' adapter: Illumina Universal


# Download FASTQ files:
echo -e "SRR941531\trep1
SRR941533\trep2
SRR941535\trep3" > Park2014_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Park2014_"$2".fastq.gz"; \
  system(cmd)}' Park2014_acc.txt

# Align to SacCer3 genome:
for file in Park*fastq.gz; do 
  echo $file && 
  STAR --genomeDir saccer3_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq AGATCGGAAGAG \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Filter by MAPQ values and sort BAM files:
for file in Park*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Merge replicates:
samtools merge Park2014_merged_mapq.bam Park2014_rep*bam

# Make stranded Bedgraph files (without strand switch):
for str in "+" "-"; do 
  [ "$str" = "-" ] && 
  n="rev" || n="fw"; 
  for file in Park*mapq.bam; do 
    sample=${file/_mapq.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 -strand $str | 
    sort -k1,1 -k2,2n > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files for the same sample:
f_str="fw"; r_str="rev"; ext=".bg"; 
for file1 in *${f_str}${ext}; 
  do file2=${file1/${f_str}/${r_str}} && 
  outfile=${file1/${f_str}${ext}/fw_rev.bedgraph.gz} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $outfile; 
done

