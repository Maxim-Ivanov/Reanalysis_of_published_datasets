# Qui 2020 (PMID 32487207)
# "Universal promoter scanning by Pol II during transcription initiation in Saccharomyces cerevisiae"
# 3' adapter: SOLiD Small RNA adapter (CGCCTTGGCCGT)


# Download FASTQ files:
echo -e "SRR8580779\trep1
SRR8580780\trep2" > Qui2020_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Qui2020_"$2".fastq.gz"; \
  system(cmd)}' Qui2020_acc.txt

# Align to SacCer3 genome:
for file in Qui*fastq.gz; do 
  echo $file && 
  STAR --genomeDir saccer3_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq CGCCTTGGCCGT \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Filter by MAPQ scores and sort:
for file in Qui*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Merge replicates:
samtools merge Qui2020_merged_mapq.bam Qui2020_rep*bam

# Make stranded Bedgraph files (without strand switch):
for str in "+" "-"; do 
  [ "$str" = "-" ] && 
  n="rev" || n="fw"; 
  for file in Qui*mapq.bam; do 
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

