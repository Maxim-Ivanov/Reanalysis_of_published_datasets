# Harrison 2014 (PMID 26092945)
# "PAT-seq: a method to study the integration of 3'-UTR dynamics with gene expression in the eukaryotic transcriptome"


# Download FASTQ files:
echo -e "SRR1054315\trep1
SRR1054316\trep2" > Harrison2014_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Harrison2014_"$2".fastq.gz"; \
  system(cmd)}' Harrison2014_acc.txt

# Trim remnants of poly(A) tails from 5' ends and Illumina adapters from 3' ends:
for file in Harrison*fastq.gz; do 
  echo $file && 
  cutadapt -j 4 -q 10 -g "T{150};e=0.05" --max-n 3 $file | 
  cutadapt -a AGATCGGAAGAGC -m 20 \
    -o ${file/.fastq.gz/_trim.fq.gz} -; 
done

# Align to the yeast genome:
for file in Harrison*trim.fq.gz; do 
  echo $file && 
  STAR --genomeDir saccer3_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/trim.fq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq AGATCGGAAGAG \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files and filter by MAPQ values:
for file in Harrison*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Merge replicates:
samtools merge Harrison2014_merged_mapq.bam Harrison2014_rep*bam

# Make stranded Bedgraph files (WITHOUT strand switch):
for str in "+" "-"; do 
  [ "$str" = "-" ] && 
  n="rev" || n="fw"; 
  for file in *mapq.bam; do 
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

