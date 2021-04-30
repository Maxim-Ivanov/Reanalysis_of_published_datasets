# Malabat et al., 2015 (PMID 25905671)
# "Quality control of transcription start site selection by nonsense-mediated-mRNA decay"
# GSE64139
# 5' adapter sequence: NNNNCGCGCGNN(G)
# 3' adapter: Illumina Universal (AGATCGGAAGAG)


# Download FASTQ files:
echo -e "SRR1708689\trep1
SRR1708697\trep2
SRR1708701\trep3
SRR1708693\trep4
SRR1708705\trep5
SRR1708709\trep6" > Malabat2015_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Malabat2015_"$2".fastq.gz"; \
  system(cmd)}' Malabat2015_acc.txt

# Process UMIs:
for file in Malabat2015*fastq.gz; do 
  echo $file && 
  umi_tools extract --stdin=${file} --bc-pattern=NNNNXXXXXXNNX \
    --stdout=${file/.fastq.gz/_UMI.fq.gz}; 
done

# Align to SacCer3 with STAR:
for file in Malabat2015*UMI.fq.gz; do 
  echo $file && 
  STAR --genomeDir ${index_sc} --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/UMI.fq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq AGATCGGAAGAG \
    --clip5pNbases 7 --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files and filter by MAPQ values:
for file in Malabat2015*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Deduplicate on UMIs:
for file in Malabat2015*mapq.bam; do 
  echo $file && 
  samtools index $file && 
  umi_tools dedup --stdin=${file} \
    --stdout=${file/.bam/_dedup.bam}; 
done

# Merge replicates:
samtools merge Malabat2015_merged_mapq_dedup.bam Malabat2015_rep*dedup.bam

# Make stranded Bedgraph files (without strand switch):
for str in "+" "-"; do 
  [ "$str" = "-" ] && 
  n="rev" || n="fw"; 
  for file in Malabat*dedup.bam; do 
    sample=${file/_mapq_dedup.bam/} && 
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
