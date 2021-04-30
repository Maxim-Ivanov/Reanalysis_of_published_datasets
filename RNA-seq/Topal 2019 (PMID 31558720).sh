# Topal 2019 (PMID 31558720)
# "Distinct transcriptional roles for Histone H3-K56 acetylation during the cell cycle in Yeast"
# Yeast NET-seq, TT-seq and 4sU-seq
# NET-seq:
# - Churchman's adapter (ATCTCGTATGCCGTCTTCTGCTTG)
# - 5' UMI (6 bp)
# - Mix of S.cerevisiae (90%) + S.pombe (10%) reads


# Download NET-seq FASTQ files:
echo -e "SRR8503038\tNETseq_rep1
SRR8503039\tNETseq_rep2" > Topal2019_SE_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Topal2019_"$2".fastq.gz"; \
  system(cmd)}' Topal2019_SE_acc.txt

# Process UMIs:
for file in Topal2019_NETseq*fastq.gz; do 
  echo $file && 
  umi_tools extract --stdin=${file} --bc-pattern=NNNNNN \
    --stdout=${file/.fastq.gz/_UMI.fq.gz}; 
done

# Align to SacCer3 + Spo2 in SE mode:
# (for details on STAR index for combined SacCer3 + Spo2 genome, see https://github.com/Maxim-Ivanov/Reanalysis_of_published_datasets/blob/main/NET-seq/Shetty_2017_(PMID_28366642).sh)

for file in Topal2019_NETseq*UMI.fq.gz; do 
  echo $file && 
  STAR --genomeDir spo2_saccer3 --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/UMI.fq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq ATCTCGTATGCCG \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab


# Sort by coordinates, filter by MAPQ values, remove Spo2 alignments (contig name do not start with "chr"):
for file in Topal2019_NETseq*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  awk '{if ($3~/chr/ || $2~/chr/) print}' | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Filter out reads from unwanted genes:
# (for details on short ncRNA annotation, see https://github.com/Maxim-Ivanov/Reanalysis_of_published_datasets/blob/main/NET-seq/Scerevisiae_NET-seq.sh)

for file in Topal2019_NETseq*mapq.bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b SacCer3_rRNA_tRNA_genes_ext100bp.bed > ${file/.bam/_clean.bam}; 
done

# Deduplicate on UMIs:
for file in Topal2019_NETseq*clean.bam; do 
  echo $file && 
  samtools index $file && 
  umi_tools dedup --stdin=${file} \
    --stdout=${file/.bam/_dedup.bam}; 
done

# Merge replicates:
samtools merge Topal2019_NETseq_merged_mapq_clean_dedup.bam Topal2019_NETseq_rep*bam

# Make stranded Bedgraph files (with strand switch):
# (whole reads + 5' bases)
for str in "+" "-"; do 
  [ "$str" = "+" ] && n="rev" || n="fw"; 
  for file in Topal2019_NETseq*dedup.bam; do 
    sample=${file/_mapq_clean_dedup.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -split\
      -strand $str > ${sample}_whole_reads_${n}.bg &&
    bedtools genomecov -ibam $file -bg -5\
      -strand $str > ${sample}_first_bases_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files:
f_str="fw"; r_str="rev"; ext=".bg"; 
for file1 in Topal2019_NETseq*${f_str}${ext}; do 
  file2=${file1/${f_str}/${r_str}} && 
  outfile=${file1/${f_str}${ext}/fw_rev.bedgraph.gz} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $outfile; 
done
