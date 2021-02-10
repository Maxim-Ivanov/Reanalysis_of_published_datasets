# Liu a. Jacobsen 2018 (PMID 29379150);
# "RNA-directed DNA methylation involves co-transcriptional small-RNA-guided slicing of polymerase V transcripts in Arabidopsis"
# TruSeq Small RNA lib prep kit;
# PE sequencing (first base of R1 = RNAPII position). Only R1 reads were aligned;

# Zhu et al. 2018 (PMID 30374093);
# "RNA polymerase II activity revealed by GRO-seq and pNET-seq in Arabidopsis"
# Homemade protocol (Wang 2011 - PMID 21572438) which gives flanks identical to Illumina Small RNA-Seq kit (adapter = TGGAATTCTCGG);
# SE sequencing;

# Download FASTQ files, remove R2 files from Liu2018:
echo -e "SRR5681055\tLiu2018_nrpd1e1_rep1
SRR5681056\tLiu2018_nrpd1e1_rep2" > Liu2018_acc.txt

awk '{cmd="fastq-dump --gzip --split-files "$1" && \
  mv "$1"_1.fastq.gz "$2"_R1.fastq.gz && \
  mv "$1"_2.fastq.gz "$2"_R2.fastq.gz"; \
  system(cmd)}' Liu2018_acc.txt

rm Liu2018*R2.fastq.gz
for file in Liu2018*R1.fastq.gz; do mv $file ${file/_R1/}; done

echo -e "SRR6661079\tZhu2018_rep1
SRR6661080\tZhu2018_rep2
SRR7518304\tZhu2018_rep3
SRR7518305\tZhu2018_rep4" > Zhu2018_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz "$2".fastq.gz"; system(cmd)}' Zhu2018_acc.txt


# Count reads in FASTQ files:
for file in Liu2018*fastq.gz Zhu2018*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Align to TAIR10 (in SE mode):
for file in Liu2018*fastq.gz Zhu2018*fastq.gz; do 
  echo $file && 
  STAR --genomeDir tair10_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq TGGAATTCTCGG \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *out *tab *STARtmp *STARgenome

# Count reads in BAM files:
for file in Liu2018*bam Zhu2018*bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}' ); 
done

# Remove intersections with rRNA and tRNA genes:
for file in Liu2018*bam Zhu2018*bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b rRNA_tRNA_U1_Araport11_ext100bp.bed > \
    ${file/.bam/_clean.bam}; 
done

# Sort and filter by MAPQ values:
for file in Liu2018*clean.bam Zhu2018*clean.bam; do 
  echo $file && 
  samtools view -hu -q 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Make stranded Bedgraph files (no strand switch; 5' terminal bases only):
for str in "+" "-"; do 
  [ "$str" = "+" ] && n="fw" || n="rev"; 
  for file in Liu2018*mapq.bam Zhu2018*mapq.bam; do 
    sample=${file/_clean_mapq.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 \
      -strand $str > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files for the same sample:
f_str="fw"; r_str="rev"; ext=".bg"; 
for file1 in Liu2018*${f_str}${ext} Zhu2018*${f_str}${ext}; do 
  file2=${file1/${f_str}/${r_str}} && 
  outfile=${file1/${f_str}${ext}/fw_rev.bedgraph.gz} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $outfile; 
done


