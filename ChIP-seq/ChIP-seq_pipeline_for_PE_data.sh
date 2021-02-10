# Sometimes ChIP-seq libraries are sequenced in PE mode. This greatly simplifies the downstream analysis, since coordinates of the original DNA fragment can be obtained directly from the BAM file (i.e. without inferring the average insert size by MACS2);

# This pipeline was designed for PE ChIP-seq data acquired by Download_PE_datasets.sh script;

# Count input read pairs:
for file in PE*_R1.fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Trim adapters and align to TAIR10:
for f1 in PE*_R1.fastq.gz; do 
  f2=${f1/_R1.f/_R2.f} && 
  echo $f1 $f2 && 
  STAR --genomeDir tair10_star --readFilesIn $f1 $f2 \
    --runThreadN 4 --outFileNamePrefix ${f1/R1.fastq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Local \
    --readFilesCommand zcat \
    --clip3pAdapterSeq "AGATCGGAAGAGC" \
    --outSAMtype BAM Unsorted; 
done

rm *out *tab; rmdir *STARtmp
for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done

# Sort BAM files, remove split alignments, MAPQ values below 5 and singletons:
for file in PE*bam; do 
  echo $file && 
  samtools view -hu -q 5 -f 2 $file | 
  awk '$6~/N/ {next}{print $0}' | 
  samtools sort - -o ${file/.bam/_sorted.bam}; 
done

# Count the filtered read pairs:
for file in PE*sorted.bam; do 
  echo $file $(samtools flagstat $file | sed -n '7p' | \
    awk '{print $1}'); 
done

# Get forward reads in proper pairs -> add TLEN (positive) to POS (start) to get the end position of the insert;
for file in PE*sorted.bam; do 
  echo $file && 
  samtools view -f 2 -F 16 $file | 
  awk '{print $3"\t"$4"\t"$4+$9}' > \
    ${file/sorted.bam/insert.bed}; 
done

# Convert BED files to Bedgraph:
for file in PE*insert.bed; do 
  echo $file && 
  bedtools genomecov -i $file \
    -g tair10_star/chrNameLength.txt -bg | 
  sort -k1,1 -k2,2n | gzip > ${file/_insert.bed/.bedgraph.gz}; 
done
