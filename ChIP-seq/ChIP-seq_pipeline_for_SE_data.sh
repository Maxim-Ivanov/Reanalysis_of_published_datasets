# Most often ChIP-seq libraries are sequenced in the SE mode. If insert size exceeds the read length, then pairs of sequencing coverage peaks flank the true origin of ChIP-seq signal. This unwanted behavior can be attenuated by extending SE reads at their 3' ends. The rate of extension should depend on the average insert size, which can be inferred from the ChIP-seq data by MACS2 (see Zhang 2008 - PMID 18798982). 

# This pipeline was designed for SE ChIP-seq data acquired by Download_SE_datasets.sh script;

# Count input reads:
for file in SE*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Trim Illumina adapters and align to TAIR10:

for file in SE*fastq.gz; do 
  echo $file && 
  STAR --genomeDir tair10_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Local \
    --readFilesCommand zcat \
    -- clip3pAdapterSeq "AGATCGGAAGAGC" \
    --outSAMtype BAM Unsorted; 
done

rm *out *tab; rmdir *STARtmp
for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done

# Sort BAM, remove split alignments and MAPQ values below 5:
for file in SE*bam; do 
  echo $file && 
  samtools view -hu -q 5 $file | 
  awk '$6~/N/ {next}{print $0}' | 
  samtools sort - -o ${file/.bam/_sorted.bam}; 
done
# (alternatively, split alignments could be suppressed by adding --alignIntronMax 1 to the STAR call)

# Count the filtered reads:
for file in SE*sorted.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Extend reads to d/2 and generate gzipped bedGraph files using MACS2:
for file in SE*sorted.bam; do 
  echo $file && 
  macs2 -t $file -n ${file/_sorted.bam/} -g 1.35e+08 -m 3,50 \
    --half-ext --bdg > ${file/_sorted.bam/_log.txt} 2>&1; 
done

rm *model.r *pvalue.bdg *qvalue.bdg *control_lambda.bdg *peaks.bed *summits.bed *peaks.xls *encodePeak *log.txt

for file in SE*pileup.bdg; do 
  mv $file ${file/_treat_pileup.bdg/.bedgraph}; 
done

# Sort, round and gzip Bedgraph files:
for file in SE*bedgraph; do 
  echo $file && 
  sort -k1,1 -k2,2n $file | 
  awk '{printf("%s\t%d\t%d\t%d\n", $1,$2,$3,$4)}' | 
  gzip > ${file}.gz; 
done





