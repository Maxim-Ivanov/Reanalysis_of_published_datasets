# Hofmann et al., 2019 (PMID 30610360)
# "The embryonic transcriptome of Arabidopsis thaliana"
# Accession GSE121236;
# Unusual lib prep protocol:
# Smart-seq2 cDNA synthesis (Picelli 2013 - PMID 24056875) -> 
# fragmentation -> ligation of Nextera adapters;
# Illumina PE sequencing;
# Unstranded data;

# Download FASTQ files:
echo -e "SRR8054359\ts01_globular_rep1
SRR8054360\ts02_globular_rep2
SRR8054361\ts03_globular_rep3
SRR8054362\ts04_earlyHeart_rep1
SRR8054363\ts05_earlyHeart_rep2
SRR8054364\ts06_earlyHeart_rep3
SRR8054365\ts07_lateHeart_rep1
SRR8054366\ts08_lateHeart_rep2
SRR8054367\ts09_lateHeart_rep3
SRR8054368\ts10_earlyTorpedo_rep1
SRR8054369\ts11_earlyTorpedo_rep2
SRR8054370\ts12_earlyTorpedo_rep3
SRR8054371\ts13_lateTorpedo_rep1
SRR8054372\ts14_lateTorpedo_rep2
SRR8054373\ts15_lateTorpedo_rep3" > Hofmann2019_acc.txt

awk '{cmd="fastq-dump --gzip --split-files "$1" && \
  mv "$1"_1.fastq.gz Hofmann2019_"$2"_R1.fastq.gz && \
  mv "$1"_2.fastq.gz Hofmann2019_"$2"_R2.fastq.gz"; \
  system(cmd)}' Hofmann2019_acc.txt

# Count raw reads:
for file in Hofmann2019*R1.fq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Trim Nextera adapters from 3' ends, trim poly(A) and poly(T) stretches from 5' ends of both R1 and R2:
for f1 in Hofmann2019*R1.fq.gz; do 
  f2=${f1/R1/R2} &&  
  echo $f1 "+" $f2 && 
  cutadapt -j 4 -m 18 -O 5 -a "CTGTCTCTTATACACATCTG" \
    -A "CTGTCTCTTATACACATCTG" -g "XA{100}" -G "XA{100}" \
    -g "XT{100}" -G "XT{100}" -o ${f1/.fq.gz/_trim.fq} \
    -p ${f2/.fq.gz/_trim.fq} <(zcat $f1) <(zcat $f2); 
done

for file in Hofmann2019*trim.fq; do 
  echo $file && gzip $file; 
done

# Download Araport11 annotation: 
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz

# Fix chromosome names:
zcat Araport11_GFF3_genes_transposons.201606.gff.gz | 
  sed 's/^ChrC/Pt/;s/^ChrM/Mt/;s/^Chr//' > Araport11.gff

# Align to TAIR10 genome in transcriptome-guided mode:
for f1 in Hofmann2019*R1.fq.gz; do 
  f2=${f1/_R1/_R2} && 
  echo $f1 $f2 && 
  STAR --genomeDir tair10_star --readFilesIn $f1 $f2 \
    --readFilesCommand zcat --runThreadN 4 \
    --outFileNamePrefix ${f1/R1.fq.gz/} --outSAMmultNmax 1 \
    --alignEndsType Local --clip3pAdapterSeq AGATCGGAAGAGC \
    --outSAMtype BAM Unsorted --sjdbGTFfile Araport11.gff \
    --sjdbGTFtagExonParentTranscript Parent \
    --outSAMstrandField intronMotif; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARgenome *STARtmp *out *tab

# Count aligned read pairs:
for file in Hofmann2019*bam; do 
  echo $file $(samtools flagstat $file | \
    sed -n '9p' | awk '{sprintf("%i", $1/2)}'); 
done

# Filter for reads with high MAPQ in proper pairs:
for file in Hofmann2019*bam; do 
  echo $file && 
  samtools view -hb -f 2 -q 10 $file -o ${file/.bam/_mapq.bam}; 
done

# Run samtools fixmate (required by samtools markdup):
for file in Hofmann2019*mapq.bam; do 
  echo $file && 
  samtools fixmate -m $file ${file/.bam/_fixmate.bam}; 
done

# Sort by coordinates:
for file in Hofmann2019*fixmate.bam; do 
  echo $file && 
  samtools sort $file -o ${file/.bam/_sorted.bam}; 
done

# Deduplicate on start coordinates of paired reads:
for file in Hofmann2019*sorted.bam; do 
  echo $file && 
  samtools markdup -r $file ${file/.bam/_dedup.bam}; 
done

# Rename deduplicated files:
for file in Hofmann2019*dedup.bam; do mv $file ${file/mapq_fixmate_sorted_dedup/final}; done

# Merge replicates:
samtools merge globular_merged_final.bam s{01..03}*final.bam
samtools merge earlyHeart_merged_final.bam s{04..06}*final.bam
samtools merge lateHeart_merged_final.bam s{07..09}*final.bam
samtools merge earlyTorpedo_merged_final.bam s{10..12}*final.bam
samtools merge lateTorpedo_merged_final.bam s{13..15}*final.bam

# Make unstranded Bedgraph files (normalized to 1M fragments):
for file in Hofmann2019*final.bam; do 
  echo $file && 
  norm=$(samtools flagstat $file | sed -n '9p' | \
    awk '{print 2000000/$1}') && 
  echo $norm && 
  bedtools genomecov -ibam $file -bg -split -scale $norm | 
  gzip > ${file/.bam/_norm1M.bedgraph.gz}; 
done