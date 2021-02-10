# Capovilla et al., 2018 (PMID 29988152)
# "PORCUPINE regulates development in response to temperature through alternative splicing"
# Accession PRJEB24412;
# PE sequencing;

# Hofmann et al., 2019 (PMID 30610360)
# "The embryonic transcriptome of Arabidopsis thaliana"
# Accession GSE121236;
# Unusual lib prep protocol:
# Smart-seq2 cDNA synthesis (Picelli 2013 - PMID 24056875) -> 
# fragmentation -> ligation of Nextera adapters;
# PE sequencing;


# Download Capovilla 2018 data:
echo -e "009/ERR2245559/ERR2245559\t9d_16C
002/ERR2245562/ERR2245562\t9d_23C
008/ERR2245568/ERR2245568\t9d_23C_3d_16C_expA_rep1
009/ERR2245569/ERR2245569\t9d_23C_3d_16C_expA_rep2
005/ERR2245575/ERR2245575\t9d_23C_3d_27C_expA_rep1
006/ERR2245576/ERR2245576\t9d_23C_3d_27C_expA_rep2
007/ERR2245577/ERR2245577\t12d_23C_expA_rep1
008/ERR2245578/ERR2245578\t12d_23C_expA_rep2" > Capovilla2018_acc.txt

awk -v BASE="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/" '{cmd="wget BASE"$1"_1.fastq.gz \
  -O PE_RNAseq_Capovilla2018_"$2"_R1.fq.gz && \
  wget BASE"$1"_2.fastq.gz \
  -O PE_RNAseq_Capovilla2018_"$2"_R2.fq.gz"; \
  system(cmd)}' Capovilla2018_acc.txt


# Download Hofmann 2019 data:
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
  mv "$1"_1.fastq.gz PE_Hofmann2019_"$2"_R1.fastq.gz && \
  mv "$1"_2.fastq.gz PE_Hofmann2019_"$2"_R2.fastq.gz"; \
  system(cmd)}' Hofmann2019_acc.txt


# Count raw reads:
for file in PE_RNAseq*R1.fq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# For Hofmann 2019 data only: trim Nextera adapters from 3' ends, trim poly(A) and poly(T) stretches from 5' ends of both R1 and R2:
for f1 in PE_RNAseq_Hofmann2019*R1.fq.gz; do 
  f2=${f1/R1/R2} &&  
  echo $f1 "+" $f2 && 
  cutadapt -j 4 -m 18 -O 5 -a "CTGTCTCTTATACACATCTG" \
    -A "CTGTCTCTTATACACATCTG" -g "A{100}" -G "A{100}" \
    -g "T{100}" -G "T{100}" -o ${f1/.fq.gz/_trim.fq} \
    -p ${f2/.fq.gz/_trim.fq} <(zcat $f1) <(zcat $f2); 
done

for file in PE_RNAseq_Hofmann2019*trim.fq; do 
  echo $file && gzip $file; 
done


# Download Araport11 annotation: 
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz

# Fix chromosome names:
zcat Araport11_GFF3_genes_transposons.201606.gff.gz | 
  sed 's/^ChrC/Pt/;s/^ChrM/Mt/;s/^Chr//' > Araport11.gff

# Align to TAIR10 genome in transcriptome-guided mode:
for f1 in PE_RNAseq*R1.fq.gz; do 
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
for file in PE_RNAseq*bam; do 
  echo $file $(samtools flagstat $file | \
    sed -n '9p' | awk '{sprintf("%i", $1/2)}'); 
done

# Filter for reads with high MAPQ in proper pairs:
for file in PE_RNAseq*bam; do 
  echo $file && 
  samtools view -hb -f 2 -q 10 $file -o ${file/.bam/_mapq.bam}; 
done

# Run samtools fixmate (required by samtools markdup):
for file in PE_RNAseq*mapq.bam; do 
  echo $file && 
  samtools fixmate -m $file ${file/.bam/_fixmate.bam}; 
done

# Sort by coordinates:
for file in PE_RNAseq*fixmate.bam; do 
  echo $file && 
  samtools sort $file -o ${file/.bam/_sorted.bam}; 
done

# Deduplicate on start coordinates of paired reads:
for file in PE_RNAseq*sorted.bam; do 
  echo $file && 
  samtools markdup -r $file ${file/.bam/_dedup.bam}; 
done

# Load deduplicated BAM files into R session and produce stranded normalized Bedgraph files using convert_RNAseq_BAM_to_normalized_bedGraph.R;