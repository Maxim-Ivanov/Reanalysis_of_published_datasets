# Capovilla et al., 2018 (PMID 29988152)
# "PORCUPINE regulates development in response to temperature through alternative splicing"
# Accession PRJEB24412;
# PE sequencing;

# Download FASTQ files:
echo -e "009/ERR2245559/ERR2245559\t9d_16C
002/ERR2245562/ERR2245562\t9d_23C
008/ERR2245568/ERR2245568\t9d_23C_3d_16C_expA_rep1
009/ERR2245569/ERR2245569\t9d_23C_3d_16C_expA_rep2
005/ERR2245575/ERR2245575\t9d_23C_3d_27C_expA_rep1
006/ERR2245576/ERR2245576\t9d_23C_3d_27C_expA_rep2
007/ERR2245577/ERR2245577\t12d_23C_expA_rep1
008/ERR2245578/ERR2245578\t12d_23C_expA_rep2" > Capovilla2018_acc.txt

awk -v BASE="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/" '{cmd="wget BASE"$1"_1.fastq.gz \
  -O Capovilla2018_"$2"_R1.fq.gz && \
  wget BASE"$1"_2.fastq.gz \
  -O Capovilla2018_"$2"_R2.fq.gz"; \
  system(cmd)}' Capovilla2018_acc.txt

# Count raw reads:
for file in Capovilla2018*R1.fq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Download Araport11 annotation: 
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz

# Fix chromosome names:
zcat Araport11_GFF3_genes_transposons.201606.gff.gz | 
  sed 's/^ChrC/Pt/;s/^ChrM/Mt/;s/^Chr//' > Araport11.gff

# Align to TAIR10 genome in transcriptome-guided mode:
for f1 in Capovilla2018*R1.fq.gz; do 
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
for file in Capovilla2018*bam; do 
  echo $file $(samtools flagstat $file | \
    sed -n '9p' | awk '{sprintf("%i", $1/2)}'); 
done

# Filter for reads with high MAPQ in proper pairs:
for file in Capovilla2018*bam; do 
  echo $file && 
  samtools view -hb -f 2 -q 10 $file -o ${file/.bam/_mapq.bam}; 
done

# Run samtools fixmate (required by samtools markdup):
for file in Capovilla2018*mapq.bam; do 
  echo $file && 
  samtools fixmate -m $file ${file/.bam/_fixmate.bam}; 
done

# Sort by coordinates:
for file in Capovilla2018*fixmate.bam; do 
  echo $file && 
  samtools sort $file -o ${file/.bam/_sorted.bam}; 
done

# Deduplicate on start coordinates of paired reads:
for file in Capovilla2018*sorted.bam; do 
  echo $file && 
  samtools markdup -r $file ${file/.bam/_dedup.bam}; 
done

# Load deduplicated BAM files into R session and produce stranded normalized Bedgraph files using convert_RNAseq_BAM_to_normalized_bedGraph.R;
