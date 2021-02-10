# Kohnen et al., 2016 (PMID 27923878)
# "Neighbor Detection Induces Organ-Specific Transcriptomes, Revealing Patterns Underlying Hypocotyl-Specific Growth"
# Accession GSE81202;
# SE sequencing;

# Download single-end FASTQ files:
prefix="SE_RNAseq_Kohnen"
acc="SRR3480142"; fastq-dump --gzip $acc && \
  mv ${acc}.fastq.gz ${prefix}_rep1.fastq.gz
acc="SRR3480144"; fastq-dump --gzip $acc && \
  mv ${acc}.fastq.gz ${prefix}_rep2.fastq.gz

# Count raw reads:
for file in SE_RNAseq*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Download Araport11 annotation: 
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz

# Fix chromosome names:
zcat Araport11_GFF3_genes_transposons.201606.gff.gz | 
  sed 's/^ChrC/Pt/;s/^ChrM/Mt/;s/^Chr//' > Araport11.gff

# Align to TAIR10 in transcriptome-guided mode:
for file in SE_RNAseq*fastq.gz; do 
  echo $file && 
  STAR --genomeDir tair10_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq AGATCGGAAGAGC \
    --outSAMtype BAM Unsorted --sjdbGTFfile Araport11.gff \
    --sjdbGTFtagExonParentTranscript Parent \
    --outSAMstrandField intronMotif; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARgenome *STARtmp *out *tab

# Sort BAM files:
for file in SE_RNAseq*bam; do 
  echo $file && 
  samtools sort $file -o ${file/.bam/_sorted.bam}; 
done

# Count aligned reads:
for file in SE_RNAseq*sorted.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Deduplicate (based on start coordinates):
for file in SE_RNAseq*sorted.bam; do 
  echo $file && 
  samtools rmdup -s $file ${file/.bam/_dedup.bam}; 
done

# Remove low MAPQ reads:
for file in SE_RNAseq*dedup.bam; do echo $file && samtools view -hb -q 10 $file > ${file/.bam/_mapq.bam}; done

# Load filtered BAM files into R session and produce stranded normalized Bedgraph files using convert_RNAseq_BAM_to_normalized_bedGraph.R;
