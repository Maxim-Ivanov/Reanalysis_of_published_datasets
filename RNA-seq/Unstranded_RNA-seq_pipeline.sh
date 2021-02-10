# Calixto et al., 2018 (PMID 29764987)
# "Rapid and Dynamic Alternative Splicing Impacts the Arabidopsis Cold Response Transcriptome"
# Accession PRJEB19974;
# PE sequencing;

# Huertas et al., 2019 (PMID 30696706)
# "Arabidopsis SME1 Regulates Plant Development and Response to Abiotic Stress by Determining Spliceosome Activity Specificity"
# Accession GSE116964;
# PE sequencing;


# Download Calixto 2018 data:
echo -e "008/ERR1886188/ERR1886188\twt_rep1
004/ERR1886284/ERR1886284\twt_rep2
003/ERR1886353/ERR1886353\twt_rep3
009/ERR1886189/ERR1886189\tcold3h_rep1
005/ERR1886285/ERR1886285\tcold3h_rep2
004/ERR1886354/ERR1886354\tcold3h_rep3
002/ERR1886192/ERR1886192\tcold12h_rep1
008/ERR1886288/ERR1886288\tcold12h_rep2
007/ERR1886357/ERR1886357\tcold12h_rep3" > Calixto2018_acc.txt

awk -v BASE="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/" '{cmd="wget BASE"$1"_1.fastq.gz \
  -O Uns_RNAseq_Calixto2018_"$2"_R1.fq.gz && \
  wget BASE"$1"_2.fastq.gz \
  -O Uns_RNAseq_Calixto2018_"$2"_R2.fq.gz"; \
  system(cmd)}' Calixto2018_acc.txt

# Download Huertas 2019 data:
echo -e "SRR7511615\twt
SRR7511616\tsme1
SRR7511617\twt_cold
SRR7511618\tsme1_cold" > Huertas2019_acc.txt

awk '{cmd="fastq-dump --gzip --split-files "$1" && \
  mv "$1"_1.fastq.gz Uns_RNAseq_Huertas2019_"$2"_R1.fq.gz && \
  mv "$1"_2.fastq.gz Uns_RNAseq_Huertas2019_"$2"_R2.fq.gz"; \
  system(cmd)}' Huertas2019_acc.txt

# Count raw reads:
for file in Uns*R1.fq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Download Araport11 annotation: 
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz

# Fix chromosome names:
zcat Araport11_GFF3_genes_transposons.201606.gff.gz | 
  sed 's/^ChrC/Pt/;s/^ChrM/Mt/;s/^Chr//' > Araport11.gff

# Align to TAIR10 genome in transcriptome-guided mode:
for f1 in Uns_RNAseq*R1.fq.gz; do 
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
for file in Uns_RNAseq*bam; do 
  echo $file $(samtools flagstat $file | \
    sed -n '9p' | awk '{sprintf("%i", $1/2)}'); 
done

# Filter for reads with high MAPQ in proper pairs:
for file in Uns_RNAseq*bam; do 
  echo $file && 
  samtools view -hb -f 2 -q 10 $file -o ${file/.bam/_mapq.bam}; 
done

# Run samtools fixmate (required by samtools markdup):
for file in Uns_RNAseq*mapq.bam; do 
  echo $file && 
  samtools fixmate -m $file ${file/.bam/_fixmate.bam}; 
done

# Sort by coordinates:
for file in Uns_RNAseq*fixmate.bam; do 
  echo $file && 
  samtools sort $file -o ${file/.bam/_sorted.bam}; 
done

# Deduplicate on start coordinates of paired reads:
for file in Uns_RNAseq*sorted.bam; do 
  echo $file && 
  samtools markdup -r $file ${file/.bam/_dedup.bam}; 
done

# Make unstranded Bedgraphs (normalized to 1M fragments):
for file in Uns_RNAseq*bam; do 
  echo $file && 
  norm=$(samtools flagstat $file | sed -n '9p' | \
    awk '{print 2000000/$1}') && 
  echo $norm && 
  bedtools genomecov -ibam $file -bg -split -scale $norm | 
  gzip > ${file/mapq_fixmate_sorted_dedup.bam/norm1M.bedgraph.gz}; 
done
