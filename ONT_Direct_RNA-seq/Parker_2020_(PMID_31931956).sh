# ONT Direct RNA-seq of A.thaliana seedlings:

# 1) WT samples:
# Parker et al., 2020 (PMID 31931956)
# "Nanopore direct RNA sequencing maps the complexity of Arabidopsis mRNA processing and m6A modification"
# https://www.ebi.ac.uk/ena/data/view/PRJEB32782

# 2) hen2-2 samples (exonuclease-deficient):
# Parker et al., 2020
# "Widespread premature transcription termination of Arabidopsis thaliana NLR genes by the spen protein FPA"
# https://doi.org/10.1101/2020.12.15.422694
# https://www.biorxiv.org/content/10.1101/2020.12.15.422694v1
# Accession PRJEB41381


##### Download FASTQ files from Parker 2020a -----------------

echo -e "ERR3764345/col0_nanopore_drs_1\trep1
ERR3764347/col0_nanopore_drs_2a\trep2a
ERR3764348/col0_nanopore_drs_2b\trep2b
ERR3764349/col0_nanopore_drs_3\trep3
ERR3764351/col0_nanopore_drs_4\trep4" > Parker2020a_acc.txt

awk '{cmd="wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/"$1".tar.gz \
  -O Parker2020a_"$2".tar.gz", system(cmd)}' Parker2020a_acc.txt

# Extract subfolders with FASTQ files:
for file in Parker2020a*tar.gz; do 
  tar xvf $file --wildcards '*fastq/pass*'; 
done

# Concatenate FASTQ files belonging to the same sample:
find ./*FAH45730* -type f -name "*fastq" \
  -exec cat {} > Parker2020a_rep1.fastq \;
find ./*FAH77434* -type f -name "*fastq" \
  -exec cat {} > Parker2020a_rep2_part1.fastq \;
find ./*FAH59362* -type f -name "*fastq" \
  -exec cat {} > Parker2020a_rep2_part2.fastq \;
cat Parker2020a_rep2_part*.fastq > Parker2020a_rep2.fastq && 
rm Parker2020a_rep2_part*.fastq
find ./*FAH83697* -type f -name "*fastq" \
  -exec cat {} > Parker2020a_rep3.fastq \;
find ./*FAH83552* -type f -name "*fastq" \
  -exec cat {} > Parker2020a_rep4.fastq \;
for file in *fastq; do gzip $file; done

##### Download FASTQ file from Parker 2020b --------------------

wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR484/ERR4844543/hen2-2_nanopore_drs.tar.gz -O Parker2020b_hen2.tar.gz

# Extract subfolders with FASTQ files:
tar xvf Parker2020b_hen2.tar.gz --wildcards '*fastq/pass*'

# Concatenate all FASTQ files for the hen2-2 sample:
cat ./cluster/gjb_lab/cdr/ON_MinION_datastore/20181114_1642_20181114_mRNA_hen2-2_2322/fastq/pass/*fastq | 
gzip > Parker2020b_hen2.fastq.gz


##### Align and post-process the long reads -------------------

# Count raw reads:
for file in Parker2020*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done


# Align to TAIR10 genome using Minimap2:
for file in Parker2020*fastq.gz; do 
  echo $file && 
  minimap2 -t 8 -ax splice -k 14 -L --cs --secondary=no \
    -G 10000 TAIR10.fa $file > ${file/fastq.gz/sam}; 
done


# Convert SAM to sorted BAM:
for file in Parker2020*sam; do 
  echo $file && 
  samtools view -hu $file | 
  samtools sort - -o ${file/sam/bam}; 
done

# Filter out unmapped reads and reads with low MAPQ:
for file in Parker2020*bam; do 
  echo $file && 
  samtools view -hb -q 10 -F 4 $file > ${file/.bam/_mapq.bam}; 
done

# Count aligned reads:
for file in Parker2020*mapq.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Convert to BED:
#for file in *mapq.bam; do echo $file && bedtools bamtobed -i $file > ${file/_mapq.bam/.bed}; done

# Post-process filtered BAM files in R session using Filter_Direct_RNA-seq_BAM_files.R;
