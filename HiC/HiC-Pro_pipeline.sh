# This script demonstrates the usage of HiC-Pro on Arabidopsis HiC data;

# Grob et al., 2014 (PMID 25132176);
# "Hi-C analysis in Arabidopsis identifies the KNOT, a structure with similarities to the flamenco locus of Drosophila"
# Accession GSE55960

# Feng et al., 2014 (PMID 25132175);
# "Genome-wide Hi-C analyses in wild-type and mutants reveal high-resolution chromatin interactions in Arabidopsis"
# Accession SRP043612;

# Download FASTQ files:
acc="SRR1197490"; fastq-dump --gzip --split-files $acc && \
  mv ${acc}_1.fastq.gz Grob2014_R1.fastq.gz && \
  mv ${acc}_2.fastq.gz Grob2014_R2.fastq.gz

acc="SRR1504819"; fastq-dump --gzip --split-files $acc && \
  mv ${acc}_1.fastq.gz Feng2014_R1.fastq.gz && \
  mv ${acc}_2.fastq.gz Feng2014_R2.fastq.gz

# Prepare the Anaconda environment (Python 2.7):
conda create --name hicpro_env python=2.7 pysam bx-python scipy
source activate hicpro
conda install -c r r-base r-ggplot2 r-rcolorbrewer

# Download HiC-Pro:
wget https://github.com/nservant/HiC-Pro/archive/v2.11.1.tar.gz
tar -zxvf v2.11.1.tar.gz
rm v2.11.1.tar.gz
cd HiC-Pro-2.11.1
make configure
make install

# Make Bowtie2 index of TAIR10 genome:
bowtie2-build TAIR10.fa tair10_bt2 

# Generate list of restriction fragments:
./bin/utils/digest_genome.py -r hindiii \
  -o ./annotation/tair10_HindIII.bed TAIR10.fa

# Copy the chromosome size file:
samtools faidx TAIR10.fa | 
cut -f1,2 > ./annotation/tair10.chrom.sizes

# Edit the configuration file:
nano config-hicpro.txt
# PAIR1_EXT = _R1
# PAIR2_EXT = _R2
# BOWTIE2_IDX_PATH = tair10_bt2
# REFERENCE_GENOME = TAIR10
# GENOME_SIZE = HiC-Pro-2.11.1/annotation/tair10.chrom.sizes
# GENOME_FRAGMENT = HiC-Pro-2.11.1/annotation/tair10_HindIII.bed

# Run the full HiC-Pro pipeline:
./bin/HiC-Pro -i ${fastq_dir} -o ${output_dir} \
  -c config-hicpro.txt

# Generate a viewpoint (e.g. "2:8122000-8123076"):
for sample in Feng2014 Grob2014; do 
  echo $sample && 
  ./bin/utils/make_viewpoints.py \
    -i ${output_dir}/hic_results/data/${sample}/ \
    ${sample}.allValidPairs \
    -f ./annotation/tair10_HindIII.bed -t target.bed -v \
    -o ${output_dir}/viewpoint_${sample}.bedgraph; 
done

# Sort and gzip Bedgraph (contains header line):
for file in viewpoint*bedgraph; do 
  header=$(cat $file | head -n 1) && 
  sed '1d' $file | sort -k1,1 -k2,2n | cat $header - | 
  gzip > ${file}.gz; 
done

