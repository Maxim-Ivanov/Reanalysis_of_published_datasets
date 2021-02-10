# Zhu et al., 2017 (PMID 28407323)
# "Proteogenomic analysis reveals alternative splicing and translation as part of the abscisic acid response in Arabidopsis seedlings"
# A.thaliana Col-0 seedlings
# Clontech lib prep protocol (SMARTer PCR cDNA Synthesis Kit)

# Install IsoSeq3.1 (the latest version which supports RSII data) and bax2bam:
conda create -n isoseq31 python=2.7
source activate isoseq31
conda install -c conda-forge biopython
conda install -c bioconda bx-python isoseq3=3.1 pbccs=3.4 pbcoretools bax2bam

# Download BAX and BAS files:
prefix="https://sra-pub-src-1.s3.amazonaws.com/SRR5298104/m161031_124550_42199_c100941222550000001823217706101620_s1_X0"

for suffix in "1.bax" "2.bax" "3.bax" "bas"; do 
  wget ${prefix}"."${suffix}.h5.1 \
    -O PacBio_Zhu2017.${suffix}.h5.1; 
done

# Convert BAX files to subread BAM:
bax2bam PacBio_Zhu2017.?.bax.h5.1 -o PacBio_Zhu2017

# Generate CCS:
ccs --noPolish --minPasses 1 PacBio_Zhu2017.subreads.bam PacBio_Zhu2017.css.bam

# Sequences of Clontech primers:
echo -e ">5p
AAGCAGTGGTATCAACGCAGAGTACATGGGG
>3p
GTACTCTGCGTTGATACCACTGCTT" > Clontech_primers.fa

# Classify full-length reads:
lima --isoseq --dump-clips --peek-guess PacBio_Zhu2017.css.bam Clontech_primers.fa PacBio_Zhu2017.fl.bam

for file in *5p--3p*; do mv $file ${file/.5p--3p/}; done

# Refine FL reads:
isoseq3 refine --require-polya PacBio_Zhu2017.fl.bam Clontech_primers.fa PacBio_Zhu2017.flnc.bam

# Cluster FLNC reads:
isoseq3 cluster PacBio_Zhu2017.flnc.bam PacBio_Zhu2017.unpolished.bam

# Polish:
isoseq3 polish PacBio_Zhu2017.unpolished.bam PacBio_Zhu2017.subreads.bam PacBio_Zhu2017.isoforms.bam

# Align clustered HQ transcripts to TAIR10:
minimap2 -t 8 -ax splice -uf --cs --secondary=no \
  --splice-flank=no -C5 -O6,24 -B4 TAIR10.fa \
  PacBio_Zhu2017.isoforms.hq.fastq.gz > \
  PacBio_Zhu2017.isoforms.sam

# Also align the unpolished FLNC reads:
bedtools bamtofastq -i PacBio_Zhu2017.flnc.bam -fq PacBio_Zhu2017.flnc.fastq

minimap2 -t 8 -ax splice -uf --cs --secondary=no \
  --splice-flank=no -C5 -O6,24 -B4 TAIR10.fa \
  PacBio_Zhu2017.flnc.fastq > \
  PacBio_Zhu2017.flnc.sam

# Convert SAM to sorted BAM:
for file in PacBio_Zhu2017*sam; do 
  echo $file && 
  samtools view -hu $file | 
  samtools sort - -o ${file/.sam/_sorted.bam} && 
  rm $file; 
done

# Filter out unmapped reads and reads with low MAPQ:
for file in PacBio_Zhu2017*sorted.bam; do 
  echo $file && 
  samtools view -hb -q 10 -F 4 $file > ${file/.bam/_mapq.bam}; 
done

# Count aligned reads:
for file in PacBio_Zhu2017*mapq.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Convert to BED:
for file in PacBio_Zhu2017*mapq.bam; do 
  echo $file && 
  bedtools bamtobed -i $file > ${file/bam/bed}; 
done
