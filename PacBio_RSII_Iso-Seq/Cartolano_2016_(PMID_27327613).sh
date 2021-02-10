# Cartolano et al., 2016 (PMID 27327613)
# "cDNA Library Enrichment of Full Length Transcripts for SMRT Long Read Sequencing"
# A.thaliana Col-0 inflorescence tissue
# Clontech lib prep protocol
# TeloPrime Full-Length cDNA Amplification Kit (Lexogen; cap-dependent linker ligation), or SMARTer PCR cDNA Synthesis Kit (Clontech)
# Size fractionation (A = 3-6Kb, B = 2-3Kb, C = 1-2Kb)

# Install IsoSeq3.1 (the latest version which supports RSII data) and bax2bam:
conda create -n isoseq31 python=2.7
source activate isoseq31
conda install -c conda-forge biopython
conda install -c bioconda bx-python isoseq3=3.1 pbccs=3.4 pbcoretools bax2bam

# Download BAX and BAS files:

echo -e "SRR3655761/m150515_010955_42207_c100759272550000001823162807221567_s1_p0\tClontech_C1
SRR3655770/m150508_080413_42207_c100759082550000001823162807221514_s1_p0\tClontech_C2
SRR3655760/m150518_170939_42207_c100759062550000001823162807221530_s1_p0\tClontech_B1
SRR3655769/m150508_034007_42207_c100759082550000001823162807221513_s1_p0\tClontech_B2
SRR3655768/m150507_231731_42207_c100759082550000001823162807221512_s1_p0\tClontech_A1
SRR3655759/m150514_204533_42207_c100759272550000001823162807221566_s1_p0\tClontech_A2
SRR3655767/m150319_170737_42207_c100759662550000001823163407221560_s1_p0\tTeloPrime_C1
SRR3655758/m150317_192824_42207_c100759072550000001823162807221520_s1_p0\tTeloPrime_C2
SRR3655766/m150317_234722_42207_c100759072550000001823162807221521_s1_p0\tTeloPrime_B1
SRR3655757/m150416_080043_42207_c100759252550000001823162807221584_s1_p0\tTeloPrime_B2
SRR3655765/m150416_034136_42207_c100759252550000001823162807221583_s1_p0\tTeloPrime_A1
SRR3655756/m150121_130920_42207_c100718232550000001823147005141590_s1_p0\tTeloPrime_A2" > Cartolano2016_acc.txt

awk '{cmd="for suffix in 1.bax 2.bax 3.bax bas; do \
  wget https://sra-pub-src-1.s3.amazonaws.com/"$1".${suffix}.h5.1 \
  -O PacBio_Cartolano2016_"$2".${suffix}.h5.1; \
  done"; system(cmd)}' Cartolano2016_acc.txt

# Convert BAX files to subread BAM:
for f1 in PacBio_Cartolano2016*.1.bax.h5.1; do 
  prefix=${f1/.1.bax.h5.1/} && 
  bax2bam $f1 ${prefix}.2.bax.h5.1 \
    ${prefix}.3.bax.h5.1 -o $prefix; 
done

# Generate CCS:
for file in PacBio_Cartolano2016*subreads.bam; do 
  echo $file && 
  ccs --noPolish --minPasses 1 \
    $file ${file/subreads.bam/ccs.bam}; 
done

# Sequences of Clontech and TeloPrime primers:
echo -e ">5p
AAGCAGTGGTATCAACGCAGAGTACATGGGG
>3p
GTACTCTGCGTTGATACCACTGCTT" > Clontech_primers.fa

echo -e ">5p
TGGATTGATATGTAATACGACTCACTATAG
>3p
CGCCTGAGA" > TeloPrime_primers.fa

# Classify full-length reads:
for file in PacBio_Cartolano2016_TeloPrime*ccs.bam; do 
  echo $file && 
  lima --isoseq --dump-clips --peek-guess \
    $file TeloPrime_primers.fa ${file/ccs.bam/fl.bam}; 
done

for file in PacBio_Cartolano2016_Clontech*ccs.bam; do 
  echo $file && 
  lima --isoseq --dump-clips --peek-guess \
    $file Clontech_primers.fa ${file/ccs.bam/fl.bam}; 
done

for file in *5p--3p*; do mv $file ${file/.5p--3p/}; done

# Refine FL reads:
for file in PacBio_Cartolano2016_TeloPrime*fl.bam; do 
  echo $file && 
  isoseq3 refine $file TeloPrime_primers.fa ${file/fl/flnc}; 
done
# (observe that --require-polya was skipped for TeloPrime)

for file in PacBio_Cartolano2016_Clontech*fl.bam; do 
  echo $file && 
  isoseq3 refine --require-polya $file \
    Clontech_primers.fa ${file/fl/flnc}; 
done

# Cluster FLNC reads:
for file in PacBio_Cartolano2016*flnc.bam; do 
  echo $file && 
  isoseq3 cluster $file ${file/flnc/unpolished}; 
done

# Polish:
for file in PacBio_Cartolano2016*unpolished.bam; do 
  echo $file && 
  isoseq3 polish $file ${file/unpolished/subreads} \
    ${file/unpolished/isoforms}; 
done

# Align clustered HQ transcripts to TAIR10:
for file in PacBio_Cartolano2016*isoforms.hq.fastq.gz; do 
  echo $file && 
  minimap2 -t 8 -ax splice -uf --cs --secondary=no \
    --splice-flank=no -C5 -O6,24 -B4 TAIR10.fa $file > \
    ${file/hq.fastq.gz/sam}; 
done

# Also align the unpolished FLNC reads:
for file in PacBio_Cartolano2016*flnc.bam; do 
  echo $file && 
  bedtools bamtofastq -i $file -fq ${file/bam/fastq} && 
  gzip ${file/bam/fastq}; 
done

for file in PacBio_Cartolano2016*flnc.fastq.gz; do 
  echo $file && 
  minimap2 -t 8 -ax splice -uf --cs --secondary=no \
    --splice-flank=no -C5 -O6,24 -B4 TAIR10.fa $file > \
    ${file/fastq.gz/sam}; 
done

# Convert SAM to sorted BAM:
for file in PacBio_Cartolano2016*sam; do 
  echo $file && 
  samtools view -hu $file | 
  samtools sort - -o ${file/.sam/_sorted.bam} && 
  rm $file; 
done

# Filter out unmapped reads and reads with low MAPQ:
for file in PacBio_Cartolano2016*sorted.bam; do 
  echo $file && 
  samtools view -hb -q 10 -F 4 $file > ${file/.bam/_mapq.bam}; 
done

# Count aligned reads:
for file in PacBio_Cartolano2016*mapq.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Convert to BED:
for file in PacBio_Cartolano2016*mapq.bam; do 
  echo $file && 
  bedtools bamtobed -i $file > ${file/bam/bed}; 
done



