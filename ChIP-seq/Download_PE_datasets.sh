##### Zhu et al. 2015 (PMID 26373455) -------------------------
# "Genome-Wide Prediction and Validation of Intergenic Enhancers in Arabidopsis Using Open Chromatin Signatures"

echo -e "SRR1509476\tbuds_H3K4me1
SRR1509479\tbuds_H3K27ac
SRR1509478\tbuds_H3K27me3
SRR1509477\tleaf_H3K4me1
SRR1509474\tleaf_H3K27ac
SRR1509472\tleaf_H3K27me3" > Zhu2015_acc.txt

awk '{cmd="fastq-dump --split-files --gzip "$1" && \
  mv "$1"_1.fastq.gz PE_Zhu2015_"$2"_R1.fastq.gz && \
  mv "$1"_2.fastq.gz PE_Zhu2015_"$2"_R2.fastq.gz"; \
  system(cmd)}' Zhu2015_acc.txt


##### Cortijo et al. 2017 (PMID 28893714) ---------------------
# "Transcriptional Regulation of the Ambient Temperature Response by H2A.Z Nucleosomes and HSF1 Transcription Factors in Arabidopsis"

echo -e "SRR3234455\tH3_rep1
SRR3234456\tH3_rep2
SRR3234414\tH2AZ_rep1
SRR3234415\tH2AZ_rep2
SRR3234416\tH2AZ_rep3
SRR3234467\tH2B" > Cortijo2017_acc.txt

awk '{cmd="fastq-dump --split-files --gzip "$1" && \
  mv "$1"_1.fastq.gz PE_Cortijo2017_"$2"_R1.fastq.gz && \
  mv "$1"_2.fastq.gz PE_Cortijo2017_"$2"_R2.fastq.gz"; \
  system(cmd)}' Cortijo2017_acc.txt


##### Dai et al. 2017 (PMID 28951178) -------------------------
# "H2A.Z Represses Gene Expression by Modulating Promoter Nucleosome Structure and Enhancer Histone Modifications in Arabidopsis"
# SRP120232

echo -e "SRR6184655\tH3K4me3
SRR6184648\tH3K27me3
SRR6184654\tH2A.Z" > Dai2017_acc.txt

awk '{cmd="fastq-dump --split-files --gzip "$1" && \
  mv "$1"_1.fastq.gz PE_Dai2017_"$2"_R1.fastq.gz && \
  mv "$1"_2.fastq.gz PE_Dai2017_"$2"_R2.fastq.gz"; \
  system(cmd)}' Dai2017_acc.txt

