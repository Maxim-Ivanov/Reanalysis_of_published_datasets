##### Inagaki et al. 2017 (PMID 28100676) ---------------------
# "Gene-body chromatin modification dynamics mediate epigenome differentiation in Arabidopsis"

echo -e "DRX066754/DRR072816\tH3_rep1
DRX066769/DRR072831\tH3_rep2
DRX066784/DRR072846\tH3_rep3
DRX066799/DRR072861\tH3_rep4
DRX066755/DRR072817\tH3K4me1_rep1
DRX066770/DRR072832\tH3K4me1_rep2
DRX066785/DRR072847\tH3K4me1_rep3
DRX066800/DRR072862\tH3K4me1_rep4
DRX066756/DRR072818\tH3K4me2_rep1
DRX066771/DRR072833\tH3K4me2_rep2
DRX066786/DRR072848\tH3K4me2_rep3
DRX066801/DRR072863\tH3K4me2_rep4
DRX066757/DRR072819\tH3K4me3_rep1
DRX066772/DRR072834\tH3K4me3_rep2
DRX066787/DRR072849\tH3K4me3_rep3
DRX066802/DRR072864\tH3K4me3_rep4
DRX066758/DRR072820\tH3K9me2_rep1
DRX066773/DRR072835\tH3K9me2_rep2
DRX066788/DRR072850\tH3K9me2_rep3
DRX066803/DRR072865\tH3K9me2_rep4" > Inagaki2017_acc.txt

awk '{cmd="wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA005/DRA005154/"$1".fastq.bz2 \
  -O SE_Inagaki2017_"$2".fastq.bz2"; \
  system(cmd)}' Inagaki2017_acc.txt

# Convert bz2 to gz:
for file in SE_Inagaki2017*bz2; do 
  echo $file && 
  bzcat $file | gzip > ${file/bz2/gz}; 
done


##### Chen et al. 2017 (PMID 28947800) ------------------------
# "Cytosolic acetyl-CoA promotes histone acetylation predominantly at H3K27 in Arabidopsis"
# GSE73972

echo -e "SRR5224481\tH4K16ac_rep1
SRR3286792\tH4K16ac_rep2
SRR5224480\tH4K12ac_rep1
SRR3286791\tH4K12ac_rep2
SRR5224479\tH4K8ac_rep1
SRR3286790\tH4K8ac_rep2
SRR5224478\tH4K5ac_rep1
SRR3286789\tH4K5ac_rep2
SRR5224477\tH3K18ac_rep1
SRR3286788\tH3K18ac_rep2
SRR5224476\tH3K14ac_rep1
SRR3286787\tH3K14ac_rep2
SRR5224475\tH3K9ac_rep1
SRR3286786\tH3K9ac_rep2
SRR5224474\tH3K27ac_rep1
SRR3286783\tH3K27ac_rep2
SRR3286801\tInput" > Chen2017_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz SE_Chen2017_"$2".fastq.gz"; \
  system(cmd)}' Chen2017_acc.txt


##### Yelagandula et al. 2014 (PMID 24995981) -----------------
# "The histone variant H2A.W defines heterochromatin and promotes chromatin condensation in Arabidopsis"

echo -e "SRR988541\tInput
SRR988542\tH3
SRR988543\tH2A
SRR988545\tH2AX
SRR988544\tH2AW
SRR988546\tH2AZ" > Yelagandula2014_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz SE_Yelagandula2017_"$2".fastq.gz"; \
  system(cmd)}' Yelagandula2014_acc.txt


##### Gomez-Zambrano et al. 2018 (PMID 29604400) --------------
# Arabidopsis SWC4 Binds DNA and Recruits the SWR1 Complex to Modulate Histone H2A.Z Deposition at Key Regulatory Genes

echo -e "SRR6785140\tInput_rep1
SRR6785143\tInput_rep2
SRR6785139\tH3_rep1
SRR6785136\tH3_rep2
SRR6785146\tH2AZ_rep1
SRR6785135\tH2AZ_rep2" > Gomez2018_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz SE_Gomez2018_"$2".fastq.gz"; \
  system(cmd)}' Gomez2018_acc.txt

##### Zhou et al. 2017 (PMID 28403905) -------------------------
# "H2A monoubiquitination in Arabidopsis thaliana is generally independent of LHP1 and PRC2 activity"

echo -e "SRR4734671\tInput
SRR4734656\tH2Aub_rep1
SRR4734658\tH2Aub_rep2
SRR4734664\tH2Aub_rep3
SRR4734666\tH2Aub_rep4
SRR4734668\tH2Aub_rep5
SRR4734670\tH2Aub_rep6
SRR4734660\tH3K27me3_rep1
SRR4734662\tH3K27me3_rep2
SRR5278091\tH3K27me3_rep3
SRR5278093\tH3K27me3_rep4" > Zhou2017_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz SE_Zhou2017_"$2".fastq.gz"; \
  system(cmd)}' Zhou2017_acc.txt


##### Nassrallah et al. 2018 (PMID 30192741) ------------------
# DET1-mediated degradation of a SAGA-like deubiquitination module controls H2Bub homeostasis

echo -e "SRR6989572\tInput_rep1
SRR6989573\tInput_rep2
SRR6989580\tH2Bub_rep1
SRR6989581\tH2Bub_rep2" > Nassrallah2018_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz SE_Nassrallah2018_"$2".fastq.gz"; \
  system(cmd)}' Nassrallah2018_acc.txt


##### Torres a. Deal 2019 (PMID 30742338) ---------------------
# "The histone variant H2A.Z and chromatin remodeler BRAHMA act coordinately and antagonistically to regulate transcription and nucleosome dynamics in Arabidopsis"

echo -e "SRR6412323\tInput_rep1
SRR6412324\tInput_rep2
SRR6412321\tH2AZ_rep1
SRR6412322\tH2AZ_rep2" > Torres2019_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz SE_Torres2019_"$2".fastq.gz"; \
  system(cmd)}' Torres2019_acc.txt


##### Wang et al. 2015 (PMID 26100864) ------------------------
# "Osmotic stress induces phosphorylation of histone H3 at threonine 3 in pericentromeric regions of Arabidopsis thaliana"
# GSE68370

echo -e "SRR2001264\tH3_rep1
SRR2001265\tH3_rep2
SRR2001266\tH3_rep3
SRR2001267\tH3_rep4
SRR2001268\tH3K4me1_rep1
SRR2001269\tH3K4me1_rep2
SRR2001270\tH3K4me1_rep3
SRR2001271\tH3K4me3_rep1
SRR2001272\tH3K4me3_rep2
SRR2001273\tH3K4me3_rep3
SRR2001274\tH3T3ph_rep1
SRR2001275\tH3T3ph_rep2
SRR2001276\tH3T3ph_rep3
SRR2001277\tH3T3ph_rep4" > Wang2015_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz SE_Wang2015_"$2".fastq.gz"; \
  system(cmd)}' Wang2015_acc.txt


##### Bewick et al. 2016 (PMID 27457936) ----------------------
# "On the origin and evolutionary consequences of gene body DNA methylation"

echo -e "SRR3087131\tH3K56ac
SRR3087130\tH3K36me3
SRR3087129\tH3K27me3
SRR3087128\tH3K9me2
SRR3087127\tH3K4me3" > Bewick2016_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz SE_Bewick2016_"$2".fastq.gz"; \
  system(cmd)}' Bewick2016_acc.txt


##### Trejo-Arellano et al. 2017 (PMID 28182313) --------------
# "H3K23me1 is an evolutionarily conserved histone modification associated with CG DNA methylation in Arabidopsis"
# GSE86498

echo -e "SRR4184990\tInput_rep1
SRR4184991\tInput_rep1
SRR4184994\tH3_rep1
SRR4184995\tH3_rep1
SRR4184998\tH3K23me1_rep1
SRR4184999\tH3K23me1_rep1" > Trejo2017_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz SE_Trejo2017_"$2".fastq.gz"; \
  system(cmd)}' Trejo2017_acc.txt


##### Zhang et al. 2016 (PMID 27208288) -----------------------
# "The Second Subunit of DNA Polymerase Delta Is Required for Genomic Stability and Epigenetic Regulation"
# GSE79259

echo -e "SRR3228527\tInput
SRR3228526\tH3
SRR3228524\tH3K4me3
SRR3228523\tH3K27me3
SRR3228525\tH3K9me2" > Zhang2016_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && '
  mv "$1".fastq.gz SE_Zhang2016_"$2".fastq.gz"; \
  system(cmd)}' Zhang2016_acc.txt


##### Wollmann et al. 2017 (PMID 28521766) --------------------
# "The histone H3 variant H3.3 regulates gene body DNA methylation in Arabidopsis thaliana"
# GSE96834

echo -e "SRR5364419\tH3
SRR5364420\tH3
SRR5364425\tH3K4me3
SRR5364426\tH3K4me3
SRR5364421\tH2A.Z
SRR5364422\tH2A.Z
SRR5364427\tH3K36me3
SRR5364428\tH3K36me3" > Wollmann2017_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz SE_Wollmann2017_"$2".fastq.gz"; \
  system(cmd)}' Wollmann2017_acc.txt
