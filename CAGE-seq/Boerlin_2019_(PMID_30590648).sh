# Boerlin et al. 2019 (PMID 30590648);
# "Saccharomyces cerevisiae displays a stable transcription start site landscape in multiple conditions"
# S.cerevisiae (industrial strain CEN.PK113-7D);
# nAnT-iCAGE protocollib prep protocol;

# Download FASTQ files:
echo -e "008/ERR2495148/ERR2495148\tAna1_rep1
009/ERR2495149/ERR2495149\tAna_rep
000/ERR2495150/ERR2495150\tEth_rep1
001/ERR2495151/ERR2495151\tEth_rep2
002/ERR2495152/ERR2495152\tGlu_rep1
003/ERR2495153/ERR2495153\tGlu_rep2
004/ERR2495154/ERR2495154\tNit_rep1
005/ERR2495155/ERR2495155\tNit_rep2" > Boerlin2019_acc.txt

awk '{cmd="wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR249/"$1".fastq.gz \
  -O Boerlin2019_"$2".fastq.gz"; \
  system(cmd)}' Boerlin2019_acc.txt

# Download CENPK113-7D genome from Supplementary data to Boerlin 2019;

# Generate STAR index for the CENPK113-7D genome:
STAR --runMode genomeGenerate \
  --genomeFastaFiles CENPK113-7D_Sequence.fasta \
  --genomeDir cenpk113_star

# Align CAGE-seq reads to the CENPK113-7D genome:
for file in Boerlin2019*fastq.gz; do 
  echo $file && 
  STAR --genomeDir cenpk_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --clip3pAdapterSeq AGATCGGAAGAGC --outSAMmultNmax 1 \
    --alignEndsType Extend5pOfRead1 --readFilesCommand zcat \
    --outSAMtype BAM Unsorted; 
done

rm *out *tab; rmdir *STARtmp
for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done

# Sort BAM files and remove low MAPQ reads:
for file in Boerlin2019*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Make stranded Bedgraph files (from 5' bases only):
for str in "+" "-"; do 
  [ "$str" = "-" ] && n="rev" || n="fw"; 
  for file in Boerlin2019*mapq.bam; do 
    sample=${file/_mapq.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 -strand $str | 
    sort -k1,1 -k2,2n > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files for the same sample:
f_str="fw"; 
r_str="rev"; 
ext=".bg"; 
for file1 in Boerlin2019*${f_str}${ext}; do 
  file2=${file1/${f_str}/${r_str}} && 
  outfile=${file1/${f_str}${ext}/fw_rev.bedgraph.gz} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $outfile; 
done
