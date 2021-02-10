# PAT-seq (poly(A) tag sequencing): first base of R1 = position of the poly(A) site;

# Wu et al., 2011 (PMID 21746925)
# "Genome-wide landscape of polyadenylation in Arabidopsis provides evidence for extensive alternative polyadenylation"
# Accession SRA028410;
# PE sequencing;
# "5' reads begin with either CATG or ACGT, and 3' reads begin with the sequence NNNCCTTTT";
# No strand switch;

# Thomas et al., 2012 (PMID 23136375)
# "Genome-wide control of polyadenylation site choice by CPSF30 in Arabidopsis"
# Accession SRA048565;
# PE sequencing;
# No strand switch;

# Yu et al., 2019 (PMID 31427469)
# "Transcriptome Analyses of FY Mutants Reveal Its Role in mRNA Alternative Polyadenylation"
# Accession SRP145554;
# SE sequencing;
# Reads have 8bp barcodes at 5' ends (followed by oligo-T stretches of variable length!) and Illumina adapters on 3' ends;
# Strand switch required;


# Download Wu 2011 data:
prefix="PATseq_Wu2011"
acc="SRR089763"; fastq-dump --gzip --split_files $acc && 
mv ${acc}_1.fastq.gz ${prefix}_rep1_R1.fastq.gz && 
mv ${acc}_2.fastq.gz ${prefix}_rep1_R2.fastq.gz

acc="SRR089776"; fastq-dump --gzip --split_files $acc && 
mv ${acc}_1.fastq.gz ${prefix}_rep2_R1.fastq.gz && 
mv ${acc}_2.fastq.gz ${prefix}_rep2_R2.fastq.gz

# Download Thomas 2012 data:
echo -e "SRR389031\trep1
SRR389028\trep2
SRR389036\trep3" > Thomas2012_acc.txt

awk '{cmd="fastq-dump --gzip --split-files "$1" && \
  mv "$1"_1.fastq.gz PATseq_Thomas2012_"$2"_R1.fastq.gz && \
  mv "$1"_2.fastq.gz PATseq_Thomas2012_"$2"_R2.fastq.gz"; \
  system(cmd)}' Thomas2012_acc.txt

# Download Yu 2019 data:
echo -e "007/SRR7160297/SRR7160297\trep1
006/SRR7160296/SRR7160296\trep2
009/SRR7160299/SRR7160299\trep3" > Yu2019_acc.txt

awk '{cmd="wget \
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR716/"$1".fastq.gz \
  -O PATseq_Yu2019_"$2".fastq.gz"; system(cmd)}' Yu2019_acc.txt

# Trim first 4 bases from R1 in Wu 2011 data:
for file in PATseq_Wu2011*R1.fastq.gz; do 
  echo $file && 
  cutadapt -u 4 -m 20 --max-n 3 -o temp $file &&
  mv temp $file; 
done

# Trim 5' barcodes and oligo-T stretches from Yu2019 data:
for file in PATseq_Yu2018*fastq.gz; do 
  echo $file && 
  cutadapt -j 4 -u 8 -q 10 -g "T{150};e=0.05" --max-n 3 \
    -m 20 -o temp $file && 
  mv temp $file; 
done

# Count reads:
for file in PATseq_Wu2011*R1.fastq.gz \
PATseq_Thomas2012*R1.fastq.gz \
PATseq_Yu2018*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Align R1 reads to the TAIR10 genome using STAR:
for file in PATseq_Wu2011*R1.fastq.gz \
PATseq_Thomas2012*R1.fastq.gz \
PATseq_Yu2018*fastq.gz; do 
  echo $file && 
  STAR --genomeDir tair10_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq AGATCGGAAGAGC \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Count aligned reads:
for file in PATseq*bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Sort BAM files and remove low MAPQ reads:
for file in PATseq*bam; do 
  echo $file && 
  samtools view -hu -q 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done



# Make stranded Bedgraph files for Wu 2011 and Thomas 2012 (no strand switch):
for str in "+" "-"; do 
  [ "$str" = "-" ] && n="rev" || n="fw"; 
  for file in PATseq_Wu2011*mapq.bam \
  PATseq_Thomas2012*mapq.bam; do 
    sample=${file/_mapq.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 \
      -strand $str > ${sample}_${n}.bg; 
  done; 
done

# Make stranded Bedgraph files for Yu2019 (with strand switch):
for str in "+" "-"; do 
  [ "$str" = "-" ] && n="fw" || n="rev"; 
  for file in PATseq_Yu2019*mapq.bam; do 
    sample=${file/_mapq.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 \
      -strand $str > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files for the same sample:
f_str="fw"; 
r_str="rev"; 
ext=".bg"; 
for file1 in *${f_str}${ext}; do 
  file2=${file1/${f_str}/${r_str}} && 
  outfile=${file1/${f_str}${ext}/fw_rev.bedgraph.gz} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | gzip > $outfile; 
done

