# MNase-seq signal is proportional to the nucleosome occupancy (i.e. MNase-seq peaks are expected to coincide with nucleosome positions). This is in contrast to DNase-seq (where increased signal corresponds to nucleosome-free regions);
# MNase data are usually PE;
# MNase fragments have to be filtered for length around 140-150 bp (to enrich for the mono-nucleosomal fraction);

# Download MNase-seq data (A.thaliana seedlings) from 3 datasets:
# Dai et al. 2017 (PMID 28951178)
# Cortijo et al. 2017 (PMID 28893714)
# Torres a. Deal 2019 (PMID 30742338)

echo -e "SRR6184650\tDai2017
SRR3234432\tCortijo2017_rep1
SRR3234433\tCortijo2017_rep2
SRR6412345\tTorres2019_rep1
SRR6412346\tTorres2019_rep2" > MNase_acc.txt

awk '{cmd="fastq-dump --split-files --gzip "$1" && \
  mv "$1"_1.fastq.gz MNase_"$2"_R1.fastq.gz && \
  mv "$1"_2.fastq.gz MNase_"$2"_R2.fastq.gz"; \
  system(cmd)}' MNase_acc.txt

# Count input read pairs:
for file in MNase*_R1.fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Trim adapters and align to TAIR10:

for f1 in MNase*_R1.fastq.gz; do 
  f2=${f1/_R1.f/_R2.f} && 
  echo $f1 $f2 && 
  STAR --genomeDir tair10_star --readFilesIn $f1 $f2 \
    --runThreadN 4 --outFileNamePrefix ${f1/R1.fastq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Local \
    --readFilesCommand zcat \
    --clip3pAdapterSeq "AGATCGGAAGAGC" \
    --outSAMtype BAM Unsorted; 
done

rm *out *tab; rmdir *STARtmp
for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done

# Sort BAM files, remove split alignments, MAPQ values below 5 and singletons:
for file in MNase*bam; do 
  echo $file && 
  samtools view -hu -q 5 -f 2 $file | 
  awk '$6~/N/ {next}{print $0}' | 
  samtools sort - -o ${file/.bam/_sorted.bam}; 
done

# Count the filtered read pairs:
for file in MNase*sorted.bam; do 
  echo $file $(samtools flagstat $file | sed -n '7p' | \
    awk '{print $1}'); 
done

# Get forward reads in proper pairs -> add TLEN (positive) to POS (start) to get the end position of the insert. Also filter for insert size between 140 and 155:
for file in PE*sorted.bam; do 
  echo $file && 
  samtools view -f 2 -F 16 $file | 
  awk '{if ($9 >= 140 && $9 <= 155) print $3"\t"$4"\t"$4+$9}' > \
  ${file/sorted.bam/insert.bed}; 
done

# Convert BED files to Bedgraph:
for file in Mnase*insert.bed; do 
  echo $file && 
  bedtools genomecov -i $file -g tair10_star/chrNameLength.txt \
    -bg | sort -k1,1 -k2,2n | 
  gzip > ${file/_insert.bed/.bedgraph.gz}; 
done
