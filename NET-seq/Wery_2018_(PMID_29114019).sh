# Fission yeast (S. pombe) NET-seq

# Wery et al., 2018 (PMID 29114019)
# Accession GSE72382;
# SE sequencing;
# Custom adapter ATCTCGTATGCCGTCTTCTGCTTG (see Churchman 2012 - PMID 21248844);
# No UMIs, no spike-in controls;

# Download FASTQ files from SRA:
acc="SRR2192581"; fastq-dump --gzip $acc && 
  mv ${acc}.fastq.gz Wery2018_rep1.fastq.gz
acc="SRR2192582"; fastq-dump --gzip $acc && 
  mv ${acc}.fastq.gz Wery2018_rep2.fastq.gz

# Download the S.pombe genome (FASTA) and gene annotation (GFF3), adjust chromosome names:
prefix="https://www.pombase.org/data/genome_sequence_and_features/genome_sequence/Schizosaccharomyces_pombe"

wget ${prefix}_chromosome_I.fa.gz -O Spo2_I.fa.gz
wget ${prefix}_chromosome_II.fa.gz -O Spo2_II.fa.gz
wget ${prefix}_chromosome_III.fa.gz -O Spo2_III.fa.gz
wget ${prefix}_mitochondrial_chromosome.fa.gz -O Spo2_MT.fa.gz

zcat Spo2_I.fa.gz | sed '1d;1i >I' | 
  gzip > temp && mv temp Spo2_I.fa.gz
zcat Spo2_II.fa.gz | sed '1d;1i >II' | 
  gzip > temp && mv temp Spo2_II.fa.gz
zcat Spo2_III.fa.gz | sed '1d;1i >III' | 
  gzip > temp && mv temp Spo2_III.fa.gz
zcat Spo2_MT.fa.gz | sed '1d;1i >MT' | 
  gzip > temp && mv temp Spo2_MT.fa.gz
zcat Spo2_chr[I|II|III|MT].fa.gz > Spo2.fa && 
rm Spo2_chr[I|II|III|MT].fa.gz

wget https://www.pombase.org/data/genome_sequence_and_features/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz -O Spo2.gff.gz

zcat Spo2.gff.gz | 
sed '/mating_type_region/d;\
/chr_II_telomeric_gap/d;\
s/mitochondrial/MT/g' | gzip > temp && 
mv temp Spo2.gff.gz

# Generate STAR index for Spo2.fa:
STAR --runMode genomeGenerate --genomeFastaFiles Spo2.fa \
  --runThreadN 4 --genomeDir spo2

# Extract coordinates of rRNA, tRNA, snRNA and snoRNA from GFF3 annotation to BED file:
zcat Spo2.gff3.gz | 
awk 'BEGIN{OFS="\t"}{match($9, /ID=[^;][^;]*/); \
  if ($3=="tRNA" || $3=="rRNA" || $3=="snRNA" || $3=="snoRNA") \
  print $1,$4,$5,substr($9,RSTART+3,RLENGTH-3),$3,$7}' > \
  Spo2_short_ncRNA.bed

# Count input reads:
for file in Wery2018*fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | \
    awk '{print $1}') / 4 )); 
done

# Align to Spo2:
for file in Wery2018*fastq.gz; do 
  echo $file && 
  STAR --genomeDir spo2 --readFilesIn $file --runThreadN 4 \
    --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat \
    --clip3pAdapterSeq ATCTCGTATGCCGT --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rmdir *STARtmp; rm *out *tab

# Sort by coordinates, filter by MAPQ values:
for file in Wery2018*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Filter out reads from unwanted genes:
for file in Wery2018*mapq.bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b Spo2_short_ncRNA.bed > ${file/.bam/_clean.bam}; 
done

# Count output reads:
for file in Wery2018*clean.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Make stranded Bedgraph files (with strand switch):
for str in "+" "-"; do 
  [ "$str" = "+" ] && n="rev" || n="fw"; 
  for file in Shetty2017*dedup.bam; do 
    sample=${file/_mapq_clean_dedup.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 \
      -strand $str > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files for the same sample:
f_str="fw"; 
r_str="rev"; 
ext=".bg"; 
for file1 in Shetty2017*${f_str}${ext}; do 
  file2=${file1/${f_str}/${r_str}} && 
  outfile=${file1/${f_str}${ext}/fw_rev.bedgraph} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' > $outfile; 
done



