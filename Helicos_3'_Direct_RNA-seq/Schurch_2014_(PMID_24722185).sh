# Schurch et al. 2014 (PMID 24722185);
# "Improved annotation of 3' untranslated regions and complex loci by combination of strand-specific direct RNA sequencing, RNA-Seq and ESTs"
# ERP003245
# Direct DNA-seq on the Helicos platform (dead since 2012);
# First bases of reads = positions of polyadenylation sites (PAS);
# FASTQ files contain only dummy quality scores;

# Download Helicos FASTQ from DDBJ:
prefix="ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/ERA223/ERA223202"
wget ${prefix}/ERX267337/ERR294004.fastq.bz2 \
  -O Schurch2014_rep1.fastq.bz2
wget ${prefix}/ERX267338/ERR294005.fastq.bz2 \
  -O Schurch2014_rep2.fastq.bz2
wget ${prefix}/ERX267339/ERR294006.fastq.bz2 \
  -O Schurch2014_rep3.fastq.bz2

# Remove trailing white spaces after read names and convert to GZ:
for file in Schurch2014*fastq.bz2; do 
  echo $file && 
  bzcat $file | sed 's/ $//' | gzip > ${file/bz2/gz} && 
  rm $file; 
done

# Align reads to TAIR10 genome by STAR:
for file in Schurch2014*fastq.gz; do 
  echo $file && 
  STAR --genomeDir tair10_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Local \
    --readFilesCommand zcat --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files and remove low MAPQ reads:
for file in Schurch2014*bam; do 
  echo $file && 
  samtools view -hu -q 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Make strand-specific Bedgraph files (use only the first base of each read and switch the strand orientation):
for str in "+" "-"; do 
  [ "$str" = "+" ] && n="rev" || n="fw"; 
  for file in Schurch2014*mapq.bam; do 
    sample=${file/_mapq.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 -strand $str | 
    sort -k1,1 -k2,2n > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files:
f_str="fw"; r_str="rev"; ext=".bg"; 
for f1 in Schurch2014*${f_str}${ext}; do 
  f2=${f1/${f_str}/${r_str}} && 
  out=${f1/${f_str}${ext}/fw_rev.bedgraph.gz} && 
  echo $f1 "+" $f2 "=" $out && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $f2 | 
  cat $f1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $out; 
done
