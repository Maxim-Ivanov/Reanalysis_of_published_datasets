# Schmid 2018 (PMID 30157438)
# "Simultaneous Measurement of Transcriptional and Post-transcriptional Parameters by 3' End RNA-Seq"
# Lib prep protocol: Lexogen QuantSeq 3' mRNA-seq

# Remapped samples:
# SRR6423292	noPap_Nab2AA_input_0_1
# SRR6423296	noPap_Nab2AA_ip_0_1
# SRR6423300	xPap_Nab2AA_input_0_1
# SRR6423304	xPap_Nab2AA_ip_0_1
# SRR6423305	xPap_Nab2AA_ip_0_2
# SRR6423306	xPap_Nab2AA_ip_0_3
# SRR6423294	noPap_Mex67AA_input_0_1
# SRR6423298	noPap_Mex67AA_ip_0_1
# SRR6423302	xPap_Mex67AA_input_0_1
# SRR6423310	xPap_Mex67AA_ip_0_1
# SRR6423311	xPap_Mex67AA_ip_0_2
# SRR6423312	xPap_Mex67AA_ip_0_3

# The "xPAP" samples were treated with poly(A)-polymerase prior to library construction (in order to detect termination sites of non-polyadenylated lncRNAs). As a side effect, products of RNA fragmentation were also polyadenylated, which resulted in smooth signal over gene bodies. Thus, it is hard to distinguish between true intragenic alternative pA sites and the artifactual RNA degradation signal in "xPAP" samples.
# The "noPAP" samples are the normal samples where the whole signal is expected to denote pA sites.


# Download and rename FASTQ files:
echo -e "SRR6423292\tstr2_noPAP_total
SRR6423296\tstr2_noPAP_4sU
SRR6423300\tstr2_withPAP_total
SRR6423304\tstr2_withPAP_4sU_rep1
SRR6423305\tstr2_withPAP_4sU_rep2
SRR6423306\tstr2_withPAP_4sU_rep3
SRR6423294\tstr1_noPAP_total
SRR6423298\tstr1_noPAP_4sU
SRR6423302\tstr1_withPAP_total
SRR6423310\tstr1_withPAP_4sU_rep1
SRR6423311\tstr1_withPAP_4sU_rep2
SRR6423312\tstr1_withPAP_4sU_rep3" > Schmid2018_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz Schmid2018_"$2".fastq.gz"; \
  system(cmd)}' Schmid2018_acc.txt

# Trim remnants of poly(A) tails from 5' ends and Illumina adapters from 3' ends:
for file in Schmid*fastq.gz; do 
  echo $file && 
  cutadapt -j 4 -q 10 -g "T{150};e=0.05" --max-n 3 $file | 
  cutadapt -a AGATCGGAAGAGC -m 20 \
    -o ${file/.fastq.gz/_trim.fq.gz} -; 
done

# Align to the yeast genome:
for file in Schmid*trim.fq.gz; do 
  echo $file && 
  STAR --genomeDir saccer3_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/trim.fq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq AGATCGGAAGAG \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARtmp *out *tab

# Sort BAM files and filter by MAPQ values:
for file in Schmid*bam; do 
  echo $file && 
  samtools view -hq 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Merge "noPAP_total" BAM files:
samtools merge Schmid2018_merged_mapq.bam Schmid*noPAP_total*bam

# Make Bedgraph files WITH strand switch:
for str in "+" "-"; do 
  [ "$str" = "+" ] && 
  n="rev" || n="fw"; 
  for file in Schmid*mapq.bam; do 
    sample=${file/_mapq.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 -strand $str | 
    sort -k1,1 -k2,2n > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files for the same sample:
f_str="fw"; r_str="rev"; ext=".bg"; 
for file1 in *${f_str}${ext}; 
  do file2=${file1/${f_str}/${r_str}} && 
  outfile=${file1/${f_str}${ext}/fw_rev.bedgraph.gz} && 
  echo $file1 "+" $file2 "=" $outfile && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $file2 | 
  cat $file1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $outfile; 
done
