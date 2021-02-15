library(GenomicAlignments)
library(tidyverse)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(BSgenome.Athaliana.TAIR.TAIR9)
bs <- BSgenome.Athaliana.TAIR.TAIR9
seqlevels(bs) <- seqlevels(txdb)

library(devtools)
scripts <- c("convert_GAlignments_to_coverage", "save_GRanges_as_bedGraph", "load_NETSeq_BAM", "find_splice_sites")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}


# Find all exons in TAIR10 and Araport11:
# (download Araport11 annotation from https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz)
ebg <- exonsBy(txdb, by = "gene")
txdb2 <- makeTxDbFromGFF("Araport11_GFF3_genes_transposons.201606.gff.gz")
ebg2 <- exonsBy(txdb2, by = "gene")
seqinfo(ebg2, new2old = c(1:5, 7, 6)) <- seqinfo(txdb)

# Find the union of TAIR10 and Araport11 splice sites:
ss1 <- find_splice_sites(ebg)
ss1_5p <- ss1[[1]]
ss1_3p <- ss1[[2]]
ss2 <- find_splice_sites(ebg2)
ss2_5p <- ss2[[1]]
ss2_3p <- ss2[[2]]
all_ss_5p <- sort(reduce(c(ss1_5p, ss2_5p)))
all_ss_3p <- sort(reduce(c(ss1_3p, ss2_3p)))

bamdir <- "."
bamfiles <- list.files(bamdir, pattern = "mapq.bam$")

# Load and post-process all BAM files, save Bedgraphs into the working directory:
df <- load_NETSeq_BAM(bamfiles, bamdir, ss_5p = all_ss_5p, ss_3p = all_ss_3p, bs, trim_names = "_sorted_dedup_clean_mapq")

# Save stats on split, SS and misprimed reads:
write.table(df, "df1.txt", sep = "\t")

