# This script loads BAM files returned by the Mammalian_NET-seq_pipeline.sh script;
# mNET-seq reads are filtered by various criteria and saved as Bedgraph and RDS files to the current working directory;
# The RDS files have to be loaded by readRDS();

library(rtracklayer)
library(GenomicAlignments)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(devtools)

scripts <- c("convert_GAlignments_to_coverage", "find_splice_sites", "load_NETSeq_BAM")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}

# Find all annotated splice sites (to skip NET-seq reads originating from splicing intermediates):
human_ebg <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")
human_ss <- find_splice_sites(human_ebg)

mouse_ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")
mouse_ss <- find_splice_sites(mouse_ebg)

# Find BAM files returned by Mammalian_NET-seq_pipeline.sh:
human_bf_1 <- list.files(pattern = "^Mayer2015.*bam$")
human_bf_2 <- list.files(pattern = "^Nojima2015.*bam$")
human_bf_3 <- list.files(pattern = "^Nojima2018_HeLa.*bam$")
human_bamfiles <- c(human_bf_1, human_bf_2, human_bf_3)

mouse_bamfiles <- list.files(pattern = "^Nojima2018_TAP.*bam$")

# Run the pipeline (silently saves RDS and Bedgraph files to the current working directory, returns data frame with statistics):
load_NETSeq_BAM(human_bamfiles, ss_5p = human_ss[[1]], ss_3p = human_ss[[2]])
load_NETSeq_BAM(mouse_bamfiles, ss_5p = mouse_ss[[1]], ss_3p = mouse_ss[[2]])



