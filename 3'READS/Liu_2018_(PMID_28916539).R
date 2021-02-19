# Input: BAM files produced by Liu_2018_(PMID_28916539).sh;
# Output: Bedgraph files showing strand-specific sequencing coverage of the first aligned bases of PASS ("PAS supporting") reads;
# PASS reads have to contain at least two consecutive T's immediately preceding the 5' aligned part;
# Observe that GenomicAlignments always converts SEQ to the forward (Watson) strand;
# Therefore, 5'-terminal TT of reads aligned to the reverse (Crick) strand appear as 3'-terminal "AA" in SEQ;

library(GenomicAlignments)
library(tidyverse)
library(devtools)

scripts <- c("save_GRanges_as_bedGraph", "merge_and_normalize_GRanges", "convert_GRanges_to_coverage")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}

process_3pReads_BAM <- function(bamfile, filter_pass = TRUE) {
  message(bamfile); flush.console()
  bam <- readGAlignments(bamfile, use.names = TRUE, param = ScanBamParam(what = "seq"))
  message("\t", length(bam), " input reads;"); flush.console()
  if (isTRUE(filter_pass)) {
    gr <- pmapToAlignments(granges(bam), bam) %>% unname()
    seq <- bam %>% as.data.frame() %>% .$seq %>% as("DNAStringSet")
    up <- subseq(seq, start = 1, end = start(gr) - 1)
    down <- subseq(seq, start = end(gr) + 1, end = width(seq))
    lgl_1 <- grepl("TT$", up)
    lgl_2 <- grepl("^AA", down)
    pass <- ifelse(strand(bam) == "+", lgl_1, lgl_2)
    message("\t", sum(pass), " PASS reads (", round(mean(pass) * 100), "%);")
    bam <- bam[pass]
  }
  out_gr <- granges(bam) %>% resize(1, "start")
  strand(out_gr) <- ifelse(strand(out_gr) == "+", "-", "+")
  out_cov <- convert_GRanges_to_coverage(out_gr)
  return(out_cov)
}

# Load BAM files, filter PASS reads, trim to first aligned bases, convert to strand-specific coverage:
bamfiles <- list.files(pattern = "^Liu2018.*dedup.bam$")
grlist <- lapply(bamfiles, process_3pReads_BAM)
names(grlist) <- bamfiles %>% str_replace("_mapq_dedup.bam", "")

# Save individual replicates as Bedgraph:
for (i in seq_along(grlist)) {
  save_GRanges_as_bedGraph(grlist[[i]], paste0(names(grlist)[[i]], ".bedgraph.gz"))
}

# Merge replicates/samples and save as Bedgraph:
grlist[1:2] %>% merge_and_normalize_GRanges(norm = FALSE) %>% save_GRanges_as_bedGraph("Liu2018_WT_Min.bedgraph.gz")
grlist[3:4] %>% merge_and_normalize_GRanges(norm = FALSE) %>% save_GRanges_as_bedGraph("Liu2018_WT_YPD.bedgraph.gz")
grlist %>% merge_and_normalize_GRanges(norm = FALSE) %>% save_GRanges_as_bedGraph("Liu2018_WT_merged.bedgraph.gz")


