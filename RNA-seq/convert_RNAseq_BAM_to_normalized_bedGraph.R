library(GenomicAlignments)
library(tidyverse)
library(devtools)

prefix <- "https://github.com/Maxim-Ivanov/Utility_functions/blob/main/"
paste0(prefix, "convert_GAlignments_to_coverage.R?raw=TRUE") %>% devtools::source_url()
paste0(prefix, "save_GRanges_as_bedGraph.R?raw=TRUE") %>% devtools::source_url()

bamdir <- "."
bamfiles <- list.files(pattern = "bam$")

RNA-seq_BAM_to_normalized_bedGraph <- function(bamfile, mode = "PE", stranded = TRUE, switch_strand = TRUE, normalize = TRUE, norm_to = 1e06) {
  stopifnot(mode %in% c("PE", "SE"))
  if (mode == "PE") {
    ga <- readGAlignmentPairs(bamfile)
  } else {
    ga <- readGAlignments(bamfile)
  }
  cov <- convert_GAlignments_to_coverage(ga, merge.strands = !stranded, flip.strands = switch_strands, normalize = normalize, norm_to = norm_to)
  return(cov)
}

# For stranded PE data:
cov_list <- lapply(bamfiles, RNA-seq_BAM_to_normalized_bedGraph)

# For stranded SE data (e.g. Kohnen et al., 2016 - PMID 27923878):
cov_list <- lapply(bamfiles, RNA-seq_BAM_to_normalized_bedGraph, mode = "SE")

# Save as bedGraph files:
names(cov_list) <- bamfiles %>% basename() %>% str_replace(".bam")
mapply(save_GRanges_as_bedGraph, col_list, paste0(names(cov_list), ".bedgraph.gz"))
