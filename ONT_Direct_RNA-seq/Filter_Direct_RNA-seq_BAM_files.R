library(GenomicAlignments)
library(rtracklayer)
library(tidyverse)
library(devtools)

scripts <- c("skip_duplicated_reads", "filter_ga_by_terminal_subalignments", "write_grl_as_bed12")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}

bamdir <- "."
bamfiles <- list.files(bamdir, pattern = "bam$")

data <- file.path(bamdir, bamfiles) %>% lapply(function(x) { readGAlignments(x) %>% sortSeqlevels() %>% sort() })
names(data) <- bamfiles %>% str_replace("_mapq.bam$", "")

# Filter out chimeric alignments:
data2 <- lapply(data, skip_duplicated_reads)

# Filter out unrealistic alignments:
data3 <- lapply(data2, filter_ga_by_terminal_subalignments)

# Convert to GRangesList and save as BED12:
data4 <- lapply(data3, grglist)
mapply(write_grl_as_bed12, data4, paste0(names(data4), ".bed"))

