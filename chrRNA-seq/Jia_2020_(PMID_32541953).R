# The regular Clontech (Takara) SMARTer kit produces double-stranded cDNA with identical flanks on both sides (see Clontech_SMARTer_oligo_system.pdf);
# Thus, sequencing such cDNA on ONT or PacBio results in long reads with unknown strand orientation. This is a common problem with ONT/PacBio full-length cDNA datasets
# However, Jia et al. 2020 (PMID 32541953) used Clontech SMARTer kit with custom primers which results in different cDNA flanks (see Clontech_SMARTer_oligo_system.pdf);
# Therefore, strand orientation of long reads can be recovered by analyzing the cDNA adapter sequences;
# The code below assigns strand info to long reads produced by the Jia_2020_(PMID_32541953).sh script;

library(GenomicAlignments)
library(tidyverse)
library(devtools)

scripts <- c("skip_duplicated_reads", "filter_ga_by_terminal_subalignments", "write_grl_as_bed12")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}

adapter <- "CTGTAGGCACCATCAAT"
adapter_rc <- "ATTGATGGTGCCTACAG"

process_stranded_Clontech_BAM <- function(bamfile) {
  message(bamfile); flush.console()
  bam <- readGAlignments(bamfile, use.names = TRUE, param = ScanBamParam(what = "seq")) # SEQ is converted to Watson strand
  # (as a result, both forward and reverse reads from genes encoded on the Watson strand always appear as if they were having a downstream adapter)
  # (similarly, genes encoded on the Crick strand always appear as having an upstream adapter)
  message("\t", length(bam), " input reads;"); flush.console()
  # Skip chimeric alignments:
  bam <- skip_duplicated_reads(bam)
  message("\t", length(bam), " non-chimeric reads;"); flush.console()
  # Filter by terminal subalignments:
  bam <- filter_ga_by_terminal_subalignments(bam)
  message("\t", length(bam), " reads after filtering by terminal subalignments;"); flush.console()
  # Detect 3' adapters:
  gr <- pmapToAlignments(granges(bam), bam) %>% unname()
  gr2 <- resize(gr, ifelse(width(gr) >= 80, width(gr) - 40, width(gr)), fix = "center")
  seq <- bam %>% as.data.frame() %>% .$seq %>% as("DNAStringSet")
  up <- subseq(seq, start = 1, end = start(gr2))
  down <- subseq(seq, start = end(gr2), end = width(seq))
  m1 <- elementNROWS(vmatchPattern(adapter_rc, up, max.mismatch = 4)) == 1
  m2 <- elementNROWS(vmatchPattern(adapter, down, max.mismatch = 4)) == 1
  rev <- m1 & !m2
  fw <- m2 & !m1
  message("\t", sum(fw), " reads on the Watson strand;"); flush.console()
  message("\t", sum(rev), " reads on the Crick strand;"); flush.console()
  # Extract stranded reads:
  bam_rev <- bam[rev]
  strand(bam_rev) <- "-"
  bam_fw <- bam[fw]
  strand(bam_fw) <- "+"
  bam2 <- c(bam_rev, bam_fw) %>% sort()
  message("\tStrand info recovered for ", round(length(bam2) / length(bam) * 100, 1), "% reads;"); flush.console()
  return(grglist(bam2))
}

bamfiles <- list.files(pattern = "Jia2020.*bam$")
data <- lapply(bamfiles, process_stranded_Clontech_BAM)
names(data) <- bamfiles %>% basename() %>% str_replace(".bam", "")

# Save GRangesList as BED12:
mapply(write_grl_as_bed12, data, paste0(names(data), ".bed"))

