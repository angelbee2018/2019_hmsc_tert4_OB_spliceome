# trimming and editing GTF files to be compatible with TranscriptCoder
# reconstructed_transript_gtf_list_path must contain THE WHOLE PATH to all the .gtf files you want to process.

reference_gtf_path <- "Z:/hg38_ensembl_reference/gtf/Homo_sapiens.GRCh38.98.gtf"

reconstructed_transcript_gtf_list_path <- "Z:/PGNEXUS_kassem_MSC/Kassem_OB_fastqc/analysis_strawberry/results_assemblyonly/merged/reconstructed_transcript_list.txt"

library(tidyverse)
library(purrr)
library(data.table)
library(dplyr)
library(rtracklayer)
library(biomaRt)

# process the reference GTF ################

reference_gtf <- rtracklayer::import(reference_gtf_path)

tibble_reference_gtf_path_startstop_ENST <- reference_gtf %>% as_tibble %>% .[, c("type", "transcript_id")] %>% .[.$type == "start_codon" | .$type == "stop_codon", ]

df_gtf_start.stop.codon._ENST_cast <- dcast(tibble_reference_gtf_path_startstop_ENST, formula = transcript_id ~ type)

vector_ENST_with_start.and.stop <- df_gtf_start.stop.codon._ENST_cast[df_gtf_start.stop.codon._ENST_cast$start_codon == 1 & df_gtf_start.stop.codon._ENST_cast$stop_codon >= 1, "transcript_id"]

reference_gtf_filtered <- reference_gtf[which(!is.na(match(reference_gtf$transcript_id, vector_ENST_with_start.and.stop)))]

reference_gtf_filtered_2 <- reference_gtf_filtered[which(reference_gtf_filtered$type == "exon" | reference_gtf_filtered$type == "start_codon" | reference_gtf_filtered$type == "stop_codon")]

reference_gtf_filtered_2@elementMetadata@listData <- reference_gtf_filtered_2@elementMetadata@listData[c("type", "transcript_id")]

rtracklayer::export(reference_gtf_filtered_2, con = paste(reference_gtf_path, ".transcriptcoder_ready", sep = ""), format = "gtf")

# list-ify the reference GTF

# list_reference_gtf <- purrr::map(.x = 1:length(reference_gtf) %>% array_tree %>% flatten, .f = ~reference_gtf[.x])

reference_gtf_transcriptcoder.ready <- reference_gtf %>% as_tibble %>%  .[.$type == "exon" | .$type == "start_codon" | .$type == "stop_codon", ] %>% 
  .[, c("seqnames", "type", "start", "end", "strand", "transcript_id")] %>%
  setNames(c("seqid", "type", "start", "end", "strand", "transcript_id"))


# process the reconstructed GTF files ################

reconstructed_transript_gtf_list <- read.delim(reconstructed_transript_gtf_list_path, stringsAsFactors = FALSE, header =  FALSE) %>% unlist

reconstructed_transript_gtf_list <- reconstructed_transript_gtf_list %>% gsub(., pattern = "/media/Ubuntu/sharedfolder/", replacement = "Z:/")

list_reconstructed_gtf_files <- purrr::map(.x = reconstructed_transript_gtf_list, .f = ~rtracklayer::import(.x))

list_reconstructed_gtf_files_transcriptcoder.ready <- purrr::map(.x = list_reconstructed_gtf_files, ~as_tibble(.x) %>% 
                                                                   .[.$type == "exon", ] %>% 
                                                                   .[, c("seqnames", "type", "start", "end", "strand", "transcript_id")] %>%
                                                                   setNames(c("seqid", "type", "start", "end", "strand", "transcript_id")))

purrr::map2(.x = list_reconstructed_gtf_files_transcriptcoder.ready, .y = paste(reconstructed_transript_gtf_list, ".transcriptcoder_ready", sep = "") %>% array_tree %>% flatten, 
           .f = ~rtracklayer::export(.x, con = .y %>% paste, format = "gtf"))

######################################################