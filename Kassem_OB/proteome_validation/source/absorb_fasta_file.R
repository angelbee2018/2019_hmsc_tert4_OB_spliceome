script_description <- "# MERGE FASTA FILE ######
Takes an unlimited number of arguments (all absolute file paths to .fasta files) and absorbs all into the first .fasta file specified.
This means any fasta entries that are subsets of any entry from the first fasta file will be deleted, and a new fasta file will be written, merging all entries together.
Final merged fasta is outputted in the directory R was run from."

# change the number of CPU cores to use here
ncores <- 16
# is this a DNA or AA fasta? options: "DNA" or "AA"
dna_or_aa <- "AA"
# do we want unique fasta sequences only?
unique_fasta_output <- FALSE

# NO NEED TO TOUCH ANYTHING BELOW HERE ###

args <- c("/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/results_database_generation/PSISigma_constitutive_ensembl_3FT.fasta",
          "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/results_database_generation/PSISigma_constitutive_strawberry_3FT.fasta")

#############

# print the arguments received by the R script
cat("Arguments input:", commandArgs(), sep = "\n")
args = 
  commandArgs(trailingOnly = TRUE)
# cat(args)
cat("number of arguments specified:", length(args))

if (length(args) < 2) {
  stop("please specify more than one .fasta file.")
}

# SET ENVIRONMENT ##########
#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")
#BiocManager::install(c("seqinr", "tidyverse", "purrr", "dplyr", "rtracklayer", "data.table", "furrr", "RhpcBLASctl", "optparse", "tictoc"))

library(seqinr)
library(tidyverse)
library(furrr)
library(data.table)

plan(multiprocess)
options(mc.cores = ncores)

library(tictoc)
# start counting execution time of the whole script
tictoc::tic("Overall execution time")

cat("create list of .fasta files\n")
list_tibbles_fasta_files <- future_map(.x = args,
                                       .f = ~seqinr::read.fasta(file = .x, seqtype = dna_or_aa, forceDNAtolower = FALSE, as.string = TRUE) %>% 
                                         unlist %>% tibble::enframe(name = "fasta_header", value = "sequence"),
                                       progress = TRUE)

tibble_first_fasta_file <- list_tibbles_fasta_files[[1]]

tibble_other_fasta_files <- list_tibbles_fasta_files[2:length(list_tibbles_fasta_files)] %>% rbindlist %>% as_tibble
# sorry, we can't have any *'s in the sequences to be absorbed.
vector_elements_with_non_alphabetical <- grep(x = tibble_other_fasta_files$sequence, pattern = "[^a-zA-Z]")
if (vector_elements_with_non_alphabetical %>% length > 0) {
  warning("Non-alphabetical sequence entries detected in one of the fasta files to be absorbed (probably a stop codon, a.k.a \"*\"). Removing...")
  tibble_other_fasta_files <- tibble_other_fasta_files[-vector_elements_with_non_alphabetical, ]
}

cat("filter substrings\n")
## add column to flag if the virtual peptide per row is a substring of any entry in the first fasta file.
# replace asterisk with \\\\\\\\* for grep
vector_substring_or_not <- future_map(.x = tibble_other_fasta_files$sequence, 
                                      .f = function(a1) {
                                        
                                        # cat(a2, "\n")
                                        
                                        # DEBUG ###
                                        # a1 <- tibble_other_fasta_files$sequence %>% .[[5]]
                                        ###########
                                        
                                        grepl(x = tibble_first_fasta_file$sequence, pattern = a1) %>% any == TRUE
                                        
                                      }, .progress = TRUE, .options = future_options(globals = c("vector_virtual_peptide_sequence"))) %>% unlist

## filter
tibble_other_fasta_files <- tibble_other_fasta_files %>% add_column("substring_or_not" = vector_substring_or_not)

# filter out substrings
tibble_other_fasta_files_no_substring <- tibble_other_fasta_files %>% dplyr::filter(substring_or_not == FALSE)

# merge all fastas and export
tibble_combined_fasta <- dplyr::bind_rows(tibble_first_fasta_file, tibble_other_fasta_files_no_substring %>% dplyr::select(-substring_or_not))

if (unique_fasta_output == TRUE) {
  tibble_combined_fasta <- tibble_combined_fasta %>% dplyr::distinct(sequence, .keep_all = TRUE)
}

cat("write the final fasta\n")
write.fasta(sequences = tibble_combined_fasta$sequence %>% array_tree %>% flatten, names = tibble_combined_fasta$fasta_header, file.out = paste(getwd(), "/merged_fasta_file.fasta", sep = ""), open = "w", nbchar = 40, as.string = TRUE)

# finish counting
tictoc::toc()

q()

