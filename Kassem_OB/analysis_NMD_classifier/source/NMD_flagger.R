script_description <- "# NMD FLAGGER ######
Accepts an input of a reconstructed/assembled transcriptome annotation (e.g. strawberry, stringtie) and flags potential candidate NMD transcripts using the 51nt rule as described in <<Hsu, M.K., Lin, H.Y. and Chen, F.C., 2017. NMD Classifier: A reliable and systematic classification tool for nonsense-mediated decay events. PloS one, 12(4)>>.

Behaviour: loads reconstructed transcriptome, then subsets each transcript. For each transcript, look-up the reference FASTA sequence and 3FT to find ALL start codons that occur more than 54nt upstream of the last exon-exon junction. Then look for any in-frame stop codons downstream. Transcript is an NMD candidate when the FIRST nucleotide of the stop codon is within 51 - 50nt of the last exon junction (the stop codon is fully encapsulated in the 51nt window).

Option is available to use already-annotated start codon to tell the program which translation frame is the correct one.

NMD candidature is flagged in the \"PTC\" entry of the GTF. First/last exon is flagged in a separate entry. Poison exons are defined on a per-transcript basis, with candidature flagged as a separate entry with the \"type\" attribute == \"poison_exon\"

Recommended system requirements: 6 threads/64GB memory. Minimum system requirements: 1 thread/16GB memory"

source("/mnt/LTS/tools/angel_suite/source/main_source.R")

# print the arguments received by the R script
cat("Arguments input:", commandArgs(), sep = "\n")
args = 
  commandArgs(trailingOnly = TRUE)
cat(args)
cat("number of arguments specified:", length(args))

# SET ENVIRONMENT ##########
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(c("seqinr", "tidyverse", "purrr", "dplyr", "rtracklayer", "data.table", "furrr", "RhpcBLASctl", "optparse", "tictoc"))

library(seqinr)
library(tidyverse)
library(purrr)
library(furrr)
library(dplyr)
library(rtracklayer)
library(data.table)
library(optparse)

library(tictoc)
# start counting execution time of the whole script
tictoc::tic("Overall execution time")

# manage arguments
list_input_arg_info = list(
  "1" = make_option(c("-R", "--reconstructed_transcript_gtf"), type = "character", default = NULL, 
                    help = "Compulsory. path to the actual reconstructed GTF file (e.g. from Cufflinks, Strawberry) that you want checked for NMD candidates. NOT THE CONTAINING DIRECTORY.", metavar = "character"),
  "2" = make_option(c("-F", "--reference_genome_fasta_dir"), type = "character", default = NULL, 
                    help = "Compulsory. path to the directory containing the genome FASTA files. Ideally from Ensembl... you need separate files by chromosomes, NOT the primary assembly. 
              FORMATTING IMPORTANT!!!! MAKE SURE THE REF. GENOME FASTA FILES ARE IN THE FORMAT: <_anything_><chr>.fa e.g. \"Homo_sapiens.GRCh38.dna.chromosome.MT.fa\" OR \"chr16.fa\" OR \"Y.fa\". What will not work: anything which does not have .fa extension e.g. \"chr16.fasta\", anything between the chromosome number and the .fa extension e.g. \"chromosome1.ensembl.fa\"", metavar = "character"),
  "3" = make_option(c("-D", "--output_dir"), type = "character", default = NULL, 
                    help = "Compulsory. output file directory. where do you want to save the newly annotated reconstructed GTF? IMPORTANT: MUST BE A FULL DIRECTORY AND NOT A FILE PATH. e.g. correct: ~/outputdir/ correct: ~/outputdir incorrect: /outputdir/annotated_sample.gtf", metavar = "character"),
  "4" = make_option(c("-O", "--output_name"), type = "character", default = NULL, 
                    help = "Compulsory. output file name, to be saved in the output directory a.k.a. what do you want to save the newly annotated reconstructed GTF as? IMPORTANT: MUST BE A STRING WITHOUT THE .GTF EXTENSION AND NOT A DIRECTORY. THE .GTF EXTENSION WILL AUTOMATICALLY BE ADDED FOR THE OUTPUT FILE. e.g. correct: annotated_sample incorrect: annotated_sample.gtf incorrect: annotated_sample/", metavar = "character"),
  "5" = make_option(c("-C", "--ncores"), type = "integer", default = 0, 
                    help = "Optional. Number of cores to use. possible inputs: numbers 1 to any integer. By default, uses all cores (ncores = 0).", metavar = "integer"),
  "6" = make_option(c("-S", "--use_start_codon"), type = "character", default = "YES", 
                    help = "Optional but you should really choose the right option. This option tells the program whether or not there are start codons provided for each transcript (where available). It's useful if you want to re-annotate e.g. the reference Ensembl GTF for the presence of PTCs, so that the program will not consider all 3 translation frames but only consider the frame containing the annotated start codon. If for any reason a transcript doesn't have an associated start codon, the program will revert to considering all 3 frames. By default, start codons are used where applicable.", metavar = "logical"),
  "7" = make_option(c("-H", "--chrmode"), type = "integer", default = 0, 
                    help = "Optional. Specifies which chromosomes to do: select what chromosomes you want translated. possible inputs: numbers 0-2. 0 (default): nuclear chromosomes only i,e, 1:22, X & Y. 1: nuclear + mitochondrial i.e. 1:22, X & Y, M. 2: everything including haplotype/fusion chromosomes etc... this is possible provided the chromosome names.", metavar = "integer"),
  "8" = make_option(c("-N", "--nonchrname"), type = "character", default = NULL, 
                    help = "Compulsory only if you have specified \"--chrmode 2\". nonchromosomal file name. if you want to consider haplotypes, please specify what the reference genome FASTA file for it is called or the script won't know. This single FASTA file must contain all the haplotype information. The program won't try to search for a second file. In ensembl, this file is called something like \"Homo_sapiens.GRCh38.dna.nonchromosomal.fa\" or more generally, \"*nonchromosomal.fa\". So if you want to use that file, then for this option, you would specify \"--nonchrname nonchromosomal\".", metavar = "character"),
  "9" = make_option(c("-W", "--checking_window_size"), type = "double", default = 51,
                    help = "ADVANCED. Use this parameter to change the nucleotide window from the default of 51. 
                    DO NOT CHANGE THIS UNLESS YOU REALLY WANT TO CHECK USING A DIFFERENT WINDOW SIZE. SUBOPTIMAL RESULTS MAY BE GENERATED AS A RESULT.", metavar = "double"),
  "10" = make_option(c("-E", "--min_exons_per_transcript"), type = "integer", default = 3,
                    help = "ADVANCED. Use this parameter to change the minimum number of exons that a transcript must have in order to be flagged for NMD. As a general rule of thumb, the false positive rate:false negative rate ratio will decrease as you increase this parameter. Default is 3, meaning 2-exon transcripts and below will not be considered for NMD. 
                    DO NOT CHANGE THIS UNLESS YOU REALLY WANT TO CHECK USING A DIFFERENT WINDOW SIZE. SUBOPTIMAL RESULTS MAY BE GENERATED AS A RESULT.", metavar = "double"),
  "11" = make_option(c("-V", "--save_workspace_when_done"), type = "character", default = FALSE,
                     help = "Turn this on if you want to save the R workspace in the same name as the --output_name. YES: saves at the end. DEBUG: saves at each critical step. NO: doesn't save.", metavar = "character")
)

# "7" = make_option(c("-R", "--return_frames"), type = "character", default = "3FT_aggregate", 
#                   help = "Optional. Changes how NMD prediction is reported. 3FT_aggregate: requires that all THREE frames have a stop codon in the required window - extremely conservative. per_frame: reports on whether there is a PTC per reading frame (relative to the chromosome) - chrframe 1: 1st nt, chrframe 2: 2nd nt, chrframe 3: 3rd nt (from the 5' end of the chromosome for + strand and from the 3' end for - strand", metavar = "logical"),

input_arg_info <- OptionParser(option_list = list_input_arg_info, description = script_description)
input_args <- input_arg_info %>% parse_args

# check if the input arguments are O.K
if ((list(input_args$reconstructed_transcript_gtf, input_args$reference_genome_fasta_dir, input_args$output_name) %>% lapply(is.null) %>% unlist %>% any == TRUE) | 
    (input_args$chrmode == 2 & is.null(input_args$nonchrname) == TRUE)) {
  
  print_help(input_arg_info)
  
  stop("Make sure you entered the arguments correctly", call. = FALSE)
  
}

# DEBUG #######

reconstructed_gtf_path <- "/mnt/LTS/reference_data/hg38_ensembl_reference/gtf/Homo_sapiens.GRCh38.98.gtf"
reference_genome_fasta_dir <- "/mnt/LTS/reference_data/hg38_ensembl_reference/raw_genome_fasta/dna_by_chr/"
output_dir <- "/mnt/LTS/projects/2019_hmsc_spliceome/Kassem_OB/analysis_NMD_classifier/results/"
output_name <- "Homo_sapiens.GRCh38.98_NMD_PTC_E4"
window_size <- 51
min_exons_per_transcript <- 3
use_start_codon <- "YES"
number_of_workers <- 32

###############

reconstructed_gtf_path <- input_args$reconstructed_transcript_gtf
reference_genome_fasta_dir <- input_args$reference_genome_fasta_dir
output_dir <- input_args$output_dir
output_name <- input_args$output_name
window_size <- input_args$checking_window_size
min_exons_per_transcript <- input_args$min_exons_per_transcript
use_start_codon <- input_args$use_start_codon
save_workspace_when_done <- input_args$save_workspace_when_done
# return_frames <- input_args$return_frames

cat("reconstructed_gtf_path:", reconstructed_gtf_path, "\n")
cat("reference_genome_fasta_dir:", reference_genome_fasta_dir, "\n")
cat("output_name:", output_name, "\n")
cat("window_size:", window_size, "\n")
cat("min_exons_per_transcript:", min_exons_per_transcript, "\n")
# cat("return_frames:", return_frames, "\n")

# if(!dir.exists(output_name) ) {
#   dir.create(output_name, recursive = TRUE)}

# manage parrallellisation rrlllRll

if (input_args$ncores != 0) {
  number_of_workers <- input_args$ncores
} 

cat(number_of_workers, "cores will be used\n")
options(future.globals.maxSize = 30000000000, mc.cores = number_of_workers, future.fork.enable = TRUE)
future::plan(multicore)

# DEFINE FUNCTIONS ##########################

# FUNCTION TO 3 FRAME TRANSLATE ONE LIST CONTAINING NUCLEOTIDE SEQUENCE AND STRAND

nt.sequence_strand_threeframetranslate <- function(vector_forward_nucleotides, strand, frame_adjust = 0) {
  
  if (strand == "+") {
    
    translation_result <- list("translation_frame_0" = seqinr::translate(vector_forward_nucleotides, frame = (0 + frame_adjust)%%3, sens = "F"),
                               "translation_frame_1" = seqinr::translate(vector_forward_nucleotides, frame = (1 + frame_adjust)%%3, sens = "F"),
                               "translation_frame_2" = seqinr::translate(vector_forward_nucleotides, frame = (2 + frame_adjust)%%3, sens = "F"))
    
  } else if (strand == "-") {
    
    translation_result <- list("translation_frame_0" = seqinr::translate(vector_forward_nucleotides, frame = (0 + frame_adjust)%%3, sens = "R"),
                               "translation_frame_1" = seqinr::translate(vector_forward_nucleotides, frame = (1 + frame_adjust)%%3, sens = "R"),
                               "translation_frame_2" = seqinr::translate(vector_forward_nucleotides, frame = (2 + frame_adjust)%%3, sens = "R"))
    
  }
  
  return(translation_result)
  
}

# END nt.sequence_strand_threeframetranslate()

# FUNCTION to test if there is a valid ORF or not. 
test_for_any_valid_ORF <- function(vector_AA_sequence) {
  
  # see if the reverse of the sequence has a methionine before any stop codons are encountered OR if the WHOLE sequence has no stop codons at all.
  validity_test <- stringr::str_detect(vector_AA_sequence %>% rev %>% paste(collapse = ""), "^[^\\*]+M|^[^\\*]+$")
  
  return(validity_test)
  
}

# END test_for_any_valid_ORF()

# FUNCTION TO LABEL THE FIRST AND LAST EXON OF A SUBSET GTF (DESCRIBING A SINGLE TRANSCRIPT) WITH "first_exon" or "last_exon"

label_first_last_exon <- function(tibble_gtf_subset_by_transcript_id) {
  
  max_exon_number <- max(tibble_gtf_subset_by_transcript_id$exon_number)
  
  # account for transcripts which only have one exon
  if (max_exon_number == 1) {
    
    tibble_gtf_subset_by_transcript_id[tibble_gtf_subset_by_transcript_id$exon_number == max_exon_number, "first_or_last_exon"] <- "only_one_exon"
    
  } else if (exon_order == "increasing" & tibble_gtf_subset_by_transcript_id$strand %>% unique == "-") {
    
    tibble_gtf_subset_by_transcript_id[tibble_gtf_subset_by_transcript_id$exon_number == max_exon_number, "first_or_last_exon"] <- "first_exon"
    
    tibble_gtf_subset_by_transcript_id[tibble_gtf_subset_by_transcript_id$exon_number == 1, "first_or_last_exon"] <- "last_exon"
    
  } else {
    
    tibble_gtf_subset_by_transcript_id[tibble_gtf_subset_by_transcript_id$exon_number == max_exon_number, "first_or_last_exon"] <- "last_exon"
    
    tibble_gtf_subset_by_transcript_id[tibble_gtf_subset_by_transcript_id$exon_number == 1, "first_or_last_exon"] <- "first_exon"
    
  }
  
  return(tibble_gtf_subset_by_transcript_id)
  
}

# END label_first_last_exon() ######

# MAIN FUNCTION TO LABEL THE FIRST/LAST EXONS OF GTF

add_first.last.exon_info <- function(input_gtf) {
  
  # add in the first_or_last_exon column
  input_gtf_2 <- input_gtf %>% add_column("first_or_last_exon" = "NA")
  
  # remove rows which have an NA exon number
  input_gtf_2 <- input_gtf_2[!is.na(input_gtf_2$exon_number), ]
  
  # remove rows which don't have strand info
  input_gtf_2 <- input_gtf_2 %>% dplyr::filter(input_gtf_2$strand != "+" | input_gtf_2$strand != "-")
  
  # get all the unique transcript IDs for looping
  list_unique_transcript.ids <- input_gtf_2$transcript_id %>% unique %>% array_tree
  
  output_gtf0 <- furrr::future_map(.x = list_unique_transcript.ids, .f = ~label_first_last_exon(input_gtf_2[input_gtf_2$transcript_id == .x, ]), .progress = TRUE) %>% rbindlist %>% as_tibble
  
  output_gtf <- output_gtf0
  
  return(output_gtf)
  
}

# BEGIN EXECUTION ###########################

cat("import reconstructed transcriptome GTF\n")
tibble_reconstructed_gtf <- rtracklayer::import(reconstructed_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)

tibble_reconstructed_gtf_no_asterisk <- tibble_reconstructed_gtf[tibble_reconstructed_gtf$strand %in% c("+", "-"), ]

vector_ref_genome_paths_by_chr <- paste(reference_genome_fasta_dir, list.files(reference_genome_fasta_dir)[list.files(reference_genome_fasta_dir) %>% grep(., pattern = ".*.fa$")], sep = "")

# automatically detect if exons are always numbered in increasing order regardless of strand (common for recon. transcripts)
## sample the first transcript on the negative strand with more than 1 exon
temp_number <- 1

first_transcript_id <- tibble_reconstructed_gtf_no_asterisk$transcript_id %>% unique %>% na.omit %>% .[temp_number]

while (tibble_reconstructed_gtf_no_asterisk[tibble_reconstructed_gtf_no_asterisk$transcript_id == first_transcript_id, "exon_number"] %>% nrow == 1 | 
       tibble_reconstructed_gtf_no_asterisk[tibble_reconstructed_gtf_no_asterisk$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "+" | 
       tibble_reconstructed_gtf_no_asterisk[tibble_reconstructed_gtf_no_asterisk$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "." | 
       tibble_reconstructed_gtf_no_asterisk[tibble_reconstructed_gtf_no_asterisk$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "*") {
  
  temp_number <- temp_number + 1
  
  first_transcript_id <- tibble_reconstructed_gtf_no_asterisk$transcript_id %>% unique %>% na.omit %>% .[temp_number]
  
}

# the test condition
max_test <- tibble_reconstructed_gtf_no_asterisk[tibble_reconstructed_gtf_no_asterisk$transcript_id == first_transcript_id, "exon_number"] %>% unlist %>% na.omit %>% max
max_exon_start_test <- tibble_reconstructed_gtf_no_asterisk[tibble_reconstructed_gtf_no_asterisk$transcript_id == first_transcript_id & tibble_reconstructed_gtf_no_asterisk$exon_number == max_test, "start"] %>% na.omit %>% paste
min_exon_start_test <- tibble_reconstructed_gtf_no_asterisk[tibble_reconstructed_gtf_no_asterisk$transcript_id == first_transcript_id & tibble_reconstructed_gtf_no_asterisk$exon_number == 1, "start"] %>% na.omit %>% paste

exon_order <- NULL

# if the exon 1 comes before the max exon, then the exon order is always increasing
if (min_exon_start_test < max_exon_start_test) {
  
  exon_order <- "increasing"
  
  # if the exon 1 cones after the max exon, then the exon order is stranded.
} else if (min_exon_start_test > max_exon_start_test) {
  
  exon_order <- "stranded"
  
}

# specify the chromosomes to be run, according to user option --chrmode
if (input_args$chrmode == 1) {
  
  # all nuclear chromosomes + mitochondria
  chr_to_run <- c(1:22, "X", "Y", "MT")
  
} else if (input_args$chrmode == 2) {
  
  chr_to_run <- c(1:22, "X", "Y", "MT", input_args$nonchrname)
  
} else {
  
  # if the user put in a stupid number then we'll assume they just want all the nuclear chromosomes.
  chr_to_run <- c(1:22, "X", "Y")
  
}

# begin looping thru each chromosome #####
for (chr in chr_to_run) {
  
  # start counting
  tictoc::tic(paste("Chromosome", chr))
  
  cat("now running chromosome ...", chr, "\n")
  
  cat("get positions of the vector where the path of the ref. genome fasta\n")
  vector_ref_genome_paths_by_chr_position <- grep(x = vector_ref_genome_paths_by_chr, pattern = paste("(\\D|^)", chr, ".fa$", sep = ""))
  
  if (length(vector_ref_genome_paths_by_chr_position) != 1) {
    
    stop("Something is wrong with the contents of the fasta file directory. Please check that it's structured in the desired format.")
    
  }
  
  cat("temporary allocation to ref genome fasta list\n")
  reference_genome_fasta_chr_temp <- seqinr::read.fasta(file = paste(vector_ref_genome_paths_by_chr[vector_ref_genome_paths_by_chr_position]), forceDNAtolower = FALSE)
  
  tibble_reconstructed_gtf_no_asterisk_tempchr <- tibble_reconstructed_gtf_no_asterisk[tibble_reconstructed_gtf_no_asterisk$seqnames == chr, ]
  
  cat("subset reconstructed transcriptome GTF by filtering for exon entries only.\n")
  # leave aside the remaining entries for later inclusion
  ## rows which do not describe an exon or start_codon entry
  row_indices_recon.gtf_non.exon.entries <- which(tibble_reconstructed_gtf_no_asterisk_tempchr$type != "exon" & tibble_reconstructed_gtf_no_asterisk_tempchr$type != "start_codon")
  ## rows where the transcript_id is NA
  row_indices_recon.gtf_na.transcript.id.entries <- which(is.na(tibble_reconstructed_gtf_no_asterisk_tempchr$transcript_id))
  
  row_indices_recon.gtf_exon.entries.only <- which(tibble_reconstructed_gtf_no_asterisk_tempchr$type == "exon" | tibble_reconstructed_gtf_no_asterisk_tempchr$type == "start_codon")
  reconstructed_gtf_exon.entries.only <- tibble_reconstructed_gtf_no_asterisk_tempchr[row_indices_recon.gtf_exon.entries.only, ]
  
  # generate list of all 'transcript_id's to loop thru
  list_transcript_ids <- tibble_reconstructed_gtf_no_asterisk_tempchr$transcript_id %>% unique %>% na.omit %>% array_tree
  
  cat("load the exon entries of the recon. GTF table into the list\n")
  list_reconstructed_gtf_exon.entries.only_by_transcript_id <- future_map(.x = list_transcript_ids, .f = ~reconstructed_gtf_exon.entries.only[reconstructed_gtf_exon.entries.only$transcript_id == .x, ], .progress = TRUE) %>% set_names(list_transcript_ids %>% unlist)
  
  # filter out transcript entries which only have one exon
  list_reconstructed_gtf_exon.entries.only_by_transcript_id <- purrr::discard(.x = list_reconstructed_gtf_exon.entries.only_by_transcript_id, .p = ~length(.x$exon_number %>% na.omit) < min_exons_per_transcript)
  
  cat("generate attributes of each transcript:
      1. all genome-relative forward coords of each nucleotide,
      2. coords of the nucleotides flanking last exon-exon junction\n")
  list_transcript_attributes <- future_imap(.x = list_reconstructed_gtf_exon.entries.only_by_transcript_id, 
                                            .f = function(a1, a2) {
                                              
                                              # cat("\nnow processing: ", a2)
                                              
                                              temp0 <- list(
                                                "recon_entries" = a1,
                                                "strand" = a1$strand %>% unique %>% paste,
                                                "forward_coords" = purrr::map2(.x = a1 %>% dplyr::filter(type == "exon") %>% .$start, .y = a1 %>% dplyr::filter(type == "exon") %>% .$end, .f = ~.x:.y) %>% unlist %>% as.numeric %>% sort, 
                                                # since the start is always less than the end in genome relative coords, for the fwd strand last_nt_of_second_last_exon will be the end coord of the second last exon.
                                                # however for rev. strand, last_nt_of_second_last_exon will be the start coord of the second last exon, and the first_nt_of_last_exon will be the end coord of the last exon.
                                                "last_nt_of_second_last_exon" = if (a1$strand %>% unique %>% paste == "+") {
                                                  a1 %>% dplyr::filter(type == "exon" & exon_number == (exon_number %>% as.numeric %>% max - 1)) %>% .$end %>% paste %>% as.numeric
                                                } else if (exon_order == "stranded" & a1$strand %>% unique %>% paste == "-") {
                                                  a1 %>% dplyr::filter(type == "exon" & exon_number == (exon_number %>% as.numeric %>% max - 1)) %>% .$start %>% paste %>% as.numeric
                                                } else if (exon_order == "increasing" & a1$strand %>% unique %>% paste == "-") {
                                                  a1 %>% dplyr::filter(type == "exon" & exon_number == (exon_number %>% as.numeric %>% min + 1)) %>% .$start %>% paste %>% as.numeric
                                                }, 
                                                "first_nt_of_last_exon" = if (a1$strand %>% unique %>% paste == "+") {
                                                  a1 %>% dplyr::filter(type == "exon" & exon_number == (exon_number %>% as.numeric %>% max)) %>% .$start %>% paste %>% as.numeric
                                                } else if (exon_order == "stranded" & a1$strand %>% unique %>% paste == "-") {
                                                  a1 %>% dplyr::filter(type == "exon" & exon_number == (exon_number %>% as.numeric %>% max)) %>% .$end %>% paste %>% as.numeric
                                                } else if (exon_order == "increasing" & a1$strand %>% unique %>% paste == "-") {
                                                  a1 %>% dplyr::filter(type == "exon" & exon_number == (exon_number %>% as.numeric %>% min)) %>% .$end %>% paste %>% as.numeric
                                                } )
                                              
                                              # only splice in the start_codon entry if use_start_codon is specified
                                              if (use_start_codon == "YES") {
                                                
                                                temp0 %>% purrr::splice("start_codon" = a1 %>% dplyr::filter(type == "start_codon")) %>% return
                                                
                                              } else {
                                                
                                                temp0 %>% purrr::splice("start_codon" = a1[0, ]) %>% return
                                                
                                              }
                                              
                                            }, .progress = TRUE
                                            
  )
  
  cat("generate all the forward nucleotides of each transcript\n")
  # also get the TRANSCRIPT RELATIVE position of the junction  (transcript-relative coord will be a half-integer value)
  # also get the transcript-relative position of the start_codon (if it's present)
  list_transcript_fwd_nts <- future_imap(.x = list_transcript_attributes, .f = function(a1, a2) {
    
    # DEBUG ###
    # a1 <- list_transcript_attributes[[1]]
    ###########
    
    cat("now looking up nucleotides of transcript number", which(names(list_transcript_attributes) == a2), "/", length(list_transcript_attributes), "\n")
    
    # if use_start_codon is called and actually exists in the GTF, get the transcript-relative position of the first nt of the start codon according to the matched entry.
    # if start_codon is not available for whatever reason, the first nt of the transcript is taken instead.
    if (a1$start_codon %>% nrow == 0) {
      
      updated_list <- purrr::splice(a1,
                                    "fwd_coordinates_relative_first_nt_of_start_codon" = if (a1$strand == "+") {1} else if (a1$strand == "-") {length(a1$forward_coords)})
      
    } else if (a1$start_codon %>% nrow > 0) {
      
      updated_list <- purrr::splice(a1,
                                    "fwd_coordinates_relative_first_nt_of_start_codon" = if (a1$strand %>% unique %>% paste == "+") {
                                      which(a1$start_codon %>% type_convert %>% dplyr::filter(exon_number == (exon_number %>% as.numeric %>% min)) %>% .$start == a1$forward_coords)
                                    } else if (exon_order == "stranded" & a1$strand %>% unique %>% paste == "-") {
                                      which(a1$start_codon %>% type_convert %>% dplyr::filter(exon_number == (exon_number %>% as.numeric %>% min)) %>% .$end == a1$forward_coords)
                                    } else if (exon_order == "increasing" & a1$strand %>% unique %>% paste == "-") {
                                      which(a1$start_codon %>% type_convert %>% dplyr::filter(exon_number == (exon_number %>% as.numeric %>% max)) %>% .$end == a1$forward_coords)
                                    } )
      
    }
    
    updated_list <- purrr::splice(updated_list,
                                  "forwards_coords_relative_junction_position" = mean(c(
                                    # transcript relative position of the last nt of second last exon
                                    which(a1$forward_coords == a1$last_nt_of_second_last_exon),
                                    # transcript relative position of the first nt of last exon
                                    which(a1$forward_coords == a1$first_nt_of_last_exon))),
                                  "forward_transcript_nucleotides" = reference_genome_fasta_chr_temp[[chr %>% paste]][a1$forward_coords],
                                  "forward_CDS_nucleotides" = if (a1$strand == "+") {
                                    reference_genome_fasta_chr_temp[[chr %>% paste]][a1$forward_coords %>% .[updated_list$fwd_coordinates_relative_first_nt_of_start_codon:length(a1$forward_coords)]]
                                  } else if (a1$strand == "-") {reference_genome_fasta_chr_temp[[chr %>% paste]][a1$forward_coords %>% .[1:updated_list$fwd_coordinates_relative_first_nt_of_start_codon]]
                                  },
                                  "frame_adjust" = if (a1$strand == "+") {
                                    (3-(((updated_list$forward_coords %>% .[updated_list$fwd_coordinates_relative_first_nt_of_start_codon])-1)%%3))%%3
                                  } else if (updated_list$strand == "-") {
                                    (3-((length(reference_genome_fasta_chr_temp[[chr %>% paste]]) - (updated_list$forward_coords %>% .[updated_list$fwd_coordinates_relative_first_nt_of_start_codon]) )%%3))%%3
                                  } )
    
    return(updated_list)
    
  }, .progress = TRUE)
  
  cat("do 3FT based on the forward nucleotides of each transcript and test for whether there is any valid ORF.\n")
  list_transcript_3FT <- future_imap(.x = list_transcript_fwd_nts, .f = function(a1, a2) {
    
    # DEBUG ###
    # a1 <- list_transcript_fwd_nts[[106]]
    # a2 <- 1
    ###########
    
    cat("now checking 3FT of transcript number", which(names(list_transcript_fwd_nts) == a2), "/", length(list_transcript_fwd_nts), "\n")
    
    # for simplicity, we exclude the last exon + 51 nt window from the nt sequence
    # hence get the effective fwd. nucleotide sequence:
    vector_all_forward_CDS_nucleotides <- a1$forward_CDS_nucleotides
    vector_all_forward_transcript_nucleotides <- a1$forward_transcript_nucleotides
    
    # if (weirdly enough) the transcript is actually shorter than 51 nt, then we skip it. Instances of this occurring are very rare indeed.
    # requirement is that the sequence length longer than the window is at least window + 5 because you need at least 5nt in order for 3FT to work succesfully in all frames
    if (((a1$strand == "+") & (floor(a1$forwards_coords_relative_junction_position) - window_size + 1 - a1$fwd_coordinates_relative_first_nt_of_start_codon + 1) < 5) |
        ((a1$strand == "-") & (a1$fwd_coordinates_relative_first_nt_of_start_codon - floor(a1$forwards_coords_relative_junction_position) - window_size + 1 < 5))) {
      
      return(NA)
      
    } else if (a1$strand == "+") {
      
      vector_forward_nucleotides_excl.last.exon <- vector_all_forward_transcript_nucleotides[(a1$fwd_coordinates_relative_first_nt_of_start_codon):(floor(a1$forwards_coords_relative_junction_position) - window_size + 1)]
      
    } else if (a1$strand == "-") {
      
      vector_forward_nucleotides_excl.last.exon <- vector_all_forward_transcript_nucleotides[(ceiling(a1$forwards_coords_relative_junction_position) + window_size - 1):(a1$fwd_coordinates_relative_first_nt_of_start_codon)]
      
    }
    
    vector_all_forward_CDS_nucleotides_excl.last.exon <- extract_common_string(vector_all_forward_CDS_nucleotides %>% paste(collapse = ""), vector_forward_nucleotides_excl.last.exon %>% paste(collapse = "")) %>% strsplit(split = "") %>% unlist
    
    # splice in the 3FT virtual amino acids
    updated_list_temp <- purrr::splice(a1, nt.sequence_strand_threeframetranslate(vector_forward_nucleotides = vector_all_forward_CDS_nucleotides_excl.last.exon, frame_adjust = a1$frame_adjust, strand = a1$strand))
    # if start codon info was available, remove the translation frames 1 and 2.
    if (updated_list_temp$start_codon %>% nrow > 0) {
      
      updated_list_temp[[paste("translation_frame_", (1 - a1$frame_adjust)%%3, sep ="")]] <- updated_list_temp[[paste("translation_frame_", (-a1$frame_adjust)%%3, sep ="")]]
      updated_list_temp[[paste("translation_frame_", (2 - a1$frame_adjust)%%3, sep ="")]] <- updated_list_temp[[paste("translation_frame_", (-a1$frame_adjust)%%3, sep ="")]]
      
    }
    
    if (updated_list_temp$strand == "+") {
      vector_stranded_CDS_nucleotides <- updated_list_temp$forward_CDS_nucleotides
    } else if (updated_list_temp$strand == "-") {
      vector_stranded_CDS_nucleotides <- updated_list_temp$forward_CDS_nucleotides %>% rev %>% seqinr::comp() %>% toupper
    }
    
    vector_stranded_CDS_nucleotides_codonised_frame_0 <- vector_stranded_CDS_nucleotides[(1 + (updated_list_temp$frame_adjust + 0)%%3):length(vector_stranded_CDS_nucleotides)] %>% seqsplitter(x = ., n = 3, keep_remainders = FALSE) %>% purrr::map(~.x %>% paste(collapse = "")) %>% unlist 
    vector_stranded_CDS_nucleotides_codonised_frame_1 <- vector_stranded_CDS_nucleotides[(1 + (updated_list_temp$frame_adjust + 1)%%3):length(vector_stranded_CDS_nucleotides)] %>% seqsplitter(x = ., n = 3, keep_remainders = FALSE) %>% purrr::map(~.x %>% paste(collapse = "")) %>% unlist 
    vector_stranded_CDS_nucleotides_codonised_frame_2 <- vector_stranded_CDS_nucleotides[(1 + (updated_list_temp$frame_adjust + 2)%%3):length(vector_stranded_CDS_nucleotides)] %>% seqsplitter(x = ., n = 3, keep_remainders = FALSE) %>% purrr::map(~.x %>% paste(collapse = "")) %>% unlist 
    
    list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon <- purrr::map(
      .x = numbers_to_intervals(setdiff(1:length(vector_stranded_CDS_nucleotides_codonised_frame_0), which(vector_stranded_CDS_nucleotides_codonised_frame_0 %in% c("TAG", "TAA", "TGA")))) %>% purrr::array_tree(), 
      .f = ~vector_stranded_CDS_nucleotides_codonised_frame_0[.x$start:.x$end]) %>% purrr::set_names(nm = 1:length(.))
    list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon <- purrr::map(
      .x = numbers_to_intervals(setdiff(1:length(vector_stranded_CDS_nucleotides_codonised_frame_1), which(vector_stranded_CDS_nucleotides_codonised_frame_1 %in% c("TAG", "TAA", "TGA")))) %>% purrr::array_tree(), 
      .f = ~vector_stranded_CDS_nucleotides_codonised_frame_1[.x$start:.x$end]) %>% purrr::set_names(nm = 1:length(.))
    list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon <- purrr::map(
      .x = numbers_to_intervals(setdiff(1:length(vector_stranded_CDS_nucleotides_codonised_frame_2), which(vector_stranded_CDS_nucleotides_codonised_frame_2 %in% c("TAG", "TAA", "TGA")))) %>% purrr::array_tree(), 
      .f = ~vector_stranded_CDS_nucleotides_codonised_frame_2[.x$start:.x$end]) %>% purrr::set_names(nm = 1:length(.))
    
    if (updated_list_temp$start_codon %>% nrow > 0) {
    
      flag_start_codon_present <- TRUE
      
      length_longest_coding_sequence_frame0_ATGonly <- NA
      genome_relative_start_coding_sequence_frame0_ATGonly <- NA
      genome_relative_end_coding_sequence_frame0_ATGonly <- NA
      has_ptc_frame0_ATGonly <- NA
      genome_relative_stop_codon_frame0_ATGonly <- NA
      length_longest_coding_sequence_frame1_ATGonly <- NA
      genome_relative_start_coding_sequence_frame1_ATGonly <- NA
      genome_relative_end_coding_sequence_frame1_ATGonly <- NA 
      has_ptc_frame1_ATGonly <- NA
      genome_relative_stop_codon_frame1_ATGonly <- NA
      length_longest_coding_sequence_frame2_ATGonly <- NA
      genome_relative_start_coding_sequence_frame2_ATGonly <- NA
      genome_relative_end_coding_sequence_frame2_ATGonly <- NA
      has_ptc_frame2_ATGonly <- NA
      genome_relative_stop_codon_frame2_ATGonly <- NA
      
      length_longest_coding_sequence_frame0_noncanonical <- length(updated_list_temp$forward_CDS_nucleotides)
      genome_relative_start_coding_sequence_frame0_noncanonical <- updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon] %>% type.convert(as.is = TRUE)
      genome_relative_end_coding_sequence_frame0_noncanonical <- if (updated_list_temp$strand == "+") {
        c(if (any(grepl(x = updated_list_temp$translation_frame_0, pattern = "\\*"))) {updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + (which(!is.na(str_locate(string = updated_list_temp$translation_frame_0, pattern = "\\*") %>% .[, 1]))[1] - 1) * 3]},
          "nostop")[1]
      } else if (updated_list_temp$strand == "-") {
        c(if (any(grepl(x = updated_list_temp$translation_frame_0, pattern = "\\*"))) {updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - (which(!is.na(str_locate(string = updated_list_temp$translation_frame_0, pattern = "\\*") %>% .[, 1]))[1] - 1) * 3]},
          "nostop")[1]
      }
      
      genome_relative_stop_codon_frame0_noncanonical <- if (updated_list_temp$strand == "+") {
        c(if (any(grepl(x = updated_list_temp$translation_frame_0, pattern = "\\*"))) {updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + (which(!is.na(str_locate(string = updated_list_temp$translation_frame_0, pattern = "\\*") %>% .[, 1]))[1] - 1) * 3)] %>% (function(x){return((x+1):(x+3))}) %>% paste(collapse = ",")},
          "nostop")[1]
      } else if (updated_list_temp$strand == "-") {
        c(if (any(grepl(x = updated_list_temp$translation_frame_0, pattern = "\\*"))) {updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - (which(!is.na(str_locate(string = updated_list_temp$translation_frame_0, pattern = "\\*") %>% .[, 1]))[1] - 1) * 3)] %>% (function(x){return((x-1):(x-3))}) %>% paste(collapse = ",")},
          "nostop")[1]
      }
      
      has_ptc_frame0_noncanonical <- any(!is.na(str_locate(string = updated_list_temp$translation_frame_0, pattern = "\\*") %>% .[, 1]))
      genome_relative_stop_codon_frame1_noncanonical <- genome_relative_stop_codon_frame0_noncanonical
      length_longest_coding_sequence_frame1_noncanonical <- length_longest_coding_sequence_frame0_noncanonical
      genome_relative_start_coding_sequence_frame1_noncanonical <- genome_relative_start_coding_sequence_frame0_noncanonical
      genome_relative_end_coding_sequence_frame1_noncanonical <- genome_relative_end_coding_sequence_frame0_noncanonical
      has_ptc_frame1_noncanonical <- has_ptc_frame0_noncanonical
      genome_relative_stop_codon_frame2_noncanonical <- genome_relative_stop_codon_frame0_noncanonical
      length_longest_coding_sequence_frame2_noncanonical <- length_longest_coding_sequence_frame0_noncanonical
      genome_relative_start_coding_sequence_frame2_noncanonical <- genome_relative_start_coding_sequence_frame0_noncanonical
      genome_relative_end_coding_sequence_frame2_noncanonical <- genome_relative_end_coding_sequence_frame0_noncanonical
      has_ptc_frame2_noncanonical <- has_ptc_frame0_noncanonical
      
    } else {
    
      flag_start_codon_present <- FALSE
      
    # find the frame that results in the longest ORF and do it for both ATG and non-ATG starts (apparently >55% of the ORFome)
    list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_ATGonly <- list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon %>% purrr::keep(.p = ~any(grepl(x = .x, pattern = "ATG")) == TRUE)
    
    if (length(list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_ATGonly) > 0) {
      
      # find the length of the longest coding sequence
      length_longest_coding_sequence_frame0_ATGonly <- list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist %>% max
      
      genome_relative_start_coding_sequence_frame0_ATGonly <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + ((updated_list_temp$frame_adjust + 0)%%3)]
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - ((updated_list_temp$frame_adjust + 0)%%3)]
      }
      genome_relative_end_coding_sequence_frame0_ATGonly <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 0)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame0_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3]
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 0)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame0_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3]
      }
      
      # this is where the logical for PTC existence is determined.
      # the location of the end of the ORF does not include the stop codon; it is the last nt before the stop codon
      if (updated_list_temp$strand == "+") {
        
        has_ptc_frame0_ATGonly <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame0_ATGonly) < (floor(updated_list_temp$forwards_coords_relative_junction_position) - window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame0_ATGonly) >= (length(updated_list_temp$forward_coords) - 2)) {
          has_ptc_frame0_ATGonly <- "nostop"
          genome_relative_stop_codon_frame0_ATGonly <- NA
        } else {
          genome_relative_stop_codon_frame0_ATGonly <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 0)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame0_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x+1):(x+3))})] %>% paste(collapse = ",")
        }
        
      } else if (updated_list_temp$strand == "-") {
        
        has_ptc_frame0_ATGonly <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame0_ATGonly) > (ceiling(updated_list_temp$forwards_coords_relative_junction_position) + window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame0_ATGonly) <= 3) {
          has_ptc_frame0_ATGonly <- "nostop"
          genome_relative_stop_codon_frame0_ATGonly <- NA
        } else {
          genome_relative_stop_codon_frame0_ATGonly <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 0)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame0_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x-1):(x-3))})] %>% paste(collapse = ",")
        }
        
      }
      
      # "valid_canonical_ORF_frame_0" = (vector_stranded_CDS_nucleotides[(updated_list_temp$frame_adjust + 1):length(vector_stranded_CDS_nucleotides)] %>% seqsplitter(x = ., n = 3, keep_remainders = FALSE) %>% purrr::map(~.x %>% paste(collapse = "")) %>% unlist %>% grep(pattern = "ATG|CTG|TTG|GTG|ACG|ATA|ATT") %>% .[1]) < (vector_stranded_CDS_nucleotides[(updated_list_temp$frame_adjust + 1):length(vector_stranded_CDS_nucleotides)] %>% seqsplitter(x = ., n = 3, keep_remainders = FALSE) %>% purrr::map(~.x %>% paste(collapse = "")) %>% unlist %>% grep(pattern = "TAG|TAA|TGA") %>% .[1])
      
    } else {
      
      has_ptc_frame0_ATGonly <- "nostart"
      
      length_longest_coding_sequence_frame0_ATGonly <- 0
      
      genome_relative_start_coding_sequence_frame0_ATGonly <- NA
      genome_relative_end_coding_sequence_frame0_ATGonly <- NA
      
      genome_relative_stop_codon_frame0_ATGonly <- NA
      
    }
    
    list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_ATGonly <- list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon %>% purrr::keep(.p = ~any(grepl(x = .x, pattern = "ATG")) == TRUE)
    
    if (length(list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_ATGonly) > 0) {
    
      length_longest_coding_sequence_frame1_ATGonly <- list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist %>% max
      
      genome_relative_start_coding_sequence_frame1_ATGonly <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + ((updated_list_temp$frame_adjust + 1)%%3)]
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - ((updated_list_temp$frame_adjust + 1)%%3)]
      }
      genome_relative_end_coding_sequence_frame1_ATGonly <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 1)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame1_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3] %>% paste(collapse = ",")
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 1)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame1_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3] %>% paste(collapse = ",")
      }
      
      if (updated_list_temp$strand == "+") {
        
        has_ptc_frame1_ATGonly <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame1_ATGonly) < (floor(updated_list_temp$forwards_coords_relative_junction_position) - window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame1_ATGonly) >= (length(updated_list_temp$forward_coords) - 2)) {
          has_ptc_frame1_ATGonly <- "nostop"
          genome_relative_stop_codon_frame1_ATGonly <- NA
        } else {
          genome_relative_stop_codon_frame1_ATGonly <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 1)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame1_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x+1):(x+3))})] %>% paste(collapse = ",")
        }
        
      } else if (updated_list_temp$strand == "-") {
        
        has_ptc_frame1_ATGonly <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame1_ATGonly) > (ceiling(updated_list_temp$forwards_coords_relative_junction_position) + window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame1_ATGonly) <= 3) {
          has_ptc_frame1_ATGonly <- "nostop"
          genome_relative_stop_codon_frame1_ATGonly <- NA
        } else {
          genome_relative_stop_codon_frame1_ATGonly <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 1)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame1_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x-1):(x-3))})] %>% paste(collapse = ",")
        }
        
      }
      
    } else {
        
      has_ptc_frame1_ATGonly <- "nostart"
      
      length_longest_coding_sequence_frame1_ATGonly <- 0
      
      genome_relative_start_coding_sequence_frame1_ATGonly <- NA
      genome_relative_end_coding_sequence_frame1_ATGonly <- NA
      
      genome_relative_stop_codon_frame1_ATGonly <- NA
      
    }
    
    list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_ATGonly <- list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon %>% purrr::keep(.p = ~any(grepl(x = .x, pattern = "ATG")) == TRUE)
    
    if (length(list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_ATGonly) > 0) {
      
      length_longest_coding_sequence_frame2_ATGonly <- list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist %>% max
      
      genome_relative_start_coding_sequence_frame2_ATGonly <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + ((updated_list_temp$frame_adjust + 2)%%3)]
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - ((updated_list_temp$frame_adjust + 2)%%3)]
      }
      genome_relative_end_coding_sequence_frame2_ATGonly <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 2)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame2_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3]
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 2)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame2_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3]
      }
      
      if (updated_list_temp$strand == "+") {
        
        has_ptc_frame2_ATGonly <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame2_ATGonly) < (floor(updated_list_temp$forwards_coords_relative_junction_position) - window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame2_ATGonly) >= (length(updated_list_temp$forward_coords) - 2)) {
          has_ptc_frame2_ATGonly <- "nostop"
          genome_relative_stop_codon_frame2_ATGonly <- NA
        } else {
          genome_relative_stop_codon_frame2_ATGonly <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 2)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame2_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x+1):(x+3))})] %>% paste(collapse = ",")
        }
        
      } else if (updated_list_temp$strand == "-") {
        
        has_ptc_frame2_ATGonly <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame2_ATGonly) > (ceiling(updated_list_temp$forwards_coords_relative_junction_position) + window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame2_ATGonly) <= 3) {
          has_ptc_frame2_ATGonly <- "nostop"
          genome_relative_stop_codon_frame2_ATGonly <- NA
        } else {
          genome_relative_stop_codon_frame2_ATGonly <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 2)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_ATGonly %>% purrr::map(~length(which(.x %in% "ATG")[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame2_ATGonly)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x-1):(x-3))})] %>% paste(collapse = ",")
        }
        
      }
      
    } else {
      
      has_ptc_frame2_ATGonly <- "nostart"
      
      length_longest_coding_sequence_frame2_ATGonly <- 0
      
      genome_relative_start_coding_sequence_frame2_ATGonly <- NA
      genome_relative_end_coding_sequence_frame2_ATGonly <- NA
      
      genome_relative_stop_codon_frame2_ATGonly <- NA
      
    }
    
    list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_noncanonical <- list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon %>% purrr::keep(.p = ~any(grepl(x = .x, pattern = "ATG|CTG|TTG|GTG|ACG|ATA|ATT")) == TRUE)
    
    if (length(list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_noncanonical) > 0) {
      
      length_longest_coding_sequence_frame0_noncanonical <- list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist %>% max
      
      genome_relative_start_coding_sequence_frame0_noncanonical <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + ((updated_list_temp$frame_adjust + 0)%%3)]
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - ((updated_list_temp$frame_adjust + 0)%%3)]
      }
      genome_relative_end_coding_sequence_frame0_noncanonical <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 0)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame0_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3]
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 0)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame0_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3]
      }
      
      if (updated_list_temp$strand == "+") {
        
        has_ptc_frame0_noncanonical <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame0_noncanonical) < (floor(updated_list_temp$forwards_coords_relative_junction_position) - window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame0_noncanonical) >= (length(updated_list_temp$forward_coords) - 2)) {
          has_ptc_frame0_noncanonical <- "nostop"
          genome_relative_stop_codon_frame0_noncanonical <- NA
        } else {
          genome_relative_stop_codon_frame0_noncanonical <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 0)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame0_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x+1):(x+3))})] %>% paste(collapse = ",")
        }
        
      } else if (updated_list_temp$strand == "-") {
        
        has_ptc_frame0_noncanonical <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame0_noncanonical) > (ceiling(updated_list_temp$forwards_coords_relative_junction_position) + window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame0_noncanonical) <= 3) {
          has_ptc_frame0_noncanonical <- "nostop"
          genome_relative_stop_codon_frame0_noncanonical <- NA
        } else {
          genome_relative_stop_codon_frame0_noncanonical <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 0)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_0_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame0_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x-1):(x-3))})] %>% paste(collapse = ",")
        }
        
      }
      
    } else {
      
      has_ptc_frame0_noncanonical <- "nostart"
      
      length_longest_coding_sequence_frame0_noncanonical <- 0
      
      genome_relative_start_coding_sequence_frame0_noncanonical <- NA
      genome_relative_end_coding_sequence_frame0_noncanonical <- NA
      
      genome_relative_stop_codon_frame0_noncanonical <- NA
      
    }
    
    list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_noncanonical <- list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon %>% purrr::keep(.p = ~any(grepl(x = .x, pattern = "ATG|CTG|TTG|GTG|ACG|ATA|ATT")) == TRUE)
    
    if (length(list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_noncanonical) > 0) {
      
      length_longest_coding_sequence_frame1_noncanonical <- list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist %>% max
      
      genome_relative_start_coding_sequence_frame1_noncanonical <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + ((updated_list_temp$frame_adjust + 1)%%3)]
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - ((updated_list_temp$frame_adjust + 1)%%3)]
      }
      genome_relative_end_coding_sequence_frame1_noncanonical <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 1)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame1_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3]
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 1)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame1_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3]
      }
      
      if (updated_list_temp$strand == "+") {
        
        has_ptc_frame1_noncanonical <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame1_noncanonical) < (floor(updated_list_temp$forwards_coords_relative_junction_position) - window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame1_noncanonical) >= (length(updated_list_temp$forward_coords) - 2)) {
          has_ptc_frame1_noncanonical <- "nostop"
          genome_relative_stop_codon_frame1_noncanonical <- NA
        } else {
          genome_relative_stop_codon_frame1_noncanonical <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 1)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame1_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x+1):(x+3))})] %>% paste(collapse = ",")
        }
        
      } else if (updated_list_temp$strand == "-") {
        
        has_ptc_frame1_noncanonical <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame1_noncanonical) > (ceiling(updated_list_temp$forwards_coords_relative_junction_position) + window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame1_noncanonical) <= 3) {
          has_ptc_frame1_noncanonical <- "nostop"
          genome_relative_stop_codon_frame1_noncanonical <- NA
        } else {
          genome_relative_stop_codon_frame1_noncanonical <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 1)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_1_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame1_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x-1):(x-3))})] %>% paste(collapse = ",")
        }
        
      }
      
    } else {
      
      has_ptc_frame1_noncanonical <- "nostart"
      
      length_longest_coding_sequence_frame1_noncanonical <- 0
      
      genome_relative_start_coding_sequence_frame1_noncanonical <- NA
      genome_relative_end_coding_sequence_frame1_noncanonical <- NA
      
      genome_relative_stop_codon_frame1_noncanonical <- NA
      
    }
    
    list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_noncanonical <- list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon %>% purrr::keep(.p = ~any(grepl(x = .x, pattern = "ATG|CTG|TTG|GTG|ACG|ATA|ATT")) == TRUE)
    
    if (length(list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_noncanonical) > 0) {
      
      length_longest_coding_sequence_frame2_noncanonical <- list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist %>% max
      
      genome_relative_start_coding_sequence_frame2_noncanonical <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + ((updated_list_temp$frame_adjust + 2)%%3)]
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - ((updated_list_temp$frame_adjust + 2)%%3)]
      }
      genome_relative_end_coding_sequence_frame2_noncanonical <- if (updated_list_temp$strand == "+") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 2)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame2_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3]
      } else if (updated_list_temp$strand == "-") {
        updated_list_temp$forward_coords[updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 2)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame2_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3]
      }
      
      if (updated_list_temp$strand == "+") {
        
        has_ptc_frame2_noncanonical <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame2_noncanonical) < (floor(updated_list_temp$forwards_coords_relative_junction_position) - window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame2_noncanonical) >= (length(updated_list_temp$forward_coords) - 2)) {
          has_ptc_frame2_noncanonical <- "nostop"
          genome_relative_stop_codon_frame2_noncanonical <- NA
        } else {
          genome_relative_stop_codon_frame2_noncanonical <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon - 1 + ((updated_list_temp$frame_adjust + 2)%%3) + (list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame2_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x+1):(x+3))})] %>% paste(collapse = ",")
        }
        
      } else if (updated_list_temp$strand == "-") {
        
        has_ptc_frame2_noncanonical <- which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame2_noncanonical) > (ceiling(updated_list_temp$forwards_coords_relative_junction_position) + window_size)
        
        if (which(updated_list_temp$forward_coords == genome_relative_end_coding_sequence_frame2_noncanonical) <= 3) {
          has_ptc_frame2_noncanonical <- "nostop"
          genome_relative_stop_codon_frame2_noncanonical <- NA
        } else {
          genome_relative_stop_codon_frame2_noncanonical <- updated_list_temp$forward_coords[(updated_list_temp$fwd_coordinates_relative_first_nt_of_start_codon + 1 - ((updated_list_temp$frame_adjust + 2)%%3) - (list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon[1:(which(list_stranded_CDS_nucleotides_codonised_frame_2_splitbystopcodon_hasstartsonly_noncanonical %>% purrr::map(~length(which(.x %in% c("ATG", "CTG", "TTG", "GTG", "ACG", "ATA", "ATT"))[1]:length(.x))) %>% unlist == length_longest_coding_sequence_frame2_noncanonical)[1] %>% names %>% type.convert(as.is = TRUE))] %>% purrr::map(~.x %>% length) %>% purrr::reduce(sum)) * 3) %>% (function(x){return((x-1):(x-3))})] %>% paste(collapse = ",")
        }
        
      }
      
    } else {
      
      has_ptc_frame2_noncanonical <- "nostart"
      
      length_longest_coding_sequence_frame2_noncanonical <- 0
      
      genome_relative_start_coding_sequence_frame2_noncanonical <- NA
      genome_relative_end_coding_sequence_frame2_noncanonical <- NA
      
      genome_relative_stop_codon_frame2_noncanonical <- NA
      
    }
    
  }
    
    updated_list_temp <- purrr::splice(
      updated_list_temp,
      
      "flag_start_codon_present" = flag_start_codon_present,
      "length_longest_coding_sequence_frame0_ATGonly" = length_longest_coding_sequence_frame0_ATGonly,
      "genome_relative_start_coding_sequence_frame0_ATGonly" = genome_relative_start_coding_sequence_frame0_ATGonly,
      "genome_relative_end_coding_sequence_frame0_ATGonly" = genome_relative_end_coding_sequence_frame0_ATGonly,
      "has_ptc_frame0_ATGonly" = has_ptc_frame0_ATGonly,
      "genome_relative_stop_codon_frame0_ATGonly" = genome_relative_stop_codon_frame0_ATGonly,
      "length_longest_coding_sequence_frame1_ATGonly" = length_longest_coding_sequence_frame1_ATGonly,
      "genome_relative_start_coding_sequence_frame1_ATGonly" = genome_relative_start_coding_sequence_frame1_ATGonly,
      "genome_relative_end_coding_sequence_frame1_ATGonly" = genome_relative_end_coding_sequence_frame1_ATGonly,
      "has_ptc_frame1_ATGonly" = has_ptc_frame1_ATGonly,
      "genome_relative_stop_codon_frame1_ATGonly" = genome_relative_stop_codon_frame1_ATGonly,
      "length_longest_coding_sequence_frame2_ATGonly" = length_longest_coding_sequence_frame2_ATGonly,
      "genome_relative_start_coding_sequence_frame2_ATGonly" = genome_relative_start_coding_sequence_frame2_ATGonly,
      "genome_relative_end_coding_sequence_frame2_ATGonly" = genome_relative_end_coding_sequence_frame2_ATGonly,
      "has_ptc_frame2_ATGonly" = has_ptc_frame2_ATGonly,
      "genome_relative_stop_codon_frame2_ATGonly" = genome_relative_stop_codon_frame2_ATGonly,
      
      "length_longest_coding_sequence_frame0_noncanonical" = length_longest_coding_sequence_frame0_noncanonical,
      "genome_relative_start_coding_sequence_frame0_noncanonical" = genome_relative_start_coding_sequence_frame0_noncanonical,
      "genome_relative_end_coding_sequence_frame0_noncanonical" = genome_relative_end_coding_sequence_frame0_noncanonical,
      "has_ptc_frame0_noncanonical" = has_ptc_frame1_noncanonical,
      "genome_relative_stop_codon_frame0_noncanonical" = genome_relative_stop_codon_frame0_noncanonical,
      "length_longest_coding_sequence_frame1_noncanonical" = length_longest_coding_sequence_frame1_noncanonical,
      "genome_relative_start_coding_sequence_frame1_noncanonical" = genome_relative_start_coding_sequence_frame1_noncanonical,
      "genome_relative_end_coding_sequence_frame1_noncanonical" = genome_relative_end_coding_sequence_frame1_noncanonical,
      "has_ptc_frame1_noncanonical" = has_ptc_frame2_noncanonical,
      "genome_relative_stop_codon_frame1_noncanonical" = genome_relative_stop_codon_frame1_noncanonical,
      "length_longest_coding_sequence_frame2_noncanonical" = length_longest_coding_sequence_frame2_noncanonical,
      "genome_relative_start_coding_sequence_frame2_noncanonical" = genome_relative_start_coding_sequence_frame2_noncanonical,
      "genome_relative_end_coding_sequence_frame2_noncanonical" = genome_relative_end_coding_sequence_frame2_noncanonical,
      "has_ptc_frame2_noncanonical" = has_ptc_frame2_noncanonical,
      "genome_relative_stop_codon_frame2_noncanonical" = genome_relative_stop_codon_frame2_noncanonical
    )
    
    return(updated_list_temp)
    
  }, .progress = TRUE)
  
  # tibblise
  list_transcript_3FT_tibblised <- purrr::map2(
    .x = list_transcript_3FT,
    .y = names(list_transcript_3FT),
    .f = function(a1, a2) {
      
      # DEBUG ###
      # a1 <- list_transcript_3FT[[1]]
      ###########
      
      list_importantdataonly <- a1[(names(a1) %in% c("sdgfsdfsdgdfahethg")) | grepl(x = names(a1), pattern = "frame._ATGonly") | grepl(x = names(a1), pattern = "frame._noncanonical") | grepl(x = names(a1), pattern = "has_ptc")]
      
      tibble_output <- list_importantdataonly %>% tibble::as_tibble() %>% tibble::add_column("transcript_id" = a2, .before = 1)
      
      # (!(a1 %>% purrr::map(~data.class(.x)) %>% unlist) %in% c("list", "tbl_df")) & unlist(a1 %>% purrr::map(~length(.x) == 1))
      
    } )
  
  # collapse and tibblise in preparation for left joins
  tibble_NMD_result_per_chr <- list_transcript_3FT_tibblised %>% rbindlist(use.names = TRUE, fill = TRUE) %>% type_convert(col_types = cols(.default = col_character()), na = c("", "na", "NA")) %>% as_tibble
  
  # assign into the temporary tibble by chromosome
  assign(x = paste("tibble_NMD_result_per_chr_", chr, sep = ""), 
         value = tibble_NMD_result_per_chr)
  
  # SAVE WORKSPACE NOW IF DEBUGGING
  if (save_workspace_when_done == "DEBUG") {
    
    save.image(file = paste(output_dir, "/", output_name, ".Rdata", sep = ""))
    
  }
  
  tictoc::toc()
  
}

# make a list combining the ORF test results of each chromosome
# trim off all those elements (per chromosome) that did not have any 3FT result.
list_NMD_results <- ls(pattern = "tibble_NMD_result_per_chr_") %>% array_tree %>% purrr::map(.x = ., .f = ~get(.x)) %>% compact
# rbindlist into a tibble
tibble_NMD_results <- list_NMD_results %>% rbindlist(use.names = TRUE, fill = TRUE) %>% as_tibble
# table join for final annotated reconstructed GTF
tibble_recon_gtf_nmd <- dplyr::left_join(x = tibble_reconstructed_gtf, y = tibble_NMD_results, by = "transcript_id")

# SAVE WORKSPACE NOW IF DEBUGGING
if (save_workspace_when_done == "DEBUG") {
  
  save.image(file = paste(output_dir, "/", output_name, ".Rdata", sep = ""))
  
}

cat("Annotate exons where stop codons reside")
tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", "contains_PTC_frame0_ATGonly"] <- furrr::future_pmap(
  .l = list(
    "a1" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$start,
    "a2" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$end,
    "a3" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$genome_relative_stop_codon_frame0_ATGonly
  ), 
  .f = function(a1, a2, a3) {
    
    any((a1 < (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE))) & (a2 > (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE)))) %>%
      return
    
  }, .progress = TRUE ) %>% unlist

tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", "contains_PTC_frame1_ATGonly"] <- furrr::future_pmap(
  .l = list(
    "a1" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$start,
    "a2" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$end,
    "a3" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$genome_relative_stop_codon_frame1_ATGonly
  ), 
  .f = function(a1, a2, a3) {
    
    any((a1 < (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE))) & (a2 > (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE)))) %>%
      return
    
  }, .progress = TRUE ) %>% unlist

tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", "contains_PTC_frame2_ATGonly"] <- furrr::future_pmap(
  .l = list(
    "a1" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$start,
    "a2" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$end,
    "a3" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$genome_relative_stop_codon_frame2_ATGonly
  ), 
  .f = function(a1, a2, a3) {
    
    any((a1 < (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE))) & (a2 > (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE)))) %>%
      return
    
  }, .progress = TRUE ) %>% unlist

tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", "contains_PTC_frame0_noncanonical"] <- furrr::future_pmap(
  .l = list(
    "a1" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$start,
    "a2" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$end,
    "a3" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$genome_relative_stop_codon_frame0_noncanonical
  ), 
  .f = function(a1, a2, a3) {
    
    any((a1 < (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE))) & (a2 > (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE)))) %>%
      return
    
  }, .progress = TRUE ) %>% unlist

tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", "contains_PTC_frame1_noncanonical"] <- furrr::future_pmap(
  .l = list(
    "a1" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$start,
    "a2" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$end,
    "a3" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$genome_relative_stop_codon_frame1_noncanonical
  ), 
  .f = function(a1, a2, a3) {
    
    any((a1 < (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE))) & (a2 > (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE)))) %>%
      return
    
  }, .progress = TRUE ) %>% unlist

tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", "contains_PTC_frame2_noncanonical"] <- furrr::future_pmap(
  .l = list(
    "a1" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$start,
    "a2" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$end,
    "a3" = tibble_recon_gtf_nmd[tibble_recon_gtf_nmd$type == "exon", ]$genome_relative_stop_codon_frame2_noncanonical
  ), 
  .f = function(a1, a2, a3) {
    
    any((a1 < (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE))) & (a2 > (a3 %>% strsplit(split = "\\,") %>% unlist %>% type.convert(as.is = TRUE)))) %>%
      return
    
  }, .progress = TRUE ) %>% unlist

# SAVE WORKSPACE NOW IF DEBUGGING
if (save_workspace_when_done == "DEBUG") {
  
  save.image(file = paste(output_dir, "/", output_name, ".Rdata", sep = ""))
  
}

cat("\nlabel first and last exons")

tibble_recon_gtf_nmd_first.last <- furrr::future_map2(
  .x = tibble_recon_gtf_nmd %>% dplyr::group_split(transcript_id),
  .y = 1:length(tibble_recon_gtf_nmd %>% dplyr::group_split(transcript_id)),
  .f = function(a1, a2) {
    
    # DEBUG ###
    # a1 <- tibble_annotated_recon_gtf %>% dplyr::group_split(transcript_id) %>% .[[1]]
    ###########
    
    # print(a2)
    
    output_gtf0 <- a1 %>% tibble::add_column("first_or_last_exon" = "NA")
    
    if (a1$exon_number %>% is.na %>% all == TRUE) {
      
      return(output_gtf0)
      
    } else {
      
      max_exon_number <- max(output_gtf0$exon_number %>% na.omit)
      
      if (max_exon_number == 1) {
        
        output_gtf0[output_gtf0$exon_number == max_exon_number & !is.na(output_gtf0$exon_number), "first_or_last_exon"] <- "only_one_exon"
        
      }
      
      if (any(!a1$strand %in% c("+", "-")) == TRUE) {
        
        return(output_gtf0)
        
      } else {
        
        # account for transcripts which only have one exon
        if (exon_order == "increasing" & output_gtf0$strand %>% unique == "-") {
          
          output_gtf0[output_gtf0$exon_number == max_exon_number & !is.na(output_gtf0$exon_number), "first_or_last_exon"] <- "first_exon"
          
          output_gtf0[output_gtf0$exon_number == 1 & !is.na(output_gtf0$exon_number), "first_or_last_exon"] <- "last_exon"
          
        } else {
          
          output_gtf0[output_gtf0$exon_number == max_exon_number & !is.na(output_gtf0$exon_number), "first_or_last_exon"] <- "last_exon"
          
          output_gtf0[output_gtf0$exon_number == 1 & !is.na(output_gtf0$exon_number), "first_or_last_exon"] <- "first_exon"
          
        }
        
        return(output_gtf0)
        
      }
      
    }
    
  }, .progress = TRUE) %>%
  data.table::rbindlist() %>%
  tibble::as_tibble()

# output_dir <- gsub(x = output_name, pattern = "(.*/)(.*)$", replacement = "\\1")
# output_name <- gsub(x = output_name, pattern = "(.*/)(.*)$", replacement = "\\2")
# 
# # set working directory to prevent screw-up
# setwd(output_dir)

if (save_workspace_when_done == "YES" | save_workspace_when_done == "DEBUG") {
  
  save.image(file = paste(output_dir, "/", output_name, ".RData", sep = ""))
  
}

# write GTF
rtracklayer::export(tibble_recon_gtf_nmd_first.last, con = paste(output_dir, "/", output_name, ".gtf", sep = ""), format = "gtf")

# finish counting
tictoc::toc()

q()


