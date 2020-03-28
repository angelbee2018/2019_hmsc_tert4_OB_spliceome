script_description <- "# NMD FLAGGER ######
Accepts an input of a reconstructed/assembled transcriptome annotation (e.g. strawberry, stringtie) and flags potential candidate NMD transcripts using the 51nt rule as described in <<Hsu, M.K., Lin, H.Y. and Chen, F.C., 2017. NMD Classifier: A reliable and systematic classification tool for nonsense-mediated decay events. PloS one, 12(4)>>.

Behaviour: loads reconstructed transcriptome, then subsets each transcript. For each transcript, look-up the reference FASTA sequence and 3FT to find ALL start codons that occur more than 54nt upstream of the last exon-exon junction. Then look for any in-frame stop codons downstream. Transcript is an NMD candidate when the FIRST nucleotide of the stop codon is within 51 - 50nt of the last exon junction (the stop codon is fully encapsulated in the 51nt window).

Recommended system requirements: 32 threads/64GB memory"

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
  "3" = make_option(c("-O", "--output_path"), type = "character", default = NULL, 
                    help = "Compulsory. output file path. what do you want to save the newly annotated reconstructed GTF as? IMPORTANT: MUST BE A FULL FILE PATH AND NOT A DIRECTORY. e.g. correct: ~/outputdir/annotated_sample.gtf incorrect: ~/outputdir/", metavar = "character"),
  "4" = make_option(c("-C", "--ncores"), type = "integer", default = 0, 
                    help = "Optional. Number of cores to use. possible inputs: numbers 1 to any integer. By default, uses all cores (ncores = 0).", metavar = "integer"),
  "5" = make_option(c("-H", "--chrmode"), type = "integer", default = 0, 
                    help = "Optional. Specifies which chromosomes to do: select what chromosomes you want translated. possible inputs: numbers 0-2. 0 (default): nuclear chromosomes only i,e, 1:22, X & Y. 1: nuclear + mitochondrial i.e. 1:22, X & Y, M. 2: everything including haplotype/fusion chromosomes etc... this is possible provided the chromosome names.", metavar = "integer"),
  "6" = make_option(c("-N", "--nonchrname"), type = "character", default = NULL, 
                    help = "Compulsory only if you have specified \"--chrmode 2\". nonchromosomal file name. if you are doing haplotypes, please specify what the reference genome FASTA file for it is called or the script won't know. This single FASTA file must contain all the haplotype information. The script won't try to search for a second file. In ensembl, this file is called \"Homo_sapiens.GRCh38.dna.nonchromosomal.fa\" or more generally, \"*nonchromosomal.fa\". So for this option, you would specify \"--nonchrname nonchromosomal\".", metavar = "character"),
  "7" = make_option(c("-W", "--checking_window_size"), type = "double", default = 51,
                    help = "ADVANCED. Use this parameter to change the nucleotide window from the default of 51. DO NOT CHANGE THIS UNLESS YOU REALLY WANT TO CHECK USING A DIFFERENT WINDOW SIZE. SUBOPTIMAL RESULTS MAY BE GENERATED AS A RESULT.", metavar = "double")
)

input_arg_info <- OptionParser(option_list = list_input_arg_info, description = script_description)
input_args <- input_arg_info %>% parse_args

# check if the input arguments are O.K
if ((list(input_args$reconstructed_transcript_gtf, input_args$reference_genome_fasta_dir, input_args$output_dir) %>% lapply(is.null) %>% unlist %>% any == TRUE) | 
    (input_args$chrmode == 2 & is.null(input_args$nonchrname) == TRUE)) {
  
  print_help(input_arg_info)
  
  stop("Make sure you entered the arguments correctly", call. = FALSE)
  
}

file_information_path <- input_args$file_information_path
file_information_table <- read_delim(file_information_path, delim = "\t")
# CHECK IF THE FILE INFO TABLE PROVIDED CONTAINS ALL THE COLUMNS WE WANT. STOP IF WE DONT SEE ANY OF THE COLNAMES WE WANT TO FIND.
if (
  any(str_detect(string = colnames(file_information_table), pattern = c("exon_table_path", "reconstructed_GTF_path", "output_database_name", "splicemode_column_name", "IR_regex_string")) == FALSE)
) {
  
  stop("Error. Please check that the file information table is formatted properly. Use --help command for more info.", call. = FALSE)
  
}

# DEBUG #######

# reconstructed_gtf_path <- "Z:/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/alltimepoints_denovo_reconstructed_stringtiemerged.gtf"
# reference_genome_fasta_dir <- "Z:/hg38_ensembl_reference/raw_genome_fasta/genome_fasta_extract2/"
# output_path <- "Z:/PGNEXUS_kassem_MSC/Kassem_OB/analysis_NMD_classifier/results/alltimepoints_denovo_reconstructed_stringtiemerged_NMDflagged.gtf"
# window_size <- 51

###############

reconstructed_gtf_path <- input_args$reconstructed_transcript_gtf
reference_genome_fasta_dir <- input_args$reference_genome_fasta_dir
output_path <- input_args$output_path
window_size <- input_args$checking_window_size

cat("reconstructed_gtf_path:", reconstructed_gtf_path, "\n")
cat("reference_genome_fasta_dir:", reference_genome_fasta_dir, "\n")
cat("output_path:", output_path, "\n")
cat("window_size:", window_size, "\n")

if(!dir.exists(output_path) ) {
  dir.create(output_path, recursive = TRUE)}

# manage parrallellisation rrlllRll

if (input_args$ncores != 0) {
  number_of_workers <- input_args$ncores
} else {
  number_of_workers <- future::availableCores()
} 

cat(future::availableCores(), "cores will be used\n")
future::plan(multiprocess)
options(future.globals.maxSize = 30000000000, mc.cores = number_of_workers)

# DEFINE FUNCTIONS ##########################

# FUNCTION TO 3 FRAME TRANSLATE ONE LIST CONTAINING NUCLEOTIDE SEQUENCE AND STRAND

nt.sequence_strand_threeframetranslate <- function(vector_forward_nucleotides, strand) {
  
  if (strand == "+") {
    
    translation_result <- list("translation_frame_0" = seqinr::translate(vector_forward_nucleotides, frame = 0, sens = "F"),
                               "translation_frame_1" = seqinr::translate(vector_forward_nucleotides, frame = 1, sens = "F"),
                               "translation_frame_2" = seqinr::translate(vector_forward_nucleotides, frame = 2, sens = "F"))
    
  } else if (strand == "-") {
    
    translation_result <- list("translation_frame_0" = seqinr::translate(vector_forward_nucleotides, frame = 0, sens = "R"),
                               "translation_frame_1" = seqinr::translate(vector_forward_nucleotides, frame = 1, sens = "R"),
                               "translation_frame_2" = seqinr::translate(vector_forward_nucleotides, frame = 2, sens = "R"))
    
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

# BEGIN EXECUTION ###########################

cat("import reconstructed transcriptome GTF")
reconstructed_gtf <- rtracklayer::import(reconstructed_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)

vector_ref_genome_paths_by_chr <- paste(reference_genome_fasta_dir, list.files(reference_genome_fasta_dir)[list.files(reference_genome_fasta_dir) %>% grep(., pattern = ".*.fa$")], sep = "")

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
  
  reconstructed_gtf_tempchr <- reconstructed_gtf[reconstructed_gtf$seqnames == chr, ]
  
  cat("subset reconstructed transcriptome GTF by filtering for exon entries only.\n")
  row_indices_recon.gtf_exon.entries.only <- which(reconstructed_gtf_tempchr$type == "exon")
  reconstructed_gtf_exon.entries.only <- reconstructed_gtf_tempchr[row_indices_recon.gtf_exon.entries.only, ]
  
  # leave aside the remaining entries for later inclusion
  row_indices_recon.gtf_non.exon.entries <- which(reconstructed_gtf_tempchr$type != "exon")
  
  # generate list of all 'transcript_id's to loop thru
  list_transcript_ids <- reconstructed_gtf_tempchr$transcript_id %>% unique %>% array_tree
  
  cat("load the exon entries of the recon. GTF table into the list\n")
  list_reconstructed_gtf_exon.entries.only_by_transcript_id <- future_map(.x = list_transcript_ids, .f = ~reconstructed_gtf_exon.entries.only[reconstructed_gtf_exon.entries.only$transcript_id == .x, ], .progress = TRUE, .options = future_options(globals = c("reconstructed_gtf_exon.entries.only"))) %>% set_names(reconstructed_gtf_tempchr$transcript_id %>% unique)
  
  # filter out transcript entries which only have one exon
  list_reconstructed_gtf_exon.entries.only_by_transcript_id <- purrr::discard(.x = list_reconstructed_gtf_exon.entries.only_by_transcript_id, .p = ~length(.x$exon_number) < 2)
  
  cat("generate attributes of each transcript:
      1. all genome-relative forward coords of each nucleotide,
      2. coords of the nucleotides flanking last exon-exon junction\n")
  list_transcript_attributes <- future_map(.x = list_reconstructed_gtf_exon.entries.only_by_transcript_id, 
                                           .f = ~list(
                                             "recon_entries" = .x,
                                             "strand" = .x$strand %>% unique %>% paste,
                                             "forward_coords" = purrr::map2(.x = .x$start, .y = .x$end, .f = ~.x:.y) %>% unlist %>% as.numeric %>% sort, 
                                             # since the start is always less than the end in genome relative coords, for the fwd strand last_nt_of_second_last_exon will be the end coord of the second last exon.
                                             # however for rev. strand, last_nt_of_second_last_exon will be the start coord of the second last exon, and the first_nt_of_last_exon will be the end coord of the last exon.
                                             "last_nt_of_second_last_exon" = if (.x$strand %>% unique %>% paste == "+") {
                                               .x[.x$exon_number == max(.x$exon_number %>% as.numeric) - 1, "end"] %>% paste %>% as.numeric
                                             } else if (.x$strand %>% unique %>% paste == "-") {
                                               .x[.x$exon_number == min(.x$exon_number %>% as.numeric) + 1, "start"] %>% paste %>% as.numeric
                                             } , 
                                             "first_nt_of_last_exon" = if (.x$strand %>% unique %>% paste == "+") {
                                               .x[.x$exon_number == max(.x$exon_number %>% as.numeric), "start"] %>% paste %>% as.numeric
                                             } else if (.x$strand %>% unique %>% paste == "-") {
                                               .x[.x$exon_number == min(.x$exon_number %>% as.numeric), "end"] %>% paste %>% as.numeric
                                             }
                                           ), .progress = TRUE
                                           
  )
  
  cat("generate all the forward nucleotides of each transcript\n")
  # also get the TRANSCRIPT RELATIVE position of the junction  (transcript-relative coord will be a half-integer value)
  list_transcript_fwd_nts <- future_imap(.x = list_transcript_attributes, .f = function(.x, .y) {
    
    # cat("now looking up nucleotides of transcript number", which(names(list_transcript_attributes) == .y), "/", length(list_transcript_attributes), "\n")
    
    updated_list <- purrr::splice(.x,
                                  
                                  "forward_nucleotides" = reference_genome_fasta_chr_temp[[chr %>% paste]][.x$forward_coords],
                                  "transcript_relative_junction_position" = mean(c(
                                    # transcript relative position of the last nt of second last exon
                                    which(.x$forward_coords == .x$last_nt_of_second_last_exon),
                                    # transcript relative position of the first nt of last exon
                                    which(.x$forward_coords == .x$first_nt_of_last_exon))))
    
    return(updated_list)
    
  }, .progress = TRUE)
  
  cat("do 3FT based on the forward nucleotides of each transcript and test for whether there is any valid ORF.\n")
  list_transcript_3FT <- future_imap(.x = list_transcript_fwd_nts, .f = function(.x, .y) {
    
    cat("now checking 3FT of transcript number", which(names(list_transcript_fwd_nts) == .y), "/", length(list_transcript_fwd_nts), "\n")
    
    # for simplicity, we exclude the last exon + 51 nt window from the nt sequence
    # hence get the effective fwd. nucleotide sequence:
    vector_all_foward_nucleotides <- .x$forward_nucleotides
    
    # if (weirdly enough) the transcript is actually shorter than 51 nt, then we dont care about the window. Instances of this occurring are very rare indeed.
    # requirement is that the sequence length longer than the window is at least window + 5 because you need at least 5nt in order for 3FT to work succesfully in all frames
    if (length(vector_all_foward_nucleotides) < (window_size + 1) |
        ((.x$strand == "+") & (floor(.x$transcript_relative_junction_position) - window_size) < 5) |
        ((.x$strand == "-") & (ceiling(.x$transcript_relative_junction_position) + window_size) > (length(vector_all_foward_nucleotides) - 5))) {
      
      vector_forward_nucleotides_excl.last.exon <- vector_all_foward_nucleotides
      
    } else if (.x$strand == "+") {
      
      vector_forward_nucleotides_excl.last.exon <- vector_all_foward_nucleotides[1:(floor(.x$transcript_relative_junction_position) - window_size)]
      
    } else if (.x$strand == "-") {
      
      vector_forward_nucleotides_excl.last.exon <- vector_all_foward_nucleotides[(ceiling(.x$transcript_relative_junction_position) + window_size):length(vector_all_foward_nucleotides)]
      
    }
    
    # splice in the 3FT data
    updated_list_temp <- purrr::splice(.x, nt.sequence_strand_threeframetranslate(vector_forward_nucleotides = vector_forward_nucleotides_excl.last.exon, strand = .x$strand))
    # test for the presence of valid ORFs
    updated_list_temp <- purrr::splice(updated_list_temp, 
                                       "any_valid_ORF_frame_0" = test_for_any_valid_ORF(updated_list_temp$translation_frame_0),
                                       "any_valid_ORF_frame_1" = test_for_any_valid_ORF(updated_list_temp$translation_frame_1),
                                       "any_valid_ORF_frame_2" = test_for_any_valid_ORF(updated_list_temp$translation_frame_2))
    # collapse the vectors into string format
    ## commas to separate numbers
    updated_list_temp <- purrr::modify_at(.x = updated_list_temp, .at = "forward_coords", .f = ~paste(.x, collapse = ","))
    ## nothing to separate nt. and AA sequences.
    updated_list <- purrr::modify_at(.x = updated_list_temp, .at = c("forward_nucleotides", "translation_frame_0", "translation_frame_1", "translation_frame_2"), .f = ~paste(.x, collapse = ""))
    
    return(updated_list)
    
  }, .progress = TRUE)
  
  # get the transcripts which have all FALSE valid ORFs
  list_all_invalid_ORFs_by_transcript <- purrr::keep(.x = list_transcript_3FT, .p = ~all(c(.x$any_valid_ORF_frame_0, .x$any_valid_ORF_frame_1, .x$any_valid_ORF_frame_2) == FALSE))
  
  # these transcripts with no valid ORFs are NMD candidates.
  list_logical_NMD_or_not <- purrr::map2(.x = list_all_invalid_ORFs_by_transcript, .y = names(list_all_invalid_ORFs_by_transcript), .f = ~tibble("transcript_id" = .y, NMD_candidate = "TRUE"))
  
  # collapse and tibblise in preparation for left joins
  tibble_logical_NMD_or_not <- list_logical_NMD_or_not %>% rbindlist %>% type_convert %>% as_tibble
  
  # assign into the temporary tibble by chromosome
  assign(x = paste("tibble_logical_NMD_or_not_temp_chr_", chr, sep = ""), 
         value = tibble_logical_NMD_or_not)
  
  tictoc::toc()
  
}

# make a list combining the ORF test results of each chromosome
# trim off all those elements (per chromosome) that did not have any 3FT result.
list_ORF_test <- ls(pattern = "tibble_logical_NMD_or_not_temp_chr_") %>% array_tree %>% purrr::map(.x = ., .f = ~get(.x)) %>% compact
# rbindlist into a tibble
tibble_ORF_test <- list_ORF_test %>% rbindlist %>% as_tibble
# table join for final annotated reconstructed GTF
tibble_annotated_recon_gtf <- dplyr::left_join(x = reconstructed_gtf, y = tibble_ORF_test, by = "transcript_id")
# turn NA in the NMD column to FALSE
tibble_annotated_recon_gtf[which(is.na(tibble_annotated_recon_gtf$NMD_candidate)), "NMD_candidate"] <- FALSE

# write GTF
rtracklayer::export(object = tibble_annotated_recon_gtf, con = output_path, format = "gtf")
  
# print stats summary

cat("number of transcripts identified as NMD:", length(tibble_ORF_test$transcript_id %>% unique), "\n")
cat("number of transcripts which were in the GTF:", length(reconstructed_gtf$transcript_id %>% unique), "\n")
cat("percentage of transcripts identified as NMD candidates:", ((length(tibble_ORF_test$transcript_id %>% unique)) / (length(reconstructed_gtf$transcript_id %>% unique))) * 100, "%\n")

# finish counting
tictoc::toc()

q()


