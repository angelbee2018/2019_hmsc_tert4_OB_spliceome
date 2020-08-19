script_description <- "# THREE-FRAME TRANSLATION ######
Attempts to translate the nucleotides of all transcripts flanking all the splice junctions given.
Uses PARALLEL PURRR (FURRR) ^___^ is much faster than the last version
Recommended CPU/RAM connsumption: 12-16 cores/60GB for ~30,000 junctions. 4-8 cores/32GB for ~4000 junctions. 2-4 cores/16GB for ~1000 junctions. 
I personally would not use more than 16 cores because you get diminishing returns due to the longer time it takes to delegate the tasks to each worker.
RAM usage  scales by the number of cores you use. I do not recommend going lower than 16GB."

# print the arguments received by the R script
cat("Arguments input:", commandArgs(), sep = "\n")
args = 
  commandArgs(trailingOnly = TRUE)
cat(args)
cat("number of arguments specified:", length(args))

# SET ENVIRONMENT ##########
#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")
#BiocManager::install(c("seqinr", "tidyverse", "purrr", "dplyr", "rtracklayer", "data.table", "furrr", "RhpcBLASctl", "optparse"))

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
  "1" = make_option(c("-J", "--junction_table_path"), type = "character", default = NULL, 
              help = "Compulsory. path to table containing junctions of interest. for example, UNION_junc_coor type files (JUM only). MUST contain columns start, end, chr, strand, and splicemode. OPTIONAL colnames: gene_name, organism, custom_identifier. one row per junction.", metavar = "character"),
  "2" = make_option(c("-I", "--intron_retention_string"), type = "character", default = "intron_retention", 
                    help = "Compulsory. A regular expression which matches to all characters in the splicemode column which are associated with IR events.", metavar = "character"),
  "3" = make_option(c("-S", "--source_tag"), type = "character", default = "three_frame_translation_junctions", 
                    help = "Compulsory. A character string that will be added to the FASTA headers to indicate the source. It is in the same position as \"sp\" for UniProt fasta files.", metavar = "character"),
  "4" = make_option(c("-G", "--reconstructed_gtf_path"), type = "character", default = NULL, 
                    help = "Compulsory. path to the actual reconstructed GTF file (e.g. from Cufflinks, Strawberry). NOT THE CONTAINING DIRECTORY. tip: for better junction matching, combine the reconstructed GTF with reference GTF beforehand e.g. using StringTie", metavar = "character"),
  "5" = make_option(c("-R", "--reference_genome_fasta_dir"), type = "character", default = NULL, 
              help = "Compulsory. path to the directory containing the genome FASTA files. Ideally from Ensembl... you need separate files by chromosomes, NOT the primary assembly. 
              FORMATTING IMPORTANT!!!! MAKE SURE THE REF. GENOME FASTA FILES ARE IN THE FORMAT: <_anything_><chr>.fa e.g. \"Homo_sapiens.GRCh38.dna.chromosome.MT.fa\" OR \"chr16.fa\" OR \"Y.fa\". What will not work: anything which does not have .fa extension e.g. \"chr16.fasta\", anything between the chromosome number and the .fa extension e.g. \"chromosome1.ensembl.fa\"", metavar = "character"),
  "6" = make_option(c("-U", "--upstream_window_size"), type = "double", default = 50,
              help = "Optional. how many nucleotides to translate upstream W.R.T. the middle of the exon for junction-based mode, the total nt. length to be translated will be arg3 + arg4. default for both is 50.", metavar = "double"),
  "7" = make_option(c("-D", "--downstream_window_size"), type = "double", default = 50,
  help = "Optional. how many nt to translate downstream W.R.T. the transcript. for exon-based mode, both will be 50 by default. instead, the window size is taken as the MINIMUM required length to be considered for translation starting from the middle of the exon.", metavar = "double"),
  "8" = make_option(c("-T", "--output_name"), type = "character", default = NULL,
                    help = "Compulsory. a character string of what the final FASTA database file name and output table will be called. a.k.a. what do you want to save the FASTA as? IMPORTANT: MUST BE A STRING WITHOUT THE EXTENSION AND NOT A DIRECTORY. THE .txt EXTENSION WILL AUTOMATICALLY BE ADDED FOR THE OUTPUT FILE. e.g. correct: custom_database incorrect: custom_database.fasta incorrect: custom_database/", metavar = "character"),
  "9" = make_option(c("-O", "--output_dir"), type = "character", default = NULL, 
              help = "Compulsory. output directory. where do you want to save the custom databases? IMPORTANT: directory must end in a \"/\". e.g. correct: ~/outputdir/ incorrect: ~/outputdir", metavar = "character"),
  "10" = make_option(c("-C", "--ncores"), type = "character", default = 0, 
                  help = "Optional. Number of cores to use. possible inputs: numbers 1 to any integer. By default, uses all cores (ncores = 0). If a single number is specified, it will just tell future to loop thru chromosomes in parallel using the specified core count. If numberxnumber for example 7x4 then 28 cores will be used. 7 for chromosomes and 4 for inside each chromosome.", metavar = "character"),
  "11" = make_option(c("-H", "--chrmode"), type = "integer", default = 0, 
                  help = "Optional. Specifies which chromosomes to do: select what chromosomes you want translated. possible inputs: numbers 0-2. 0 (default): nuclear chromosomes only i,e, 1:22, X & Y. 1: nuclear + mitochondrial i.e. 1:22, X & Y, M. 2: everything including haplotype/fusion chromosomes etc... this is possible provided the chromosome names.", metavar = "integer"),
  "12" = make_option(c("-N", "--nonchrname"), type = "character", default = NULL, 
                    help = "Compulsory only if you have specified \"--chrmode 2\". nonchromosomal file name. if you are doing haplotypes, please specify what the reference genome FASTA file for it is called or the script won't know. This single FASTA file must contain all the haplotype information. The script won't try to search for a second file. In ensembl, this file is called \"Homo_sapiens.GRCh38.dna.nonchromosomal.fa\" or more generally, \"*nonchromosomal.fa\". So for this option, you would specify \"--nonchrname nonchromosomal\".", metavar = "character")
)

input_arg_info <- OptionParser(option_list = list_input_arg_info, description = script_description)
input_args <- input_arg_info %>% parse_args

# check if the input arguments are O.K
if ((list(input_args$junction_table_path, input_args$reconstructed_gtf_path, input_args$reference_genome_fasta_dir, input_args$output_name, input_args$output_dir) %>% lapply(is.null) %>% unlist %>% any == TRUE) | 
    (input_args$chrmode == 2 & is.null(input_args$nonchrname) == TRUE)) {
  
  print_help(input_arg_info)
  
  stop("Make sure you entered the arguments correctly", call. = FALSE)
  
}

junction_table_path <- input_args$junction_table_path
intron_retention_string <- input_args$intron_retention_string
source_tag <- input_args$source_tag
reconstructed_gtf_path <- input_args$reconstructed_gtf_path
reference_genome_fasta_dir <- input_args$reference_genome_fasta_dir
upstream_window_size <- input_args$upstream_window_size
downstream_window_size <- input_args$downstream_window_size
output_name <- input_args$output_name
output_dir <- input_args$output_dir
ncores <- input_args$ncores
chrmode <- input_args$chrmode
nonchrname <- input_args$nonchrname

# DEBUG ########
# tibble_JUM_diff_table <- read.delim("/media/Ubuntu/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_JUM/run_2_PGNEXUS_OBseries_readlength100/R_processing_results/wide_table_of_7855_constitutive_VSRs_dPSI_OB_diff_qvalue0.01_dPSI0.15_no_na.txt", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>% as_tibble
# # strsplit into a chr start end strand tibble
# list_constituent_junctions <- tibble_JUM_diff_table[, c("Gene", "splicemode", "chr", "start", "end", "strand")] %>% 
#   array_tree %>%
#   # strsplit the chr, strand, start and end, and leave the rest alone.
#   purrr::map(.f = function(a1) {
#     
#     # DEBUG ###
#     # a1 <- tibble_JUM_diff_table[, c("Gene", "splicemode", "chr", "start", "end", "strand")] %>% 
#       # array_tree %>%
#       # .[[1]]
#     ###########
#     
#     # strsplit the chr, start, end, strand elements.
#     tibble_split_chr_start_end_strand <- a1[c("chr", "start", "end", "strand")] %>% 
#       purrr::map(~strsplit(.x, split = ";") %>% unlist)
#     
#     # combine the gene and splicemode elements, tibblise and return
#     purrr::splice(a1[c("Gene", "splicemode")],
#                   tibble_split_chr_start_end_strand) %>%
#       as_tibble %>% 
#       return
#   })
# 
# # rbind and tibblise
# tibble_constituent_junctions <- list_constituent_junctions %>% rbindlist(use.names = TRUE) %>% as_tibble
# 
# tibble_junction_table <- tibble_constituent_junctions %>% 
#   dplyr::rename("gene_name" = "Gene") %>%
#   add_column("source" = "JUM_strawberry",
#              "organism" = "Homo sapiens",
#              "custom_identifier" = NA)
# 
# junction_table_path <- "/media/Ubuntu/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_JUM/run_2_PGNEXUS_OBseries_readlength100/R_processing_results/wide_table_of_983_differential_VSRs_qvalue0.01_dPSI0.15_with_na_constituent_junctions.txt"
# intron_retention_string <- "intron_retention"
# reconstructed_gtf_path <- "/media/Ubuntu/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/GRAND_OBseries_ref_denovo_reconstructed_stringtiemerged.gtf"
# reference_genome_fasta_dir <- "/media/Ubuntu/sharedfolder/hg38_ensembl_reference/raw_genome_fasta/genome_fasta_extract2/"
# upstream_window_size <- 50
# downstream_window_size <- 50
# output_name <- "test_differential_JUM_strawberry"
# output_dir <- "/media/Ubuntu/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/results_database_generation/"
# ncores <- 4
# chrmode <- 1
# nonchrname <- NULL
##################################

cat("junction_table_path:", junction_table_path, "\n")
cat("intron_retention_string:", intron_retention_string, "\n")
cat("source_tag:", source_tag, "\n")
cat("reconstructed_gtf_path:", reconstructed_gtf_path, "\n")
cat("reference_genome_fasta_dir:", reference_genome_fasta_dir, "\n")
cat("upstream_window_size:", upstream_window_size, "\n")
cat("downstream_window_size:", downstream_window_size, "\n")
cat("output_name:", output_name, "\n")
cat("output_dir:", output_dir, "\n")
cat("ncores:", ncores, "\n")
cat("chrmode:", chrmode, "\n")
cat("nonchrname:", nonchrname, "\n")

if(!dir.exists(output_dir) ) {
  dir.create(output_dir, recursive = TRUE)}

# Open a file to send messages to
message_divert_path <- file(paste(output_dir, "/", output_name, "_messages.txt", sep = ""), open = "wt")
# Divert messages to that file
sink(message_divert_path, type = "message")

# manage parrallellisation rrlllRll

if (grepl(x = ncores, pattern = "x") == FALSE) {
  
  if (ncores != 0) {
    number_of_workers <- ncores
    cat(future::availableCores(), "cores will be used\n")
  } else {
    number_of_workers <- future::availableCores()
    cat(future::availableCores(), "cores will be used\n")
  } 
  
} else if (grepl(x = ncores, pattern = "x") == TRUE) {
  
  plan(list(tweak(multiprocess, workers = ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert), 
            tweak(multiprocess, workers = ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert))
  )
  
  cat((ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert) * (ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert), "cores will be used in total\n")
  cat("first layer:", ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert, "cores\n")
  cat("second layer:", ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert, "cores\n")
  
}

options(future.globals.maxSize = 30000000000, future.fork.enable = TRUE)

# set layered parallelisation in furrr

# if ((number_of_workers - length(vector_chr_in_commmon)) %/% number_of_workers > 1) {
#   plan(list(tweak(multiprocess, workers = min(number_of_workers, length(vector_chr_in_commmon))), 
#             tweak(multiprocess, workers = (number_of_workers - length(vector_chr_in_commmon)) %/% number_of_workers))
# )
# } else {
#   plan(multiprocess, workers = number_of_workers)
# }

# activate the below if linux R starts acting up

# if(Sys.info()["sysname"] == "Windows") {
#   
#   cat("number of workers:", number_of_workers, "\n")
#   
#   future::plan(multiprocess)
#   options(future.globals.maxSize = 30000000000, mc.cores = input_args$ncores)
#   
# } else {
#   
#   cat("number of workers:", number_of_workers, "\n")
#   
#   future::plan(multisession, workers = number_of_workers)
#   # library(RhpcBLASctl)
#   # RhpcBLASctl::omp_set_num_threads(number_of_workers)
#   # RhpcBLASctl::blas_set_num_threads(number_of_workers)
#   options(future.globals.maxSize = 30000000000, mc.cores = number_of_workers)
#   
# }

# specify the chromosomes to be run, according to user option --chrmode
if (chrmode == 1) {
  
  # all nuclear chromosomes + mitochondria
  chr_to_run <- c(1:22, "X", "Y", "MT")
  
} else if (chrmode == 2) {
  
  chr_to_run <- c(1:22, "X", "Y", "MT", nonchrname)
  
} else {
  
  # if the user put in a stupid number then we'll assume they just want all the nuclear chromosomes.
  chr_to_run <- c(1:22, "X", "Y")
  
}

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

# END nt.sequence_strand_threeframetranslate

# FUNCTION TO EXTRACT TRANSCRIPTS WITH JUNCTION-FLANKING EXONS.
# NOTE: to be used with purrr
# input: spliceregion_list: a list containing details of ONE junction: $diff_exon_chr, $diff_exon_start, $diff_exon_end
# tibble_gtf_table: rtracklayer::import(...) %>% as_tibble. can be any GTF. reconstructed or reference.
extract_junction.flanking.exons <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, tolerance_left = 1, tolerance_right = 1, tolerance_inside = 1, tolerance_outside = 0, match_consecutive = TRUE, return_type = "exon") {
  
  # DEBUG ###################
  
  # query_chr = a1$chr %>% type.convert
  # query_start = a1$event_region_start %>% type.convert
  # query_end = a1$event_region_end %>% type.convert
  # query_strand = "*"
  # tibble_gtf_table = tibble_ref_gtf
  # tolerance_left = 0
  # tolerance_right = 0
  # tolerance_inside = 0
  # tolerance_outside = 0
  # match_consecutive = FALSE
  # return_type = "exon"
  
  ###########################
  
  # print(paste("now processing junction number", index))
  
  if (query_strand == "." | query_strand == 0 | query_strand == "*") {
    
    tibble_gtf_subset_flanking_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% .[.$start <= ((query_end %>% as.numeric) + 1 + tolerance_outside + tolerance_left) & .$end >= ((query_start %>% as.numeric) - 1 - tolerance_outside - tolerance_left), ] %>% .[!(.$start <= ((query_end %>% as.numeric) - tolerance_inside - tolerance_right) & .$end >= ((query_start %>% as.numeric) + tolerance_inside + tolerance_left)), ] %>% .[.$type %in% return_type, ]
    
  } else if (query_strand == "+" | query_strand == "-") {
    
    tibble_gtf_subset_flanking_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% .[.$strand == query_strand %>% trimws, ] %>% .[.$start <= ((query_end %>% as.numeric) + 1 + tolerance_right) & .$end >= ((query_start %>% as.numeric) - 1 - tolerance_left), ] %>% .[!(.$start <= ((query_end %>% as.numeric) - tolerance_inside - tolerance_right) & .$end >= ((query_start %>% as.numeric) + tolerance_inside + tolerance_left)), ] %>% .[.$type %in% return_type, ]
    
  } else {
    
    stop("Could not match the strand information in the transposed differential-only UNION_junc_coor_table. Make sure that the \"strand\" column in the UNION_junc_coor_table contains only +, - or .")
    
  }
  
  list_of_junction_associated_transcripts <- tibble_gtf_subset_flanking_exons$transcript_id %>% unique %>% array_tree %>% flatten
  
  # make a list for each transcript that directly flanks a junction.
  # then filter so that there are only a) exon PAIRS which b) are directly connected in the mature (spliced) transcript
  
  if (match_consecutive == TRUE) {
    
    list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_gtf_subset_flanking_exons[tibble_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts) %>% keep(.x = ., .p = ~nrow(.x) == 2) %>% keep(.x = ., .p = ~abs((.x[2, "exon_number"] %>% paste %>% as.numeric) - (.x[1, "exon_number"] %>% paste %>% as.numeric)) == 1)
    
  } else if (match_consecutive == FALSE) {
    
    list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_gtf_subset_flanking_exons[tibble_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts) %>% keep(.x = ., .p = ~nrow(.x) == 2)
    
  }
  
  # for those consecutive exons which were found to flank a junction, get all the entries of the parent transcript
  list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_tibbles_flanking_exon_gtf.entries_per_transcript, .f = ~list(
    
    "matched_flanking_exons" = .x, 
    "parent_transcript" = tibble_gtf_table[tibble_gtf_table$transcript_id == .x$transcript_id %>% unique %>% paste, ] %>% .[-which(is.na(.$exon_number)), ] %>%
      dplyr::arrange(exon_number %>% as.numeric)))
  
  return(list_of_tibbles_flanking_exon_gtf.entries_per_transcript)
  
}

# END extract_junction.flanking.exons_JUM() ###

# FUNCTION to calculate where the exon's amino acid sequence will be within a whole stretch of translated transcript
# which is the ceiling of the distance /3 - translation frame + 1. reverse for reverse strand.
calculate_translation_frame_relative_start_end_position <- function(ES, EE, TL, strand, frame) {
  # ES: exon start (transcript-relative nucleotide position), EE: exon end, TL: transcript length, frame: 0-2
  
  if (strand == "+") {
    
    exon_start_AA_position <- ceiling((ES - 1 - frame) / 3) + 1
    exon_end_AA_position <- floor((EE - 1 - frame) / 3)
    
  } else if (strand == "-") {
    
    exon_start_AA_position <- ceiling((TL - EE - frame) / 3) + 1
    exon_end_AA_position <- floor((TL - ES - frame) / 3)
    
  }
  
  return(list("exon_start_AA_position" = exon_start_AA_position, "exon_end_AA_position" = exon_end_AA_position))
  
}

# FUNCTIONS to test if the left/right side of stop codons in an exon are translatable or not. (i.e. whether uORF or dORF exists or not)
find_valid_uORF <- function(list) {
  
  AA_sequence <- list[[1]] %>% unlist
  window_start_AA_position <- list[[2]] %>% paste %>% as.numeric
  window_end_AA_position <- list[[3]] %>% paste %>% as.numeric
  junction_AA_position <- list[[4]] %>% paste %>% as.numeric 
  
  validity_test <- stringr::str_detect(AA_sequence[1:(window_start_AA_position - 1)] %>% rev %>% paste(collapse = ""), "^[^\\*]+M|^[^\\*]+$")
  
  exonic_uORF_sequence <- AA_sequence[window_start_AA_position:window_end_AA_position] %>% paste(collapse = "") %>% strsplit(., split = "\\*") %>% unlist %>% first
  
  # if there is indeed a valid uORF then translate the first part within the exon
  if (validity_test == TRUE & nchar(exonic_uORF_sequence) > (junction_AA_position - window_start_AA_position + 1)) {
    
    return(exonic_uORF_sequence)
    
  } else {
    
    return("NONE_VALID")
    
  }
  
}


find_valid_dORF <- function(list) {
  
  AA_sequence <- list[[1]] %>% unlist
  window_start_AA_position <- list[[2]] %>% paste %>% as.numeric
  window_end_AA_position <- list[[3]] %>% paste %>% as.numeric
  junction_AA_position <- list[[4]] %>% paste %>% as.numeric 
  
  validity_test <- stringr::str_detect(AA_sequence[window_start_AA_position:window_end_AA_position] %>% rev %>% paste(collapse = ""), "^[^\\*]+M|^[^\\*]+$")
  
  exonic_dORF_sequence <- AA_sequence[window_start_AA_position:window_end_AA_position] %>% paste(collapse = "") %>% strsplit(., split = "\\*") %>% unlist %>% last
  
  # if there is indeed a valid dORF then translate the first part within the exon
  if (validity_test == TRUE & nchar(exonic_dORF_sequence) > (junction_AA_position - window_start_AA_position + 1)) {
    
    return(exonic_dORF_sequence)
    
  } else {
    
    return("NONE_VALID")
    
  }
  
}

# END find_valid_uORF() and find_valid_dORF()

# MAIN FUNCTION TO DO 3FT OF JUNCTIONS
# Behaviour: for each junction, look up every single transcript that it's associated with. Translate in three frames and find which frame can validly translate the exon.
# These transcripts are matched according to both reference and reconstructed GTF.
# The smallest protein in humans is 44 AA so the smallest valid translatable nucleotide length is 132 nt. 
# Translated region included in the upstream/downstream window must cross the splice junction.

# END FUNCTIONS #####################################################################################################
#####################################################################################################################
#####################################################################################################################

cat("# BEGIN EXECUTION #################################\n")

vector_ref_genome_paths_by_chr <- paste(reference_genome_fasta_dir, list.files(reference_genome_fasta_dir)[list.files(reference_genome_fasta_dir) %>% grep(., pattern = ".*.fa$")], sep = "")

# reference_genome_fasta <- seqinr::read.fasta(file = reference_genome_fasta_path, forceDNAtolower = FALSE)

# total_peptide_window_size <- upstream_window_size + downstream_window_size
# cat("total translation window size:", total_peptide_window_size, "\n")

tibble_reconstructed_gtf <- rtracklayer::import(reconstructed_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character) %>% type_convert
# list-ify the recon GTF by chromosome
list_recon_gtf_subset_by_chr <- tibble_reconstructed_gtf %>% dplyr::group_split(seqnames)

names(list_recon_gtf_subset_by_chr) <- list_recon_gtf_subset_by_chr %>% purrr::map(.f = ~.x$seqnames %>% unique) %>% unlist

# filter chr for only user specified chr
list_recon_gtf_subset_by_chr <- list_recon_gtf_subset_by_chr[chr_to_run]

cat("GTF importing done\n")

tibble_junction_table <- read.delim(junction_table_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>% as_tibble %>% unique
# %>% .[sample(1:nrow(.), size = 100),]
# list-ify the junction table by chromosome
list_junction_table_by_chr <- tibble_junction_table %>% dplyr::group_split(chr)

# filter chr for only user specified chr
names(list_junction_table_by_chr) <- list_junction_table_by_chr %>% purrr::map(~.$chr %>% unique) %>% unlist
list_junction_table_by_chr <- list_junction_table_by_chr[chr_to_run]

# subset the recon GTF and junction table by chromosomes in common so that we can map2 over them
vector_chr_in_commmon <- intersect(names(list_recon_gtf_subset_by_chr), names(list_junction_table_by_chr))

list_recon_gtf_subset_by_chr <- list_recon_gtf_subset_by_chr[vector_chr_in_commmon]
list_junction_table_by_chr <- list_junction_table_by_chr[vector_chr_in_commmon]

cat("get positions of the vector where the path of the ref. genome fasta\n")
vector_ref_genome_paths_by_chr_position <- vector_chr_in_commmon %>% purrr::map(.f = function(a1) {
  
  # DEBUG ###
  # a1 <- 1
  ###########
  
  # cat(a1, "\n")
  
  ref_genome_path_by_chr_position <- grep(x = vector_ref_genome_paths_by_chr, pattern = paste("(\\D|^)", a1, ".fa$", sep = ""))
  
  if (length(ref_genome_path_by_chr_position) != 1) {
    
    stop("Something is wrong with the contents of the fasta file directory. Please check that it's structured in the desired format.")
    
  } else {
    
    return(ref_genome_path_by_chr_position)
    
  }
  
} ) %>% unlist

# fetch full path of ref. genome fasta
vector_ref_genome_fasta_path <- paste(vector_ref_genome_paths_by_chr[vector_ref_genome_paths_by_chr_position])

cat("Map over the junction table, recon GTF and ref. genome assembly fasta\n")

list_junction_3FT_result <- future_pmap(
  .l = list("a1" = list_junction_table_by_chr,
            "a2" = list_recon_gtf_subset_by_chr,
            "a3" = vector_ref_genome_fasta_path,
            "a4" = vector_chr_in_commmon),
  .f = function(a1, a2, a3, a4) {
    
    # DEBUG ###
    # a1 <- list_junction_table_by_chr[[1]]
    # a2 <- list_recon_gtf_subset_by_chr[[1]]
    # a3 <- vector_ref_genome_fasta_path[[1]]
    # a4 <- chr_to_run[[1]]
    ###########
    
    cat("temporary allocation to ref genome fasta list\n")
    reference_genome_fasta_chr_temp <- seqinr::read.fasta(file = a3, forceDNAtolower = FALSE)
    
    # Strategy: match junctions to transcripts with directly flanking exons
    # add artificial entry alongside the parent transcript entries indicating the junction
    # later, we use the junction coords to determine the start/end coords according to the specified window using the transcript topology.
    cat("get junction-flanking exon matches from the GTF\n")
    list_GTF_matching_junction_entries <- future_imap(.x = a1 %>% array_tree, .f = function(b1, b2) {
      
      # DEBUG ###
      # b1 <- a1 %>% array_tree %>% .[[1]]
      ###########
      
      cat(b2, "/", nrow(b2), "\n")
      
      matching_GTF_entries <- extract_junction.flanking.exons(query_chr = b1$chr %>% type.convert,
                                                              query_start = b1$start %>% type.convert,
                                                              query_end = b1$end %>% type.convert,
                                                              query_strand = b1$strand %>% type.convert,
                                                              tibble_gtf_table = a2, 
                                                              tolerance_left = 1, 
                                                              tolerance_right = 1, 
                                                              tolerance_inside = 1, 
                                                              tolerance_outside = 0, 
                                                              match_consecutive = TRUE, 
                                                              return_type = "exon")
      
      matching_GTF_entries_with_artificial_entry <- purrr::map(.x = matching_GTF_entries, .f = function(c1) {
        
        # DEBUG ###
        # c1 <- matching_GTF_entries[[1]]
        ###########
        
        # splice in the artificial overlapping exon element using the flanking exons.
        # this entry describes the matched junction in the GTF, which is practically a magnetisation.
        # it will be a single GTF row, with median exon number, spanning the intron junction of the matched transcript.
        list_matched_GTF_entries_with_junction_specification <- purrr::splice(c1,
          "junction_specifications" = c1$matched_flanking_exons %>% .[1, ] %>% dplyr::select(., -start, -end, -width, -exon_number) %>% add_column("start" = (c1$matched_flanking_exons %>% .[1, "end"] %>% paste %>% as.numeric + 1), "end" = (c1$matched_flanking_exons %>% .[2, "start"] %>% paste %>% as.numeric - 1), "exon_number" = mean(c1$matched_flanking_exons %>% .$exon_number %>% as.numeric)) %>% add_column("width" = .$end %>% as.numeric - .$start %>% as.numeric + 1))
        
        # if intronic, then splice in the intronic region into the parent transcript.
        if (grepl(x = b1$splicemode, pattern = intron_retention_string) == TRUE) {
          
          list_matched_GTF_entries_with_junction_specification$parent_transcript <- dplyr::bind_rows(list_matched_GTF_entries_with_junction_specification$parent_transcript, 
                                                                                                     list_matched_GTF_entries_with_junction_specification$junction_specifications) %>% dplyr::arrange(exon_number)
          
        }
        
        return(list_matched_GTF_entries_with_junction_specification)
        
      } ) # L3
      
      final_identifier <- if (b1$custom_identifier %>% is.na != TRUE) {b1$custom_identifier
      } else {
          paste(b1$chr, ":", b1$start %>% type.convert, "-", b1$end%>% type.convert, 
                if ((b1$strand == "+" | b1$strand == "-") & b1$strand %>% is.na != TRUE) {
                  paste(":", b1$strand, sep = "")
                } else {
                    ""
                }, sep = "")
      }
      
      # create fasta header
      fasta_header <- paste(source_tag, 
                            "|", 
                            final_identifier, 
                            "|", 
                            matching_GTF_entries %>% names %>% paste(collapse = ","), 
                            " OS=",
                            b1$organism, 
                            " GN=",
                            b1$gene_name, sep = "")
      
      return(purrr::splice(b1 %>% flatten,
                           "matching_GTF_entries" = matching_GTF_entries_with_artificial_entry %>% list,
                           "final_identifier" = final_identifier %>% list,
                           "fasta_header" = fasta_header %>% list))
      
    }, .progress = TRUE, .options = future_options(globals = c("reconstructed_gtf_temp", "extract_junction.flanking.exons", "splicemode_column_name", "data.table", "type_convert", "dplyr", "tibble", "intron_retention_string", "a2", "a3", "a4", "source_tag")) ) # L2
    
    
    
    cat("get co-ordinates for all nucleotides in the associated transcripts, then add in the forward nucleotide sequence for the whole transcript\n")
    cat("lookup reference FASTA for transcript co-ordinates\n")
    list_matched_coords_temp <- future_imap(.x = list_GTF_matching_junction_entries, .f = function(b1, b2) {
      
      # DEBUG ###
      # b1 <- list_GTF_matching_junction_entries[[269]]
      ###########
      
      cat("now processing entry number", b2, "/", length(list_GTF_matching_junction_entries), "\n")
      
      # extract the magnetised junction start and end coords.
      result <- list("3FT_info" = purrr::map(.x = b1$matching_GTF_entries, .f = ~list(
        "matched_junction_chr" = .x$junction_specifications$seqnames %>% paste,
        "matched_junction_start" = .x$junction_specifications$start %>% paste,
        "matched_junction_end" = .x$junction_specifications$end %>% paste,
        "matched_junction_strand" = .x$junction_specifications$strand %>% paste,
        # generate vector of all nucleotide coors from $start to $end
        "all_parent_transcript_coords" = purrr::map2(.x = .x[["parent_transcript"]]$start, .y = .x[["parent_transcript"]]$end, .f = ~.x:.y) %>% unlist %>% sort) %>%
          # add in the TRANSCRIPT-RELATIVE POSITIONS of the JUNCTION ACCORDING TO THE SPECIFIED WINDOW
          purrr::splice(., 
                        "translation_window_start_transcript.relative" = if (.$matched_junction_strand == "+") {
                          max(1, which(.$all_parent_transcript_coords == (.$matched_junction_start %>% as.numeric - 1)) - upstream_window_size + 1) 
                        } else if (.$matched_junction_strand == "-") {
                          max(1, which(.$all_parent_transcript_coords %>% rev == (.$matched_junction_end %>% as.numeric + 1)) - upstream_window_size + 1) 
                        }, 
                        "translation_window_end_transcript.relative" = if (.$matched_junction_strand == "+") {
                          min(length(.$all_parent_transcript_coords), which(.$all_parent_transcript_coords == (.$matched_junction_end %>% as.numeric + 1)) + downstream_window_size - 1) 
                        } else if (.$matched_junction_strand == "-") {
                          min(length(.$all_parent_transcript_coords), which(.$all_parent_transcript_coords %>% rev == (.$matched_junction_start %>% as.numeric - 1)) + downstream_window_size - 1) 
                        },
                        "last_nt_before_splice_junction_transcript_relative" = if (.$matched_junction_strand == "+") {
                          which(.$all_parent_transcript_coords == (.$matched_junction_start %>% as.numeric - 1))
                        } else if (.$matched_junction_strand == "-") {
                          which(.$all_parent_transcript_coords %>% rev == (.$matched_junction_end %>% as.numeric + 1))
                        },
                        "first_nt_after_splice_junction_transcript_relative" = if (.$matched_junction_strand == "+") {
                          which(.$all_parent_transcript_coords == (.$matched_junction_end %>% as.numeric + 1))
                        } else if (.$matched_junction_strand == "-") {
                          which(.$all_parent_transcript_coords %>% rev == (.$matched_junction_start %>% as.numeric - 1)) 
                        },
                        # add in all the nucleotide coords of the PARENT transcript
                        "parent_transcript_forward_nucleotides" = reference_genome_fasta_chr_temp[[b1$chr %>% paste]][.$all_parent_transcript_coords])))
      
      return(splice(b1,
                    result))
      
    }, .progress = TRUE)
    
    
    cat("doing three-frame translate using the fetched co-ordinates\n")
    list_3FT_result_temp <- future_imap(.x = list_matched_coords_temp, .f = function(b1, b2) {
      
      # DEBUG ###
      # b1 <- list_matched_coords_temp[[269]]
      ###########
      
      cat("now processing entry number", b2, "/", length(list_matched_coords_temp), "\n")
      
      updated_3FT_info <- purrr::imap(.x = b1$`3FT_info`, .f = function(c1, c2) {
        
        # DEBUG ###
        # c1 <- b1$`3FT_info`$MSTRG.241.10
        ###########
        
        # cat("now processing entry number", c2, "/", length(list_matched_coords_temp), "\n")
        
        # three frame translation, add translation frame-relative coordinates of window start and end
        purrr::splice(c1, nt.sequence_strand_threeframetranslate(vector_forward_nucleotides = c1$parent_transcript_forward_nucleotides, strand = c1$matched_junction_strand),
                      "window_start_AA_position_frame_0" = calculate_translation_frame_relative_start_end_position(ES = c1$translation_window_start_transcript.relative, EE = c1$translation_window_end_transcript.relative, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 0)$exon_start_AA_position,
                      "window_start_AA_position_frame_1" = calculate_translation_frame_relative_start_end_position(ES = c1$translation_window_start_transcript.relative, EE = c1$translation_window_end_transcript.relative, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 1)$exon_start_AA_position,
                      "window_start_AA_position_frame_2" = calculate_translation_frame_relative_start_end_position(ES = c1$translation_window_start_transcript.relative, EE = c1$translation_window_end_transcript.relative, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 2)$exon_start_AA_position,
                      "window_end_AA_position_frame_0" = calculate_translation_frame_relative_start_end_position(ES = c1$translation_window_start_transcript.relative, EE = c1$translation_window_end_transcript.relative, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 0)$exon_end_AA_position,
                      "window_end_AA_position_frame_1" = calculate_translation_frame_relative_start_end_position(ES = c1$translation_window_start_transcript.relative, EE = c1$translation_window_end_transcript.relative, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 1)$exon_end_AA_position,
                      "window_end_AA_position_frame_2" = calculate_translation_frame_relative_start_end_position(ES = c1$translation_window_start_transcript.relative, EE = c1$translation_window_end_transcript.relative, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 2)$exon_end_AA_position,
                      "junction_AA_position_frame_0_upstream" = calculate_translation_frame_relative_start_end_position(ES = c1$last_nt_before_splice_junction_transcript_relative + 1, EE = c1$last_nt_before_splice_junction_transcript_relative + 1, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 0)$exon_start_AA_position,
                      "junction_AA_position_frame_1_upstream" = calculate_translation_frame_relative_start_end_position(ES = c1$last_nt_before_splice_junction_transcript_relative + 1, EE = c1$last_nt_before_splice_junction_transcript_relative + 1, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 1)$exon_start_AA_position,
                      "junction_AA_position_frame_2_upstream" = calculate_translation_frame_relative_start_end_position(ES = c1$last_nt_before_splice_junction_transcript_relative + 1, EE = c1$last_nt_before_splice_junction_transcript_relative + 1, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 2)$exon_start_AA_position,
                      "junction_AA_position_frame_0_downstream" = calculate_translation_frame_relative_start_end_position(ES = c1$first_nt_after_splice_junction_transcript_relative - 1, EE = c1$first_nt_after_splice_junction_transcript_relative - 1, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 0)$exon_start_AA_position,
                      "junction_AA_position_frame_1_downstream" = calculate_translation_frame_relative_start_end_position(ES = c1$first_nt_after_splice_junction_transcript_relative - 1, EE = c1$first_nt_after_splice_junction_transcript_relative - 1, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 1)$exon_start_AA_position,
                      "junction_AA_position_frame_2_downstream" = calculate_translation_frame_relative_start_end_position(ES = c1$first_nt_after_splice_junction_transcript_relative - 1, EE = c1$first_nt_after_splice_junction_transcript_relative - 1, TL = length(c1$all_parent_transcript_coords), strand = "+", frame = 2)$exon_start_AA_position) %>%
          # check for upstream ORF as well as downstream ORF (starting from within the exon) by seeing if you reverse the AA sequence, do you see a methionine always before the first stop codon
          # will enter logical indicating whether the exon has a valid uORF or not in the given translation frame
          purrr::splice(
            # uORF
            "uORF_valid_frame_0" = find_valid_uORF(.[c("translation_frame_0", "window_start_AA_position_frame_0", "window_end_AA_position_frame_0", "junction_AA_position_frame_0_upstream")]),
            "uORF_valid_frame_1" = find_valid_uORF(.[c("translation_frame_1", "window_start_AA_position_frame_1", "window_end_AA_position_frame_1", "junction_AA_position_frame_1_upstream")]),
            "uORF_valid_frame_2" = find_valid_uORF(.[c("translation_frame_2", "window_start_AA_position_frame_2", "window_end_AA_position_frame_2", "junction_AA_position_frame_2_upstream")]),
            # dORF
            "dORF_valid_frame_0" = find_valid_dORF(.[c("translation_frame_0", "window_start_AA_position_frame_0", "window_end_AA_position_frame_0", "junction_AA_position_frame_0_downstream")]),
            "dORF_valid_frame_1" = find_valid_dORF(.[c("translation_frame_1", "window_start_AA_position_frame_1", "window_end_AA_position_frame_1", "junction_AA_position_frame_1_downstream")]),
            "dORF_valid_frame_2" = find_valid_dORF(.[c("translation_frame_2", "window_start_AA_position_frame_2", "window_end_AA_position_frame_2", "junction_AA_position_frame_2_downstream")])
          ) %>%
          # collapse the forward nucleotide and 3FT sequences into string from a vector
          purrr::modify_at(.x = ., .at = c("parent_transcript_forward_nucleotides", "translation_frame_0", "translation_frame_1", "translation_frame_2"), .f = ~.x %>% paste(collapse = "")) %>% 
          # collapse the parent transcript coords into string from a vector
          purrr::modify_at(.x = ., .at = c("all_parent_transcript_coords"), .f = ~.x %>% paste(collapse = ",")) %>% 
          
          return
        
      } ) # L3
      
      updated_list <- b1
      
      updated_list$`3FT_info` <- updated_3FT_info
      
      return(updated_list)
      
    }, .progress = TRUE) # L2
    
    # flatten and distribute the fasta header into all of its child 3FT results
    list_3FT_result_unnest_temp <- list_3FT_result_temp %>% purrr::discard(.p = ~.x %>% length == 0) %>% purrr::discard(.p = ~.x$`3FT_info` %>% length == 0) %>% purrr::map(.x = list_3FT_result_temp, .f = ~purrr::cross2(.x[c("chr", "start", "end", "strand", "fasta_header", "final_identifier")] %>% list, .x$`3FT_info`)) %>% flatten %>% purrr::map(.f = ~.x %>% flatten)
    
    # rearrange info into a tibble. one translation frame per row. %>% flatten
    # first, tibblise the translation frame elements of the list, then as_tibble the rest.
    cat("rearrange the three-frame translate results into a tibble\n")
    list_3FT_result_unnest_temp_2 <- future_imap(.x = list_3FT_result_unnest_temp, .f = function(b1, b2) {
      
      # DEBUG #####
      # b1 <- list_3FT_result_unnest_temp[[1]]
      #############
      
      cat("now processing entry number", b2, "/", length(list_3FT_result_unnest_temp), "\n")
      
      # define the element indices inside the list which contains info for each frame.
      element.indices_frame_info <- grep(x = names(b1), pattern = "frame_\\d$")
      # element indices containing the virtual peptides
      element.indices_virtual_peptides <- grep(x = names(b1), pattern = "translation_frame_\\d")
      # element indices of the frame info but not the virtual peptides
      element.indices_not_virtual_peptides <- setdiff(element.indices_frame_info, element.indices_virtual_peptides)
      
      # make the tibble of the virtual peptides first.
      tibble_frame_info <- b1[c("translation_frame_0", "translation_frame_1", "translation_frame_2")] %>% as.data.frame %>% t %>% as_tibble(rownames = "translation_frame", .name_repair = "unique")
      colnames(tibble_frame_info)[2] <- "parent_transcript_virtual_peptide_sequence"
      # tidy the frame number
      tibble_frame_info[, "translation_frame"] <- gsub(x = tibble_frame_info$translation_frame, pattern = "translation_frame_", replacement = "")
      
      # then make the tibble of the rest. join it onto the existing table
      tibble_frame_info_temp <- b1[element.indices_not_virtual_peptides] %>% as.data.frame %>% t %>% matrix(nrow = 3, byrow = FALSE, dimnames = list(0:2, gsub(x = names(b1[element.indices_not_virtual_peptides]), pattern = "(_frame_\\d)", replacement = "") %>% unique)) %>% as_tibble(rownames = "translation_frame", .name_repair = "unique")
      # table join
      tibble_frame_info <- dplyr::full_join(tibble_frame_info, tibble_frame_info_temp, by = "translation_frame")
      
      # finally join the tibblised parts of the list to the non-tibblised.
      list_half_tibblised <- splice("frame_info" = tibble_frame_info, b1[-element.indices_frame_info])
      
      fully_tibblised_list <- list_half_tibblised %>% flatten %>% as_tibble
      
      return(fully_tibblised_list)
      
    }, .progress = TRUE) # L2
    
    return(list_3FT_result_unnest_temp_2)
    
  }, .progress = TRUE ) # L1

# rbind into a summary tibble
tibble_three_frame_translate_result <- list_junction_3FT_result %>% flatten %>% rbindlist %>% as_tibble

# remove entries with no valid translations. 
tibble_three_frame_translate_result <- tibble_three_frame_translate_result[-which((tibble_three_frame_translate_result$uORF_valid == "NONE_VALID" | tibble_three_frame_translate_result$uORF_valid == "") & (tibble_three_frame_translate_result$dORF_valid == "NONE_VALID" | tibble_three_frame_translate_result$dORF_valid == "")), ]

# remove all rows with duplicated u/dORF_valid columns from the table
tibble_three_frame_translate_result <- tibble_three_frame_translate_result[-which(duplicated(tibble_three_frame_translate_result[, c("uORF_valid", "dORF_valid")])), ]

# bind the uORF and dORF translation results together.
tibble_three_frame_translate_result_unique.entries <- dplyr::bind_rows(tibble_three_frame_translate_result %>% dplyr::select(-dORF_valid) %>% dplyr::rename("virtual_peptide_sequence" = "uORF_valid"), tibble_three_frame_translate_result %>% dplyr::select(-uORF_valid) %>% dplyr::rename("virtual_peptide_sequence" = "dORF_valid")) %>% 
  # remove all "NONE_VALID" sequences
  .[-which(.$virtual_peptide_sequence == "NONE_VALID" | .$virtual_peptide_sequence == ""), ] %>% dplyr::distinct(., virtual_peptide_sequence, .keep_all = TRUE) 

tibble_three_frame_translate_result_unique.entries <- tibble_three_frame_translate_result_unique.entries %>% 
  # add column to indicate if the virtual peptide per row is a substring of another row
  add_column("substring_or_not" = future_imap(.x = tibble_three_frame_translate_result_unique.entries$virtual_peptide_sequence, .f = ~grepl(x = tibble_three_frame_translate_result_unique.entries$virtual_peptide_sequence %>% .[-.y], pattern = .x) %>% any == TRUE, .progress = TRUE) %>% unlist)

# filter out substrings
tibble_three_frame_translate_result_no_substring <- tibble_three_frame_translate_result_unique.entries %>% dplyr::filter(substring_or_not == FALSE)

# tally up the number of valid frames we ended up with
tibble_junctions_frame_tally <- tibble_three_frame_translate_result_no_substring %>% dplyr::distinct(translation_frame, fasta_header) %>% dplyr::group_by(fasta_header) %>% dplyr::summarise("tally" = n())

cat("\nnumber of junctions input: ", tibble_junction_table %>% dplyr::distinct(chr, start, end, strand) %>% nrow, "\n")
cat("\nnumber of junctions translated: ", tibble_three_frame_translate_result_no_substring$final_identifier %>% unique %>% length, "\n")
cat("\naverage number of translation frames for UNIQUE junction sequences: ", mean(tibble_junctions_frame_tally$tally), "\n")

# WE WRITE A TIBBLE CONTAINING THE GENOME COORD-PEPTIDE MAPPING
# write a table
write.table(x = tibble_three_frame_translate_result_no_substring, file = paste(output_dir, "/", output_name, "_3FT.summary.info.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# FINALLY! WE WRITE THE FASTA!
write.fasta(sequences = tibble_three_frame_translate_result_no_substring$virtual_peptide_sequence %>% array_tree %>% flatten, names = tibble_three_frame_translate_result_no_substring$fasta_header, file.out = paste(output_dir, "/", output_name, ".fasta", sep = ""), open = "w", nbchar = 40, as.string = TRUE)

# write final junction table as .bed file
cat(paste("track name=\"Alternative junctions\" description=\"", output_name, "\" graphType=junctions\n", sep = ""), file = paste(output_dir, "/", output_name, "_junctions.bed", sep = ""))

junc_bed_table <- tibble_three_frame_translate_result_no_substring[, c("chr", "start", "end", "final_identifier", "strand")] %>% setNames(c("chr", "start", "end", "name", "strand")) %>% add_column(., "score" = 1000, .after = "name") %>% type_convert

write.table(junc_bed_table, file = paste(output_dir, "/", output_name, "_junctions.bed", sep = ""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

# finish counting
tictoc::toc()
  
q()

