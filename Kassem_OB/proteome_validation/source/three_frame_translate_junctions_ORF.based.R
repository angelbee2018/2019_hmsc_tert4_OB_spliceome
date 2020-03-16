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
  "1" = make_option(c("-F", "--file_information_path"), type = "character", default = NULL, 
              help = "Compulsory. path to a table containing the junction, reconstructed GTF file and desired output name information. should be in the format of a tab separated tibble to minimise the risk of error. each row describes A PATH the file to the junction table and reconstr. GTF to be translated. This program will automatically loop thru them. The column names must be: c(\"reconstructed_GTF_path\", \"junction_table_path\" and \"output_database_name\").
              - junction_table_path: path to table containing junctions of interest. for example, UNION_junc_coor type files. MUST contain columns start, end, chr, strand, and fasta_header. one row per junction. should be in the form of a tibble. the fasta header you want included in the final FASTA file needs to be in the column called \"fasta_header\".
              - reconstructed_GTF_path: path to the actual reconstructed GTF file (e.g. from Cufflinks, Strawberry). NOT THE CONTAINING DIRECTORY. tip: for better junction matching, combine the reconstructed GTF with reference GTF beforehand e.g. using StringTie
              - output_database_name: a character string of what the final FASTA database file name and output table will be called.", metavar = "character"),
  "2" = make_option(c("-R", "--reference_genome_fasta_dir"), type = "character", default = NULL, 
              help = "Compulsory. path to the directory containing the genome FASTA files. Ideally from Ensembl... you need separate files by chromosomes, NOT the primary assembly. 
              FORMATTING IMPORTANT!!!! MAKE SURE THE REF. GENOME FASTA FILES ARE IN THE FORMAT: <_anything_><chr>.fa e.g. \"Homo_sapiens.GRCh38.dna.chromosome.MT.fa\" OR \"chr16.fa\" OR \"Y.fa\". What will not work: anything which does not have .fa extension e.g. \"chr16.fasta\", anything between the chromosome number and the .fa extension e.g. \"chromosome1.ensembl.fa\"", metavar = "character"),
  "3" = make_option(c("-U", "--upstream_window_size"), type = "double", default = 50,
              help = "Optional. how many nucleotides to translate upstream W.R.T. the middle of the exon for junction-based mode, the total nt. length to be translated will be arg3 + arg4. default for both is 50.", metavar = "double"),
  "4" = make_option(c("-D", "--downstream_window_size"), type = "double", default = 50,
  help = "Optional. how many nt to translate downstream W.R.T. the transcript. for exon-based mode, both will be 50 by default. instead, the window size is taken as the MINIMUM required length to be considered for translation starting from the middle of the exon.", metavar = "double"),
  "5" = make_option(c("-O", "--output_dir"), type = "character", default = NULL, 
              help = "Compulsory. output directory. where do you want to save the custom databases? IMPORTANT: directory must end in a \"/\". e.g. correct: ~/outputdir/ incorrect: ~/outputdir", metavar = "character"),
  "6" = make_option(c("-C", "--ncores"), type = "integer", default = 0, 
                  help = "Optional. Number of cores to use. possible inputs: numbers 1 to any integer. By default, uses all cores (ncores = 0).", metavar = "integer"),
  "7" = make_option(c("-H", "--chrmode"), type = "integer", default = 0, 
                  help = "Optional. Specifies which chromosomes to do: select what chromosomes you want translated. possible inputs: numbers 0-2. 0 (default): nuclear chromosomes only i,e, 1:22, X & Y. 1: nuclear + mitochondrial i.e. 1:22, X & Y, M. 2: everything including haplotype/fusion chromosomes etc... this is possible provided the chromosome names.", metavar = "integer"),
  "8" = make_option(c("-N", "--nonchrname"), type = "character", default = NULL, 
                    help = "Compulsory only if you have specified \"--chrmode 2\". nonchromosomal file name. if you are doing haplotypes, please specify what the reference genome FASTA file for it is called or the script won't know. This single FASTA file must contain all the haplotype information. The script won't try to search for a second file. In ensembl, this file is called \"Homo_sapiens.GRCh38.dna.nonchromosomal.fa\" or more generally, \"*nonchromosomal.fa\". So for this option, you would specify \"--nonchrname nonchromosomal\".", metavar = "character")
)

input_arg_info <- OptionParser(option_list = list_input_arg_info, description = script_description)
input_args <- input_arg_info %>% parse_args

# check if the input arguments are O.K
if ((list(input_args$file_information_path, input_args$reference_genome_fasta_dir, input_args$output_dir) %>% lapply(is.null) %>% unlist %>% any == TRUE) | 
    (input_args$chrmode == 2 & is.null(input_args$nonchrname) == TRUE)) {
  
  print_help(input_arg_info)
  
  stop("Make sure you entered the arguments correctly", call. = FALSE)
  
}

file_information_path <- input_args$file_information_path
file_information_table <- read_delim(file_information_path, delim = "\t")
# CHECK IF THE FILE INFO TABLE PROVIDED CONTAINS ALL THE COLUMNS WE WANT. STOP IF WE DONT SEE ANY OF THE COLNAMES WE WANT TO FIND.
if (
  any(str_detect(string = colnames(file_information_table), pattern = c("junction_table_path", "reconstructed_GTF_path", "output_database_name")) == FALSE)
    ) {
  
  stop("Error. Please check that the file information table is formatted properly. Use --help command for more info.", call. = FALSE)
  
}

reference_genome_fasta_dir <- input_args$reference_genome_fasta_dir
upstream_window_size <- input_args$upstream_window_size
downstream_window_size <- input_args$downstream_window_size
output_dir <- input_args$output_dir

if(!dir.exists(output_dir) ) {
  dir.create(output_dir, recursive = TRUE)}

# manage parrallellisation rrlllRll

if (input_args$ncores != 0) {
  number_of_workers <- input_args$ncores
} else {
  number_of_workers <- future::availableCores()
} 

future::plan(multisession)
options(future.globals.maxSize = 30000000000, mc.cores = number_of_workers)

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

cat(future::availableCores(), "cores will be used\n")

# testing ########

# file_information_path <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/input_files/file_info_test.txt"
# reconstructed_gtf_path <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/GRAND_OBseries_ref_reconstructed_stringtiemerged.gtf"
# reference_genome_fasta_path <- "/srv/scratch/z3463471/hg38_ensembl_reference/raw_genome_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# output_dir <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/results_database_generation/"
# output_database_name <- "three_frame_translation_text"

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
# index: loop progress marker to be used with imap

extract_junction.flanking.exons <- function(spliceregion_list, tibble_gtf_table, index) {
  
  # DEBUG ###################
  
  # index <- 1
  # spliceregion_list <- wide_tibble_of_all_unique_exon_coords_array.tree[[3]]
  # tibble_gtf_table <- tibble_ref_gtf
  # tibble_gtf_table <- tibble_alltimepoints_reconstructed_gtf
  
  ###########################
  
  print(paste("now processing junction number", index))
  
  if (spliceregion_list$strand %>% trimws == ".") {
    
    # +/- 1 nt tolerance
    tibble_gtf_subset_flanking_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == spliceregion_list$chr %>% trimws, ] %>% 
      .[.$start <= ((spliceregion_list$end %>% as.numeric) + 2) & .$end >= ((spliceregion_list$start %>% as.numeric) - 2), ] %>% 
      .[!(.$start <= ((spliceregion_list$end %>% as.numeric) - 1) & .$end >= ((spliceregion_list$start %>% as.numeric) + 1)), ] %>% 
      .[.$type == "exon", ]
    
  } else if (spliceregion_list$strand %>% trimws == "+" | spliceregion_list$strand %>% trimws == "-") {
    
    tibble_gtf_subset_flanking_exons <- tibble_gtf_table[(tibble_gtf_table$seqnames == spliceregion_list$chr %>% trimws) & (tibble_gtf_table$strand == spliceregion_list$strand %>% trimws), ] %>% 
      .[.$start <= ((spliceregion_list$end %>% as.numeric) + 2) & .$end >= ((spliceregion_list$start %>% as.numeric) - 2), ] %>% 
      .[!(.$start <= ((spliceregion_list$end %>% as.numeric) - 1) & .$end >= ((spliceregion_list$start %>% as.numeric) + 1)), ] %>% 
      .[.$type == "exon", ]
    
  } else {
    
    stop("Could not match the strand information in the transposed differential-only UNION_junc_coor_table. Make sure that the \"strand\" column in the UNION_junc_coor_table contains only +, - or .")
    
  }
  
  list_of_junction_associated_transcripts <- tibble_gtf_subset_flanking_exons$transcript_id %>% unique %>% array_tree %>% flatten
  
  # make a list for each transcript that directly flanks the junction.
  # then filter so that there are only a) exon PAIRS which b) are directly connected in the mature (spliced) transcript
  list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_gtf_subset_flanking_exons[tibble_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts) %>% keep(.x = ., .p = ~nrow(.x) == 2) %>% keep(.x = ., .p = ~abs((.x[2, "exon_number"] %>% paste %>% as.numeric) - (.x[1, "exon_number"] %>% paste %>% as.numeric)) == 1)
  
  # for those consecutive exons which were found to flank a junction, get all the entries of the parent transcript
  list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_tibbles_flanking_exon_gtf.entries_per_transcript, .f = ~list(
    
    "matched_flanking_exons" = .x, 
    "parent_transcript" = tibble_gtf_table[tibble_gtf_table$transcript_id == .x$transcript_id %>% unique %>% paste, ] %>% .[-which(is.na(.$exon_number)), ] %>%
      dplyr::arrange(exon_number %>% as.numeric)))
  
  return(list_of_tibbles_flanking_exon_gtf.entries_per_transcript)
  
}

# END extract_junction.flanking.exons ########

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

##############################################################################

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

# END find_valid_uORF() and find_valid_dORF() #######################################################################

# MAIN FUNCTION TO DO 3FT OF JUNCTIONS
# Behaviour: for each junction, look up every single transcript that it's associated with. Translate in three frames and find which frame can validly translate the exon.
# These transcripts are matched according to both reference and reconstructed GTF.
# The smallest protein in humans is 44 AA so the smallest valid translatable nucleotide length is 132 nt. 
# Translated region included in the upstream/downstream window must cross the splice junction.

# END FUNCTIONS #####################################################################################################
#####################################################################################################################
#####################################################################################################################

cat("# BEGIN EXECUTION #################################\n")

# DEBUG ################

# reconstructed_gtf_path <- "Z:/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/GRAND_OBseries_ref_denovo_reconstructed_stringtiemerged.gtf"
# reconstructed_gtf <- rtracklayer::import(reconstructed_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
# 
# reference_genome_fasta_dir <- "Z:/hg38_ensembl_reference/raw_genome_fasta/genome_fasta_extract2/"
# 
# junction_table_path <- "Z:/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/results_database_generation/angel_3FT_junctions/junction_table_OBseries_SOM_1663_junctions_any_qvalue0.01_any_deltaPSI_greaterthan_0.2.txt"
# 
# tibble_junction_table <- read.delim(junction_table_path, sep = "\t", stringsAsFactors = FALSE) %>% as_tibble
# upstream_window_size <- 50
# downstream_window_size <- 50

########################


vector_ref_genome_paths_by_chr <- paste(reference_genome_fasta_dir, list.files(reference_genome_fasta_dir)[list.files(reference_genome_fasta_dir) %>% grep(., pattern = ".*.fa$")], sep = "")

# reference_genome_fasta <- seqinr::read.fasta(file = reference_genome_fasta_path, forceDNAtolower = FALSE)

# total_peptide_window_size <- upstream_window_size + downstream_window_size
# cat("total translation window size:", total_peptide_window_size, "\n")

## loop thru reconstructed GTF path slowly, loop thru junction files quickly.
## the reason why we loop is to save importing the reference_genome_fasta many times for many combinations of reconstructed GTF and junction files.

for (i in 1:nrow(file_information_table)) {
  
  reconstructed_gtf_path <- file_information_table[i, "reconstructed_GTF_path"] %>% paste
  cat("reconstructed_gtf_path:", reconstructed_gtf_path, "\n")
  junction_table_path <- file_information_table[i, "junction_table_path"] %>% paste
  cat("junction_table_path:", junction_table_path, "\n")
  output_database_name <- file_information_table[i, "output_database_name"] %>% paste
  cat("output_database_name:", output_database_name, "\n")
    
  reconstructed_gtf <- rtracklayer::import(reconstructed_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)

  cat("GTF importing done\n")

  tibble_junction_table <- read.delim(junction_table_path, sep = "\t", stringsAsFactors = FALSE) %>% as_tibble

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
  
  # start counting execution time of the 3FT portion
  tictoc::tic("3FT execution time")
  
  # begin looping thru each chromosome
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
    
    cat("temporary subsetting of reconstructed GTF to reduce overhead when exporting\n")
    reconstructed_gtf_temp <- reconstructed_gtf[reconstructed_gtf$seqnames == chr, ]
    
    cat("prepare to loop thru each junction entry.\n")
    tibble_junction_table_temp_array.tree <- tibble_junction_table[tibble_junction_table$chr == chr, ] %>% array_tree
    
    # Strategy: match junctions to transcripts with directly flanking exons
    # add artificial entry alongside the parent transcript entries indicating the junction
    # later, we use the junction coords to determine the start/end coords according to the specified window using the transcript topology.
    cat("get junction-flanking exon matches from the GTF\n")
    list_GTF_matching_junction_entries <- future_imap(.x = tibble_junction_table_temp_array.tree, .f = ~list(
      # "chr" = .x$chr %>% trimws,
      # "start" = .x$start %>% trimws,
      # "end" = .x$end %>% trimws,
      # "strand" = .x$strand %>% trimws,
      # "splicemode" = .x[[splicemode_column_name]] %>% trimws,
      "fasta_header" = .$fasta_header,
      "matching_GTF_entries" = extract_junction.flanking.exons(.x, reconstructed_gtf_temp, .y) %>% 
        # add artificial overlapping exon element using the flanking exons.
        # it will be a single GTF row, with median exon number, spanning the intron junction of the matched transcript.
        purrr::map(.x = ., .f = ~.x %>% purrr::splice(
          "junction_specifications" = .x$matched_flanking_exons %>% .[1, ] %>% dplyr::select(., -start, -end, -width, -exon_number) %>% add_column("start" = (.x$matched_flanking_exons %>% .[1, "end"] %>% paste %>% as.numeric + 1), "end" = (.x$matched_flanking_exons %>% .[2, "start"] %>% paste %>% as.numeric - 1), "exon_number" = mean(.x$matched_flanking_exons %>% .$exon_number %>% as.numeric)) %>% add_column("width" = .$end %>% as.numeric - .$start %>% as.numeric + 1)
        ))),
      .progress = TRUE, .options = future_options(globals = c("reconstructed_gtf_temp", "extract_junction.flanking.exons", "splicemode_column_name", "data.table", "type_convert", "dplyr", "tibble")))
    
    cat("get co-ordinates for all nucleotides in the associated transcripts, then add in the forward nucleotide sequence for the whole transcript\n")
    cat("lookup reference FASTA for transcript co-ordinates\n")
    list_matched_coords_temp <- purrr::imap(.x = list_GTF_matching_junction_entries, .f = function(.x, .y) {
      
      cat("now processing entry number", .y, "/", length(list_GTF_matching_junction_entries), "\n")
      
      result <- list("fasta_header" = .x[["fasta_header"]], "3FT_info" = purrr::map(.x = .x[["matching_GTF_entries"]], .f = ~list(
        "matched_junction_chr" = chr,
        "matched_junction_start" = .x[["junction_specifications"]]$start %>% paste,
        "matched_junction_end" = .x[["junction_specifications"]]$end %>% paste,
        "matched_junction_strand" = .x[["junction_specifications"]]$strand %>% paste,
        # generate vector of all nucleotide coors from $start to $end
        "all_parent_transcript_coords" = purrr::map2(.x = .x[["parent_transcript"]]$start, .y = .x[["parent_transcript"]]$end, .f = ~.x:.y) %>% unlist %>% sort) %>%
          # add in the TRANSCRIPT-RELATIVE POSITIONS of the JUNCTION ACCORDING TO THE SPECIFIED WINDOW
          purrr::splice(., "translation_window_start_transcript.relative" = max(1, which(.$all_parent_transcript_coords == (.$matched_junction_start %>% as.numeric - 1)) - upstream_window_size + 1), 
                        "translation_window_end_transcript.relative" = min(length(.$all_parent_transcript_coords), which(.$all_parent_transcript_coords == (.$matched_junction_end %>% as.numeric + 1)) + downstream_window_size - 1),
                        # add in all the nucleotide coords of the PARENT transcript
                        "parent_transcript_forward_nucleotides" = reference_genome_fasta_chr_temp[[chr %>% paste]][.$all_parent_transcript_coords])))
      
      return(result)
      
    })
    
    cat("doing three-frame translate using the fetched co-ordinates\n")
    list_3FT_result_temp <- future_map(.x = list_matched_coords_temp, .f = ~list("fasta_header" = .x[["fasta_header"]], "3FT_info" = purrr::map(.x = .x[["3FT_info"]], .f = ~.x %>%
                                                                                                                                                                                                        # three frame translation, add translation frame-relative coordinates of window start and end
                                                                                                                                                                                                        purrr::splice(., nt.sequence_strand_threeframetranslate(vector_forward_nucleotides = .$parent_transcript_forward_nucleotides, strand = .$matched_junction_strand),
                                                                                                                                                                                                                      "window_start_AA_position_frame_0" = calculate_translation_frame_relative_start_end_position(ES = .$translation_window_start_transcript.relative, EE = .$translation_window_end_transcript.relative, TL = length(.$all_parent_transcript_coords), strand = .$matched_junction_strand, frame = 0)$exon_start_AA_position,
                                                                                                                                                                                                                      "window_start_AA_position_frame_1" = calculate_translation_frame_relative_start_end_position(ES = .$translation_window_start_transcript.relative, EE = .$translation_window_end_transcript.relative, TL = length(.$all_parent_transcript_coords), strand = .$matched_junction_strand, frame = 1)$exon_start_AA_position,
                                                                                                                                                                                                                      "window_start_AA_position_frame_2" = calculate_translation_frame_relative_start_end_position(ES = .$translation_window_start_transcript.relative, EE = .$translation_window_end_transcript.relative, TL = length(.$all_parent_transcript_coords), strand = .$matched_junction_strand, frame = 2)$exon_start_AA_position,
                                                                                                                                                                                                                      "window_end_AA_position_frame_0" = calculate_translation_frame_relative_start_end_position(ES = .$translation_window_start_transcript.relative, EE = .$translation_window_end_transcript.relative, TL = length(.$all_parent_transcript_coords), strand = .$matched_junction_strand, frame = 0)$exon_end_AA_position,
                                                                                                                                                                                                                      "window_end_AA_position_frame_1" = calculate_translation_frame_relative_start_end_position(ES = .$translation_window_start_transcript.relative, EE = .$translation_window_end_transcript.relative, TL = length(.$all_parent_transcript_coords), strand = .$matched_junction_strand, frame = 1)$exon_end_AA_position,
                                                                                                                                                                                                                      "window_end_AA_position_frame_2" = calculate_translation_frame_relative_start_end_position(ES = .$translation_window_start_transcript.relative, EE = .$translation_window_end_transcript.relative, TL = length(.$all_parent_transcript_coords), strand = .$matched_junction_strand, frame = 2)$exon_end_AA_position,
                                                                                                                                                                                                                      "junction_AA_position_frame_0" = calculate_translation_frame_relative_start_end_position(ES = mean(c(.$translation_window_start_transcript.relative, .$translation_window_end_transcript.relative)), EE = mean(c(.$translation_window_start_transcript.relative, .$translation_window_end_transcript.relative)), TL = length(.$all_parent_transcript_coords), strand = .$matched_junction_strand, frame = 1)$exon_start_AA_position,
                                                                                                                                                                                                                      "junction_AA_position_frame_1" = calculate_translation_frame_relative_start_end_position(ES = mean(c(.$translation_window_start_transcript.relative, .$translation_window_end_transcript.relative)), EE = mean(c(.$translation_window_start_transcript.relative, .$translation_window_end_transcript.relative)), TL = length(.$all_parent_transcript_coords), strand = .$matched_junction_strand, frame = 1)$exon_end_AA_position,
                                                                                                                                                                                                                      "junction_AA_position_frame_2" = calculate_translation_frame_relative_start_end_position(ES = mean(c(.$translation_window_start_transcript.relative, .$translation_window_end_transcript.relative)), EE = mean(c(.$translation_window_start_transcript.relative, .$translation_window_end_transcript.relative)), TL = length(.$all_parent_transcript_coords), strand = .$matched_junction_strand, frame = 1)$exon_end_AA_position) %>%
                                                                                                                                                                                                        # check for upstream ORF as well as downstream ORF (starting from within the exon) by seeing if you reverse the AA sequence, do you see a methionine always before the first stop codon
                                                                                                                                                                                                        # will enter logical indicating whether the exon has a valid uORF or not in the given translation frame
                                                                                                                                                                                                        purrr::splice(
                                                                                                                                                                                                          # uORF
                                                                                                                                                                                                          "uORF_valid_frame_0" = find_valid_uORF(.[c("translation_frame_0", "window_start_AA_position_frame_0", "window_end_AA_position_frame_0", "junction_AA_position_frame_0")]),
                                                                                                                                                                                                          "uORF_valid_frame_1" = find_valid_uORF(.[c("translation_frame_0", "window_start_AA_position_frame_1", "window_end_AA_position_frame_1", "junction_AA_position_frame_1")]),
                                                                                                                                                                                                          "uORF_valid_frame_2" = find_valid_uORF(.[c("translation_frame_0", "window_start_AA_position_frame_2", "window_end_AA_position_frame_2", "junction_AA_position_frame_2")]),
                                                                                                                                                                                                          # dORF
                                                                                                                                                                                                          "dORF_valid_frame_0" = find_valid_dORF(.[c("translation_frame_0", "window_start_AA_position_frame_0", "window_end_AA_position_frame_0", "junction_AA_position_frame_0")]),
                                                                                                                                                                                                          "dORF_valid_frame_1" = find_valid_dORF(.[c("translation_frame_1", "window_start_AA_position_frame_1", "window_end_AA_position_frame_1", "junction_AA_position_frame_1")]),
                                                                                                                                                                                                          "dORF_valid_frame_2" = find_valid_dORF(.[c("translation_frame_2", "window_start_AA_position_frame_2", "window_end_AA_position_frame_2", "junction_AA_position_frame_2")])
                                                                                                                                                                                                        ) %>%
                                                                                                                                                                                                        # collapse the forward nucleotide and 3FT sequences into string from a vector
                                                                                                                                                                                                        purrr::modify_at(.x = ., .at = c("parent_transcript_forward_nucleotides", "translation_frame_0", "translation_frame_1", "translation_frame_2"), .f = ~.x %>% paste(collapse = "")) %>% 
                                                                                                                                                                                                        # collapse the parent transcript coords into string from a vector
                                                                                                                                                                                                        purrr::modify_at(.x = ., .at = c("all_parent_transcript_coords"), .f = ~.x %>% paste(collapse = ",")))), .progress = TRUE)
    
    # flatten and distribute the fasta header into all of its child 3FT results
    list_3FT_result_unnest_temp <- purrr::map(.x = list_3FT_result_temp, .f = ~purrr::cross2(.x$`fasta_header`, .$`3FT_info`)) %>% flatten %>%
      ## rename the fasta header entry (first element of each list) then flatten
      purrr::map(.x = ., .f = ~set_names(.x, c("fasta_header", "3FT_result")) %>% flatten)
    
    # rearrange info into a tibble. one translation frame per row.
    # first, tibblise the translation frame elements of the list, then as_tibble the rest.
    cat("rearrange the three-frame translate results into a tibble\n")
    list_3FT_result_unnest_temp_2 <- future_map(.x = list_3FT_result_unnest_temp, .f = function(.x) {
      
      # DEBUG #####
      
      # .x <- list_matched_exon_coords_3FT_result_combined_unnest_temp[[1]]
      
      #############
      
      # cat("now processing entry number", .y, "/", length(list_matched_exon_coords_3FT_result_combined_unnest_temp), "\n")
      
      # define the element indices inside the list which contains info for each frame.
      element.indices_frame_info <- grep(x = names(.x), pattern = "frame_\\d$")
      # element indices containing the virtual peptides
      element.indices_virtual_peptides <- grep(x = names(.x), pattern = "translation_frame_\\d")
      # element indices of the frame info but not the virtual peptides
      element.indices_not_virtual_peptides <- setdiff(element.indices_frame_info, element.indices_virtual_peptides)
      
      # make the tibble of the virtual peptides first.
      tibble_frame_info <- .x[c("translation_frame_0", "translation_frame_1", "translation_frame_2")] %>% as.data.frame %>% t %>% as_tibble(rownames = "translation_frame")
      colnames(tibble_frame_info)[2] <- "parent_transcript_virtual_peptide_sequence"
      # tidy the frame number
      tibble_frame_info[, "translation_frame"] <- gsub(x = tibble_frame_info$translation_frame, pattern = "translation_frame_", replacement = "")
      
      # then make the tibble of the rest. join it onto the existing table
      tibble_frame_info_temp <- .x[element.indices_not_virtual_peptides] %>% as.data.frame %>% t %>% matrix(nrow = 3, byrow = FALSE, dimnames = list(0:2, gsub(x = names(.x[element.indices_not_virtual_peptides]), pattern = "(_frame_\\d)", replacement = "") %>% unique)) %>% as_tibble(rownames = "translation_frame")
      # table join
      tibble_frame_info <- dplyr::full_join(tibble_frame_info, tibble_frame_info_temp, by = "translation_frame")
      
      # finally join the tibblised parts of the list to the non-tibblised.
      list_half_tibblised <- splice("frame_info" = tibble_frame_info, .x[-element.indices_frame_info])
      
      fully_tibblised_list <- list_half_tibblised %>% flatten %>% as_tibble
      
      return(fully_tibblised_list)
      
    }, .progress = TRUE)
    
    # rbind into a summary tibble
    tibble_3FT_result_temp <- list_3FT_result_unnest_temp_2 %>% rbindlist %>% as_tibble
    
    # remove entries with no valid translations. 
    tibble_3FT_result_temp <- tibble_3FT_result_temp[-which(tibble_3FT_result_temp$uORF_valid == "NONE_VALID" & tibble_3FT_result_temp$dORF_valid == "NONE_VALID"), ]
    
    # allocate to temporary tibble by chromosome
    assign(x = paste("list_three_frame_translate_result_temp_chr_", chr, sep = ""), 
           value = tibble_3FT_result_temp)
    
    tictoc::toc()
    
  }
  
  tictoc::toc()
  
  # make a list combining the 3FT results of each chromosome
  # trim off all those elements (per chromosome) that did not have any 3FT result.
  list_three_frame_translate_result <- ls(pattern = "list_three_frame_translate_result_temp_chr_") %>% array_tree %>% purrr::map(.x = ., .f = ~get(.x)) %>% compact
  # rbindlist into a tibble
  tibble_three_frame_translate_result <- list_three_frame_translate_result %>% rbindlist %>% as_tibble
  # remove all rows with duplicated u/dORF_valid columns from the table
  tibble_three_frame_translate_result <- tibble_three_frame_translate_result[-which(duplicated(tibble_three_frame_translate_result[, c("uORF_valid", "dORF_valid")])), ]
  
  tibble_three_frame_translate_result_unique.entries <- dplyr::bind_rows(tibble_three_frame_translate_result[, c("uORF_valid", "fasta_header")] %>% setNames(c("virtual_peptide_sequence", "fasta_header")), tibble_three_frame_translate_result[, c("dORF_valid", "fasta_header")] %>% setNames(c("virtual_peptide_sequence", "fasta_header"))) %>% 
    # remove all "NONE_VALID" sequences
    .[-which(.$virtual_peptide_sequence == "NONE_VALID"), ] %>% dplyr::distinct(., virtual_peptide_sequence, .keep_all = TRUE)
  
  # WE WRITE A TIBBLE CONTAINING THE GENOME COORD-PEPTIDE MAPPING
  # write a table
  write.table(x = tibble_three_frame_translate_result, file = paste(output_dir, output_database_name, "_3FT.summary.info.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # FINALLY! WE WRITE THE FASTA!
  write.fasta(sequences = tibble_three_frame_translate_result_unique.entries$virtual_peptide_sequence %>% array_tree %>% flatten, names = tibble_three_frame_translate_result_unique.entries$fasta_header, file.out = paste(output_dir, output_database_name, ".fasta", sep = ""), open = "w", nbchar = 40, as.string = TRUE)
  
}
  
# finish counting
tictoc::toc()
  
q()

