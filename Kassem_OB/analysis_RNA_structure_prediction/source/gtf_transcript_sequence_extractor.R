script_description <- "# EXTRACTS TRANSCRIPT SEQUENCE FROM GTF ######
Takes a GTF as input and returns a .fasta file of nucleotide sequences of all transcripts. User can choose whether to output in cDNA format or transcript format (i.e. thymine or uracil-based sequences).
This is particularly useful for use with downstream RNA secondary structure prediction software."

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
library(regioneR)

library(tictoc)
# start counting execution time of the whole script
tictoc::tic("Overall execution time")

# manage arguments
# manage arguments 
list_input_arg_info = list(
  "1" = make_option(c("-G", "--organism_name"), type = "character", default = "default_organism", 
                    help = "Optional. It's the string that ends up in the \"OS\" field in the final fasta file", metavar = "character"),
  "2" = make_option(c("-T", "--source_tag"), type = "character", default = "three_frame_translation_junctions", 
                    help = "Compulsory. A character string that will be added to the FASTA headers to indicate the source. It is in the same position as \"sp\" for UniProt fasta files.", metavar = "character"),
  "3" = make_option(c("-R", "--reconstructed_transcript_gtf_path"), type = "character", default = NULL, 
                    help = "Compulsory. path to the actual reference GTF file (e.g. from Cufflinks, Strawberry) that you want to extract nucleotide sequences from. NOT THE CONTAINING DIRECTORY.", metavar = "character"),
  "4" = make_option(c("-F", "--reference_genome_fasta_dir"), type = "character", default = NULL, 
                    help = "Compulsory. path to the directory containing the genome FASTA files. Ideally from Ensembl... you need separate files by chromosomes, NOT the primary assembly. 
              FORMATTING IMPORTANT!!!! MAKE SURE THE REF. GENOME FASTA FILES ARE IN THE FORMAT: <_anything_><chr>.fa e.g. \"Homo_sapiens.GRCh38.dna.chromosome.MT.fa\" OR \"chr16.fa\" OR \"Y.fa\". What will not work: anything which does not have .fa extension e.g. \"chr16.fasta\", anything between the chromosome number and the .fa extension e.g. \"chromosome1.ensembl.fa\"", metavar = "character"),
  "5" = make_option(c("-M", "--rna_mode"), type = "character", default = "rna", 
                    help = "Optional. Specify whether you want RNA or cDNA sequences (i.e. thymine or uracil-based nucleotides?). Default is \"rna\", but you can also specify \"cdna\".", metavar = "character"),
  "6" = make_option(c("-D", "--output_dir"), type = "character", default = NULL, 
                    help = "Compulsory. output file directory. where do you want to save the annotated exon table? IMPORTANT: MUST BE A FULL DIRECTORY AND NOT A FILE PATH. e.g. correct: ~/outputdir/ correct: ~/outputdir incorrect: /outputdir/annotated_exons.txt", metavar = "character"),
  "7" = make_option(c("-O", "--output_name"), type = "character", default = NULL, 
                    help = "Compulsory. output file name, to be saved in the output directory a.k.a. what do you want to save the annotated exon table as? IMPORTANT: MUST BE A STRING WITHOUT THE EXTENSION AND NOT A DIRECTORY. THE .txt EXTENSION WILL AUTOMATICALLY BE ADDED FOR THE OUTPUT FILE. e.g. correct: annotated_sample incorrect: annotated_exons.txt incorrect: annotated_sample/", metavar = "character"),
  "8" = make_option(c("-C", "--ncores"), type = "character", default = 0, 
                    help = "Optional. Number of cores to use. possible inputs: numbers 1 to any integer. By default, uses all cores (ncores = 0). If a single number is specified, it will just tell future to loop thru chromosomes in parallel using the specified core count. If numberxnumber for example 7x4 then 28 cores will be used. 7 for chromosomes and 4 for inside each chromosome.", metavar = "character"), 
  "9" = make_option(c("-H", "--chrmode"), type = "integer", default = 0, 
                     help = "Optional. Specifies which chromosomes to do: select what chromosomes you want considered. possible inputs: numbers 0-2. 0 (default): nuclear chromosomes only i,e, 1:22, X & Y. 1: nuclear + mitochondrial i.e. 1:22, X & Y, M. 2: everything including haplotype/fusion chromosomes etc... this is possible provided the chromosome names.", metavar = "integer"),
  "10" = make_option(c("-N", "--nonchrname"), type = "character", default = NULL, 
                     help = "Compulsory only if you have specified \"--chrmode 2\". nonchromosomal file name. if you want to consider haplotypes, please specify what the reference genome FASTA file for it is called or the script won't know. This single FASTA file must contain all the haplotype information. The program won't try to search for a second file. In ensembl, this file is called something like \"Homo_sapiens.GRCh38.dna.nonchromosomal.fa\" or more generally, \"*nonchromosomal.fa\". So if you want to use that file, then for this option, you would specify \"--nonchrname nonchromosomal\".", metavar = "character"),
  "11" = make_option(c("-V", "--save_workspace_when_done"), type = "character", default = "NO",
                     help = "Turn this on if you want to save the R workspace in the same name as the --output_name. YES: saves at the end. DEBUG: saves at each critical step. NO: doesn't save.", metavar = "character"),
  "12" = make_option(c("-W", "--line_wrap"), type = "integer", default = 0, 
                    help = "Optional. Specifies whether you want to wrap the sequence lines in the final FASTA output file. HIGHLY not recommended if using the FASTA for RNA structure prediction for example, which requires a single line per sequence. Defaults to 0 (no wrapping)", metavar = "integer")
)

input_arg_info <- OptionParser(option_list = list_input_arg_info, description = script_description)
input_args <- input_arg_info %>% parse_args

# check if the input arguments are O.K
if ((list(input_args$reconstructed_transcript_gtf_path, input_args$reference_genome_fasta_dir, input_args$output_name) %>% lapply(is.null) %>% unlist %>% any == TRUE) | 
    (input_args$chrmode == 2 & is.null(input_args$nonchrname) == TRUE)) {
  
  print_help(input_arg_info)
  
  stop("Make sure you entered the arguments correctly", call. = FALSE)
  
}


# DEBUG #######

# organism_name <- "default_organism"
# source_tag <- "debug_alltimepoints_stringtiemerged"
# reconstructed_transcript_gtf_path <- "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/GRAND_OBseries_ref_denovo_reconstructed_stringtiemerged.gtf"
# reconstructed_transcript_gtf_path <- "/mnt/Tertiary/sharedfolder/hg38_ensembl_reference/gtf/Homo_sapiens.GRCh38.98.gtf"
# reference_genome_fasta_dir <- "/mnt/Tertiary/sharedfolder/hg38_ensembl_reference/raw_genome_fasta/dna_by_chr/"
# rna_mode <- "rna"
# output_dir <- "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_RNA_structure_prediction/results_gtf_sequence_extractor/debug/"
# output_name <- "debug_alltimepoints_stringtiemerged"
# ncores <- "30x4"
# chrmode <- 1
# save_workspace_when_done <- "YES"

###############

organism_name <- input_args$organism_name
source_tag <- input_args$source_tag
reconstructed_transcript_gtf_path <- input_args$reconstructed_transcript_gtf_path
reference_genome_fasta_dir <- input_args$reference_genome_fasta_dir
rna_mode <- input_args$rna_mode %>% tolower
output_dir <- input_args$output_dir
output_name <- input_args$output_name
ncores <- input_args$ncores
chrmode <- input_args$chrmode
nonchrname <- input_args$nonchrname
save_workspace_when_done <- input_args$save_workspace_when_done %>% toupper
wrap_lines_in_fasta <- if (input_args$line_wrap < 1) {9e10} else {input_args$line_wrap}

cat("organism_name:", organism_name, "\n")
cat("source_tag:", source_tag, "\n")
cat("reconstructed_transcript_gtf_path:", reconstructed_transcript_gtf_path, "\n")
cat("reference_genome_fasta_dir:", reference_genome_fasta_dir, "\n")
cat("rna_mode:", rna_mode, "\n")
cat("output_dir:", output_dir, "\n")
cat("output_name:", output_name, "\n")
cat("chrmode:", chrmode, "\n")
cat("nonchrname:", nonchrname, "\n")
cat("save_workspace_when_done:", save_workspace_when_done, "\n")
cat("wrap_lines_in_fasta:", wrap_lines_in_fasta, "\n")

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
  
  plan(list(tweak(multiprocess, workers = ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert, gc = TRUE), 
            tweak(multiprocess, workers = ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert, gc = TRUE))
  )
  
  cat((ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert) * (ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert), "cores will be used in total\n")
  cat("first layer:", ncores %>% strsplit(split = "x") %>% unlist %>% .[1] %>% type.convert, "cores\n")
  cat("second layer:", ncores %>% strsplit(split = "x") %>% unlist %>% .[2] %>% type.convert, "cores\n")
  
}

options(future.globals.maxSize = 30000000000, future.fork.enable = TRUE)

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

# FUNCTION TO EXTRACT TRANSCRIPTS WITH JUNCTION-FLANKING EXONS.
# NOTE: to be used with purrr
# details of ONE junction: $chr, $start, $end, $strand
# tibble_gtf_table: rtracklayer::import(...) %>% as_tibble. can be any GTF. reconstructed or reference.
# index: loop progress marker to be used with imap

extract_junction.flanking.exons <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, tolerance_left = 1, tolerance_right = 1, tolerance_inside = 1, tolerance_outside = 0, match_consecutive = TRUE, return_type = "exon") {
  
  # DEBUG ###################
  
  # query_chr = a1$VSR_chr
  # query_start = a1$VSR_start
  # query_end = a1$VSR_end
  # query_strand = a1$strand
  # tibble_gtf_table = tibble_recon_gtf
  # tolerance_left = 1
  # tolerance_right = 1
  # tolerance_inside = 1
  # tolerance_outside = 0
  # match_consecutive = FALSE
  # return_type = "exon"
  
  ###########################
  
  # print(paste("now processing junction number", index))
  
  if (query_strand == "." | query_strand == 0 | query_strand == "*") {
    
    tibble_gtf_subset_flanking_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% .[which(.$start <= ((query_end %>% as.numeric) + 1 + tolerance_outside + tolerance_left) & .$end >= ((query_start %>% as.numeric) - 1 - tolerance_outside - tolerance_left)), ] %>% .[!(.$start <= ((query_end %>% as.numeric) - tolerance_inside - tolerance_right) & .$end >= ((query_start %>% as.numeric) + tolerance_inside + tolerance_left)), ] %>% .[.$type %in% return_type, ]
    
  } else if (query_strand == "+" | query_strand == "-") {
    
    tibble_gtf_subset_flanking_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% .[.$strand == query_strand %>% trimws, ] %>% .[which(.$start <= ((query_end %>% as.numeric) + 1 + tolerance_right) & .$end >= ((query_start %>% as.numeric) - 1 - tolerance_left)), ] %>% .[!(.$start <= ((query_end %>% as.numeric) - tolerance_inside - tolerance_right) & .$end >= ((query_start %>% as.numeric) + tolerance_inside + tolerance_left)), ] %>% .[.$type %in% return_type, ]
    
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
  
  return(list_of_tibbles_flanking_exon_gtf.entries_per_transcript)
  
}

# END extract_junction.flanking.exons()

# FUNCTION TO EXTRACT REFERENCE EXONS WHICH OVERLAP EXACTLY WITH QUERY EXONS
# NOTE: to be used with purrr
# input: spliceregion_list: a list containing details of ONE junction: $chr, $diff_exon_start, $diff_exon_end
# tibble_gtf_table: rtracklayer::import(...) %>% as_tibble. can be any GTF. reconstructed or reference.
# index: loop progress marker to be used with imap

extract_overlapping.exons <- function(query_chr, query_start, query_end, query_strand, tibble_gtf_table, tolerance_left = 0, tolerance_right = 0, tolerance_inside = 0, tolerance_outside = 0, return_type = "exon") {
  
  # DEBUG ###################
  # index <- 1
  # spliceregion_list <- wide_tibble_of_all_unique_VSR_and_exon_coords_array.tree_not_IR[[index]]
  # # tibble_gtf_table <- tibble_ref_gtf
  # tibble_gtf_table <- tibble_recon_gtf
  # stranded = FALSE
  ###########################
  
  # print(paste("now processing junction number", index))
  
  if (query_strand == "." | query_strand == 0 | query_strand == "*") {
    
    # +/- 1 nt tolerance
    tibble_gtf_subset_overlapping_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws, ] %>% 
      .[which(.$start > ((query_start %>% as.numeric) - 1 - tolerance_left - tolerance_outside) & .$end < ((query_end %>% as.numeric) + 1 + tolerance_right + tolerance_outside)), ] %>% 
      .[which((.$start < ((query_start %>% as.numeric) + 1 + tolerance_left + tolerance_inside) & .$end > ((query_end %>% as.numeric) - 1 - tolerance_right - tolerance_inside))), ] %>% 
      .[which(.$type == return_type), ]
    
  } else if (query_strand == "+" | query_strand == "-") {
    
    # +/- 1 nt tolerance
    tibble_gtf_subset_overlapping_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == query_chr %>% trimws &
                                                              tibble_gtf_table$strand == query_strand %>% trimws, ] %>% 
      .[which(.$start > ((query_start %>% as.numeric) - 1 - tolerance_left - tolerance_outside) & .$end < ((query_end %>% as.numeric) + 1 + tolerance_right + tolerance_outside)), ] %>% 
      .[which((.$start < ((query_start %>% as.numeric) + 1 + tolerance_left + tolerance_inside) & .$end > ((query_end %>% as.numeric) - 1 - tolerance_right - tolerance_inside))), ] %>% 
      .[which(.$type == return_type), ]
    
  } else {
    
    stop("Could not match the strand information in the transposed differential-only UNION_junc_coor_table. Make sure that the \"strand\" column in the UNION_junc_coor_table contains only +, - or .")
    
  }
  
  return(tibble_gtf_subset_overlapping_exons)
  
}

# END extract_overlapping.exons()

# FUNCTION to calculate where the exon's amino acid sequence will be within a whole stretch of translated transcript
# which is the ceiling of the distance /3 - translation frame + 1. reverse for reverse strand.
calculate_translation_frame_relative_start_end_position <- function(TL, AUG, ES, EE, frame, greedy = FALSE) {
  # TL: transcript length, AUG: transcript relative first nt. of start codon, ES: exon start (transcript-relative nucleotide position), EE: exon end, frame: 0-2
  # greedy flag means that we will include the AAs which span a junction
  
  if (greedy == FALSE) {
    exon_start_AA_position <- ceiling((ES - frame - AUG) / 3) + 1
    exon_end_AA_position <- floor((EE - frame - AUG) / 3)
  } else if (greedy == TRUE) {
    exon_start_AA_position <- floor((ES - frame - AUG) / 3) + 1
    exon_end_AA_position <- ceiling((EE - frame - AUG) / 3)
  }
  
  return(list("exon_start_AA_position" = exon_start_AA_position, 
              "exon_end_AA_position" = exon_end_AA_position))
  
}

# FUNCTIONS to test if the left/right side of stop codons in an exon are translatable or not. (i.e. whether uORF or dORF exists or not)
find_valid_uORF <- function(AA_sequence, exon_start_AA_position, exon_end_AA_position) {
  
  exon_start_AA_position <- exon_start_AA_position %>% type.convert
  exon_end_AA_position <- exon_end_AA_position %>% type.convert
  
  # add an extra "N" to account for exons starting at the start of the translated sequence
  validity_test <- stringr::str_detect(c("N", AA_sequence[1:(exon_start_AA_position - 1)]) %>% rev %>% paste(collapse = ""), "^[^\\*]+M|^[^\\*]+$")
  
  exonic_uORF_sequence <- AA_sequence[exon_start_AA_position:exon_end_AA_position] %>% paste(collapse = "") %>% strsplit(., split = "\\*") %>% unlist %>% first
  
  # if there is indeed a valid uORF then translate the first part within the exon
  if (validity_test == TRUE & nchar(exonic_uORF_sequence) >= 7) {
    
    return(exonic_uORF_sequence)
    
  } else {
    
    return("NONE_VALID")
    
  }
  
}

find_valid_dORF <- function(AA_sequence, exon_start_AA_position, exon_end_AA_position) {
  
  exon_start_AA_position <- exon_start_AA_position %>% type.convert
  exon_end_AA_position <- exon_end_AA_position %>% type.convert
  
  validity_test <- stringr::str_detect(c(AA_sequence[exon_start_AA_position:exon_end_AA_position]) %>% rev %>% paste(collapse = ""), "^[^\\*]+M|^[^\\*]+$")
  
  exonic_dORF_sequence <- AA_sequence[exon_start_AA_position:exon_end_AA_position] %>% paste(collapse = "") %>% strsplit(., split = "\\*") %>% unlist %>% last
  
  # if there is indeed a valid dORF then translate the first part within the exon
  if (validity_test == TRUE & nchar(exonic_dORF_sequence) >= 7) {
    
    return(exonic_dORF_sequence)
    
  } else {
    
    return("NONE_VALID")
    
  }
  
}

# END find_valid_uORF() and find_valid_dORF()


# BEGIN EXECUTION #################################

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

cat("import reference transcriptome GTF\n")
# extract only protein_coding transcripts
# if user has specified to output the FASTA, then we consider all transcripts. 
# if on poison exon detection mode only, then we go for only protein_coding transcript_biotype.
tibble_recon_gtf <- rtracklayer::import(reconstructed_transcript_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# list-ify by chromosome
list_recon_gtf_subset_by_chr <- tibble_recon_gtf %>% 
  dplyr::group_split(seqnames) %>%
  set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$seqnames %>% unique) %>% unlist)

list_recon_gtf_subset_by_chr <- list_recon_gtf_subset_by_chr[chr_to_run]

cat("checking reference transcriptome GTF\n")
# automatically detect if exons are always numbered in increasing order regardless of strand (common for ref. transcripts)
## sample the first transcript on the negative strand with more than 1 exon
temp_number <- 1

first_transcript_id <- tibble_recon_gtf$transcript_id %>% unique %>% na.omit %>% .[temp_number]

while (tibble_recon_gtf[tibble_recon_gtf$transcript_id == first_transcript_id, "exon_number"] %>% nrow == 1 | 
       tibble_recon_gtf[tibble_recon_gtf$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "+" | 
       tibble_recon_gtf[tibble_recon_gtf$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "." | 
       tibble_recon_gtf[tibble_recon_gtf$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "*") {
  
  temp_number <- temp_number + 1
  
  first_transcript_id <- tibble_recon_gtf$transcript_id %>% unique %>% na.omit %>% .[temp_number]
  
}

# the test condition
max_test <- tibble_recon_gtf[tibble_recon_gtf$transcript_id == first_transcript_id, "exon_number"] %>% unlist %>% na.omit %>% max
max_exon_start_test <- tibble_recon_gtf[tibble_recon_gtf$transcript_id == first_transcript_id & tibble_recon_gtf$exon_number == max_test, "start"] %>% na.omit %>% paste
min_exon_start_test <- tibble_recon_gtf[tibble_recon_gtf$transcript_id == first_transcript_id & tibble_recon_gtf$exon_number == 1, "start"] %>% na.omit %>% paste

exon_order <- NULL

# if the exon 1 comes before the max exon, then the exon order is always increasing
if (min_exon_start_test < max_exon_start_test) {
  
  exon_order <- "increasing"
  
  # if the exon 1 cones after the max exon, then the exon order is stranded.
} else if (min_exon_start_test > max_exon_start_test) {
  
  exon_order <- "stranded"
  
}

cat("checking reference genome directory\n")
vector_ref_genome_paths_by_chr <- paste(reference_genome_fasta_dir, list.files(reference_genome_fasta_dir)[list.files(reference_genome_fasta_dir) %>% grep(., pattern = ".*.fa$")], sep = "")

cat("set up the reference genome assemblies\n")
vector_ref_genome_paths_by_chr_position <- chr_to_run %>% purrr::map(.f = function(a1) {
  
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

cat("retrieve the transcript-relative nucleotide sequence\n")
list_extracted_nucleotide_sequences <- future_map2(
  .x = list_recon_gtf_subset_by_chr,
  .y = vector_ref_genome_fasta_path,
  .f = function(a1, a2) {
    
    if (save_workspace_when_done == "DEBUG") {
      cat(a2, "\n")
    }
    
    # DEBUG ###
    # a1 <- list_recon_gtf_subset_by_chr[[1]]
    # a2 <- vector_ref_genome_fasta_path[[1]]
    ##########
    
    # temporary allocation to ref genome fasta list
    reference_genome_fasta_chr_temp <- seqinr::read.fasta(file = a2, forceDNAtolower = FALSE)
    
    # list-ify the (now) single-chromosome GTF by transcript_id
    list_GTF_by_transcript_id <- a1 %>% 
      dplyr::filter(type == "exon" & strand %in% c("+", "-")) %>%
      dplyr::group_split(transcript_id) %>%
      set_names(x = ., nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist)
    
    list_transcript_relative_nucleotides_per_transcript <- future_map2(
      .x = list_GTF_by_transcript_id,
      .y = names(list_GTF_by_transcript_id),
      .f = function(b1, b2) {
        
        # DEBUG ###
        # b1 <- list_GTF_by_transcript_id %>% .[[19]]
        # b2 <- list_GTF_by_transcript_id %>% names %>% .[[19]]
        ###########
        
        if (save_workspace_when_done == "DEBUG") {
          cat(b2, "\n")
        }
        
        transcript_chr <- b1$seqnames %>% unique
        transcript_strand <- b1$strand %>% unique
        suppressWarnings(gene_name <- b1$gene_name %>% unique)
        if (gene_name %>% is.null == TRUE) {
          gene_name <- b1$gene_id %>% unique
        }
        
        # generate all forward genome-relative coords for each matched transcript
        vector_all_forward_genome_relative_coords_of_transcript <- purrr::map2(.x = b1$start, .y = b1$end, .f = ~.x:.y) %>% unlist %>% unique %>% sort
        # generate stranded coords
        
        if (transcript_strand == "+") {
          vector_all_stranded_genome_relative_coords_of_transcript <- vector_all_forward_genome_relative_coords_of_transcript
          vector_transcript_relative_nucleotides <- reference_genome_fasta_chr_temp[[transcript_chr %>% as.character]][vector_all_stranded_genome_relative_coords_of_transcript]
        } else if (transcript_strand == "-") {
          vector_all_stranded_genome_relative_coords_of_transcript <- vector_all_forward_genome_relative_coords_of_transcript %>% rev
          vector_transcript_relative_nucleotides <- reference_genome_fasta_chr_temp[[transcript_chr %>% as.character]][vector_all_stranded_genome_relative_coords_of_transcript] %>% seqinr::comp(forceToLower = FALSE)
        }
        
        if (rna_mode == "rna") {
          vector_transcript_relative_nucleotides[vector_transcript_relative_nucleotides == "T"] <- "U"
        }
        
        # string-ify the nucleotides
        string_transcript_relative_nucleotides <- vector_transcript_relative_nucleotides %>% paste(collapse = "")
        
        # CREATE FASTA HEADER ###
        # use transcript_id as final identifier
        final_identifier <- b2
        
        # create FASTA header
        fasta_header <- paste(source_tag, 
                              "|", 
                              final_identifier, 
                              "|", 
                              b1$score %>% unique %>% paste(collapse = ","), 
                              " OS=",
                              organism_name, 
                              " GN=",
                              gene_name, sep = "")
        
        # create tibble of the extracted nucleotide sequence along with the fasta header
        tibble_fasta_header_and_nt_sequence <- tibble::tibble("fasta_header" = fasta_header, 
                                                              "transcript_relative_nts" = string_transcript_relative_nucleotides)
        
        return(tibble_fasta_header_and_nt_sequence)
        
      }, .progress = TRUE)
    
  }, .progress = TRUE)

if (save_workspace_when_done == "DEBUG") {
  cat("saving workspace...\n")
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

cat("cleanup\n")

# rbind and tibblise
tibble_all_combined <- list_extracted_nucleotide_sequences %>% flatten %>% rbindlist %>% as_tibble

cat("write .fasta file\n")

# write the fasta
seqinr::write.fasta(sequences = tibble_all_combined$transcript_relative_nts %>% array_tree %>% flatten, 
                    names = tibble_all_combined$fasta_header, 
                    file.out = paste(output_dir, "/", output_name, ".fasta", sep = ""), open = "w", nbchar = wrap_lines_in_fasta, as.string = TRUE)

if (save_workspace_when_done == "YES" | save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# finish counting
tictoc::toc()

q()

