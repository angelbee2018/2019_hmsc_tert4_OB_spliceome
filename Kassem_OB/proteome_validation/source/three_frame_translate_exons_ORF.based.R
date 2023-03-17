script_description <- "# THREE-FRAME TRANSLATION ######
Attempts to translate the nucleotides of all transcripts flanking all the splice junctions given.
Uses PARALLEL PURRR (FURRR) ^___^ is much faster than the last version

BEHAVIOR:
1. Requires input of both Variable Splice Region (VSR) and Differential Regions (DRs) as called by PSI-Sigma. For PSI-Sigma, coordinates are pre-magnetised to the reference genome so there's no need to do nonspecific matching.
2. Automatically detects whether DRs are exon skipping events, partial extensions (A3/5SS) or full extensions (IR).
3. Matching rules are different for exon skipping and exon extension events."

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
library(furrr)
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
  "1" = make_option(c("-E", "--exon_table_path"), type = "character", default = NULL, 
                    help = "Compulsory. path to the actual exon table file (e.g. from PSI-Sigma) that you want to generate a custom database from. NOT THE CONTAINING DIRECTORY. Must be tab-separated, with two options of column names. You can either have three columns with identifier coords like: 1:230124:230245 or chr, start, end and (optional) strand. If the former, the first column must be called \"VSR_coords\", and the second column \"alternative_exon_coords\". Strand info is placed in columns called \"VSR_strand\" or \"alternative_exon_strand\", however the program will only look at one column because it assumes the VSR is on the same strand as the alternative exon. As long as it contains the word \"strand\". If the latter, then VSR_chr, VSR_start, VSR_end etc..., and alternative_exon_chr, alternative_exon_start etc...
                    The final compulsory column is \"splicemode\". These have to at least demarcate IR and non-IR events.
                    If a fasta is to be outputted (three-frame translation mode), then the columns necessary to create the fasta header can be optionally be supplied i.e. gene_name, organism,  custom_identifier.", metavar = "character"),
  "2" = make_option(c("-I", "--intron_retention_string"), type = "character", default = "IR", 
                    help = "Compulsory. A regular expression which matches to all characters in the splicemode column which are associated with IR events.", metavar = "character"),
  "3" = make_option(c("-T", "--source_tag"), type = "character", default = "three_frame_translation_junctions", 
                    help = "Compulsory. A character string that will be added to the FASTA headers to indicate the source. It is in the same position as \"sp\" for UniProt fasta files.", metavar = "character"),
  "4" = make_option(c("-R", "--reconstructed_transcript_gtf_path"), type = "character", default = NULL, 
                    help = "Compulsory. path to the actual reference GTF file (e.g. from Cufflinks, Strawberry) that you want to use as a guide for translation. NOT THE CONTAINING DIRECTORY.", metavar = "character"),
  "5" = make_option(c("-F", "--reference_genome_fasta_dir"), type = "character", default = NULL, 
                    help = "Compulsory. path to the directory containing the genome FASTA files. Ideally from Ensembl... you need separate files by chromosomes, NOT the primary assembly. 
              FORMATTING IMPORTANT!!!! MAKE SURE THE REF. GENOME FASTA FILES ARE IN THE FORMAT: <_anything_><chr>.fa e.g. \"Homo_sapiens.GRCh38.dna.chromosome.MT.fa\" OR \"chr16.fa\" OR \"Y.fa\". What will not work: anything which does not have .fa extension e.g. \"chr16.fasta\", anything between the chromosome number and the .fa extension e.g. \"chromosome1.ensembl.fa\"", metavar = "character"),
  "6" = make_option(c("-D", "--output_dir"), type = "character", default = NULL, 
                    help = "Compulsory. output file directory. where do you want to save the annotated exon table? IMPORTANT: MUST BE A FULL DIRECTORY AND NOT A FILE PATH. e.g. correct: ~/outputdir/ correct: ~/outputdir incorrect: /outputdir/annotated_exons.txt", metavar = "character"),
  "7" = make_option(c("-O", "--output_name"), type = "character", default = NULL, 
                    help = "Compulsory. output file name, to be saved in the output directory a.k.a. what do you want to save the annotated exon table as? IMPORTANT: MUST BE A STRING WITHOUT THE EXTENSION AND NOT A DIRECTORY. THE .txt EXTENSION WILL AUTOMATICALLY BE ADDED FOR THE OUTPUT FILE. e.g. correct: annotated_sample incorrect: annotated_exons.txt incorrect: annotated_sample/", metavar = "character"),
  "8" = make_option(c("-C", "--ncores"), type = "character", default = 0, 
                    help = "Optional. Number of cores to use. possible inputs: numbers 1 to any integer. By default, uses all cores (ncores = 0). If a single number is specified, it will just tell future to loop thru chromosomes in parallel using the specified core count. If numberxnumber for example 7x4 then 28 cores will be used. 7 for chromosomes and 4 for inside each chromosome.", metavar = "character"), 
  "9" = make_option(c("-S", "--use_start_codon"), type = "character", default = "YES", 
                    help = "Optional but you should really choose the right option. This option tells the program whether or not there are start codons provided for each transcript (where available). It's useful if you want to re-annotate e.g. the reference Ensembl GTF for the presence of PTCs, so that the program will not consider all 3 translation frames but only consider the frame containing the annotated start codon. If for any reason a transcript doesn't have an associated start codon, the program will revert to considering all 3 frames. By default, start codons are used where applicable, and 3-frame translation is used whenever a start codon annotation is not available (YES). NO: always 3-frame translation. ALWAYS: discard transcripts without start codon annotation.", metavar = "character"),
  "10" = make_option(c("-H", "--chrmode"), type = "integer", default = 0, 
                     help = "Optional. Specifies which chromosomes to do: select what chromosomes you want considered. possible inputs: numbers 0-2. 0 (default): nuclear chromosomes only i,e, 1:22, X & Y. 1: nuclear + mitochondrial i.e. 1:22, X & Y, M. 2: everything including haplotype/fusion chromosomes etc... this is possible provided the chromosome names.", metavar = "integer"),
  "11" = make_option(c("-N", "--nonchrname"), type = "character", default = NULL, 
                     help = "Compulsory only if you have specified \"--chrmode 2\". nonchromosomal file name. if you want to consider haplotypes, please specify what the reference genome FASTA file for it is called or the script won't know. This single FASTA file must contain all the haplotype information. The program won't try to search for a second file. In ensembl, this file is called something like \"Homo_sapiens.GRCh38.dna.nonchromosomal.fa\" or more generally, \"*nonchromosomal.fa\". So if you want to use that file, then for this option, you would specify \"--nonchrname nonchromosomal\".", metavar = "character"),
  "12" = make_option(c("-V", "--save_workspace_when_done"), type = "character", default = FALSE,
                     help = "Turn this on if you want to save the R workspace in the same name as the --output_name. YES: saves at the end. DEBUG: saves at each critical step. NO: doesn't save.", metavar = "character")
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
# 
# exon_table_path <- "/mnt/LTS/projects/2020_RNA_atlas/results/R_processing_results_PSISigma/atlas_polya_psisigma_LIS_export_for_3FT.txt"
# intron_retention_string <- "IR"
# source_tag <- "atlas_polya_psisigma_exons"
# reference_genome_fasta_dir <- "/mnt/LTS/reference_data/hg38_ensembl_reference/raw_genome_fasta/dna_by_chr/"
# output_dir <- "/mnt/LTS/projects/2020_RNA_atlas/results/results_proteome_validation/"
# output_name <- "atlas_polya_psisigma_exons_3FT"
# # reconstructed_transcript_gtf_path <- "/mnt/LTS/reference_data/hg38_ensembl_reference/gtf/Homo_sapiens.GRCh38.98.gtf"
# reconstructed_transcript_gtf_path <- "/mnt/LTS/projects/2020_RNA_atlas/results/analysis_strawberry_polya/atlas_polya_stringtiemerged.gtf"
# ncores <- "2x16"
# use_start_codon <- "YES"
# chrmode <- 1
# save_workspace_when_done <- "YES"

# edit our psi-sigma exon tables ###

# psisigma_result_path <- "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_PSIsigma/results/run_1.9g_in_parallel_with_denominator_sorted_GTF/R_processing_results/long_tibble_of_psisigma_results_allcomparisons_differential_info_exons1810_dpsi15_DEXSeq_padj0.01_anysig_with_na.txt"

# create an identifier and a chr start end table for debugging
# tibble_psisigma_exon_tables <- read.delim(file = psisigma_result_path, stringsAsFactors = FALSE, sep = "\t", header = TRUE, row.names = NULL) %>% as_tibble
# identifier table
# tibble_alternative_exons_identifier <- tibble_psisigma_exon_tables %>% dplyr::select(event_region_coords, diff_exon_coords, splicemode, matched_gene_names) %>% setNames(c("VSR_coords", "alternative_exon_coords", "splicemode", "gene_name")) %>% add_column("organism" = "Homo sapiens", "source" = "all_PSISigma_results", "custom_identifier" = NA) %>% unique
# write.table(x = tibble_alternative_exons_identifier, file = paste(output_dir, "tibble_alternative_exons_identifier_all.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
# chr start end strand table
# vector_VSR_chr <- gsub(x = tibble_alternative_exons_identifier$VSR_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\1")
# vector_VSR_start <- purrr::map2(.x = gsub(x = tibble_alternative_exons_identifier$VSR_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\2") %>% type.convert,
#                                 .y = gsub(x = tibble_alternative_exons_identifier$VSR_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\3") %>% type.convert,
#                                 .f = ~min(.x, .y)) %>% unlist
# vector_VSR_end <- purrr::map2(.x = gsub(x = tibble_alternative_exons_identifier$VSR_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\2") %>% type.convert,
#                               .y = gsub(x = tibble_alternative_exons_identifier$VSR_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\3") %>% type.convert,
#                               .f = ~max(.x, .y)) %>% unlist
# vector_VSR_strand <- tibble_alternative_exons_identifier$VSR_strand
# 
# vector_alternative_exon_chr <- gsub(x = tibble_alternative_exons_identifier$alternative_exon_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\1")
# vector_alternative_exon_start <- purrr::map2(.x = gsub(x = tibble_alternative_exons_identifier$alternative_exon_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\2") %>% type.convert,
#                                              .y = gsub(x = tibble_alternative_exons_identifier$alternative_exon_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\3") %>% type.convert,
#                                              .f = ~min(.x, .y)) %>% unlist
# vector_alternative_exon_end <- purrr::map2(.x = gsub(x = tibble_alternative_exons_identifier$alternative_exon_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\2") %>% type.convert,
#                                            .y = gsub(x = tibble_alternative_exons_identifier$alternative_exon_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\3") %>% type.convert,
#                                            .f = ~max(.x, .y)) %>% unlist
# vector_alternative_exon_strand <- tibble_alternative_exons_identifier$alternative_exon_strand

# tibble_alternative_exons_chr_start_end_strand <- tibble("VSR_chr" = vector_VSR_chr,
# "VSR_start" = vector_VSR_start,
# "VSR_end" = vector_VSR_end,
# "VSR_strand" = if (vector_VSR_strand %>% is.null == TRUE) {"*"} else {vector_VSR_strand},
# "alternative_exon_chr" = vector_alternative_exon_chr,
# "alternative_exon_start" = vector_alternative_exon_start,
# "alternative_exon_end" = vector_alternative_exon_end,
# "alternative_exon_strand" = if (vector_alternative_exon_strand %>% is.null == TRUE) {"*"} else {vector_VSR_strand},
# "splicemode" = tibble_alternative_exons_identifier$splicemode,
# "gene_name" = tibble_alternative_exons_identifier$gene_name,
# "organism" = tibble_alternative_exons_identifier$organism,
# "custom_identifier" = tibble_alternative_exons_identifier$custom_identifier,)

# write.table(x = tibble_alternative_exons_chr_start_end_strand, file = paste(output_dir, "tibble_alternative_exons_chr.start.end.strand_PSIsigma_differential.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)

# exon_table_path <- paste(output_dir, "tibble_alternative_exons_chr.start.end.strand_PSIsigma_differential.txt", sep = "")

###############

exon_table_path <- input_args$exon_table_path
intron_retention_string <- input_args$intron_retention_string
source_tag <- input_args$source_tag
reconstructed_transcript_gtf_path <- input_args$reconstructed_transcript_gtf_path
reference_genome_fasta_dir <- input_args$reference_genome_fasta_dir
output_dir <- input_args$output_dir
output_name <- input_args$output_name
ncores <- input_args$ncores
use_start_codon <- input_args$use_start_codon
chrmode <- input_args$chrmode
nonchrname <- input_args$nonchrname
save_workspace_when_done <- input_args$save_workspace_when_done

cat("exon_table_path:", exon_table_path, "\n")
cat("intron_retention_string:", intron_retention_string, "\n")
cat("source_tag:", source_tag, "\n")
cat("reconstructed_transcript_gtf_path:", reconstructed_transcript_gtf_path, "\n")
cat("reference_genome_fasta_dir:", reference_genome_fasta_dir, "\n")
cat("output_dir:", output_dir, "\n")
cat("output_name:", output_name, "\n")
cat("use_start_codon:", use_start_codon, "\n")
cat("chrmode:", chrmode, "\n")
cat("nonchrname:", nonchrname, "\n")
cat("save_workspace_when_done:", save_workspace_when_done, "\n")

# if(!dir.exists(output_dir) ) {
#   dir.create(output_dir, recursive = TRUE)}

if (!dir.exists(paste(output_dir, "temp/", sep = "")) ) {
  dir.create(paste(output_dir, "temp/", sep = ""), recursive = TRUE)}

# Open a file to send messages to
# message_divert_path <- file(paste(output_dir, "/", output_name, "_messages.txt", sep = ""), open = "wt")
# Divert messages to that file
# sink(message_divert_path, type = "message")

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

cat("checking alternative exon table\n")
tibble_alternative_exons <- data.table::fread(file = exon_table_path, stringsAsFactors = FALSE, sep = "\t", header = TRUE, check.names = FALSE) %>% as_tibble %>% unique

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# check for whether we have all the complete VSR/alternative exon identifiers or whether we have complete chr start end info.
flag_identifier_exon_tibble <- all(c("VSR_coords", "alternative_exon_coords", "splicemode") %in% colnames(tibble_alternative_exons) == TRUE)
flag_chr_start_end_strand_exon_tibble <- all(c("VSR_start", "VSR_end", "alternative_exon_start", "alternative_exon_end", "splicemode") %in% colnames(tibble_alternative_exons) == TRUE) & any(grepl(x = colnames(tibble_alternative_exons), pattern = "chr") == TRUE)

if (flag_identifier_exon_tibble == TRUE) {
  cat("exon table format detected: genomic identifier coords\n")
  tibble_alternative_exons <- tibble_alternative_exons %>% dplyr::distinct(VSR_coords, alternative_exon_coords, .keep_all = TRUE)
  }
if (flag_chr_start_end_strand_exon_tibble == TRUE) {
  cat("exon table format detected: chr start end strand\n")
  tibble_alternative_exons <- tibble_alternative_exons %>% dplyr::distinct(VSR_start, VSR_end, alternative_exon_start, alternative_exon_end, .keep_all = TRUE)
  }

# if neither are complete, then throw a warning and die.
if (any(c(flag_identifier_exon_tibble, flag_chr_start_end_strand_exon_tibble) == TRUE) != TRUE) {
  stop("Error: make sure the exon table has the correct headers")
}

# load the coords into a master tibble
## if identifier table is in identifier form, then we are going to strsplit.
## if the table is already in chr start end form, then we just load it.

if (flag_identifier_exon_tibble == TRUE) {
  
  vector_VSR_chr <- gsub(x = tibble_alternative_exons$VSR_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\1")
  vector_VSR_start <- future_map2(.x = gsub(x = tibble_alternative_exons$VSR_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\2") %>% type.convert, 
                                  .y = gsub(x = tibble_alternative_exons$VSR_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\3") %>% type.convert, 
                                  .f = ~min(.x, .y), .progress = TRUE, .options = future_options(globals = FALSE)) %>% unlist
  vector_VSR_end <- future_map2(.x = gsub(x = tibble_alternative_exons$VSR_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\2") %>% type.convert, 
                                .y = gsub(x = tibble_alternative_exons$VSR_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\3") %>% type.convert, 
                                .f = ~max(.x, .y), .progress = TRUE, .options = future_options(globals = FALSE)) %>% unlist
  
  vector_alternative_exon_chr <- gsub(x = tibble_alternative_exons$alternative_exon_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\1")
  list_alternative_exon_starts <- tibble_alternative_exons$alternative_exon_coords %>% strsplit(split = "\\;") %>% purrr::map(
    .f = function(a1) {
      
      purrr::map2(.x = gsub(x = a1, pattern = "(.*):(\\d+)-(.*)", replacement = "\\2") %>% type.convert, 
                  .y = gsub(x = a1, pattern = "(.*):(\\d+)-(.*)", replacement = "\\3") %>% type.convert,
                  .f = function(b1, b2) {
                    
                    min(b1, b2) %>% unlist
                    
                  } ) %>% return
      
    } )
  
  list_alternative_exon_ends <- tibble_alternative_exons$alternative_exon_coords %>% strsplit(split = "\\;") %>% purrr::map(
    .f = function(a1) {
      
      purrr::map2(.x = gsub(x = a1, pattern = "(.*):(\\d+)-(.*)", replacement = "\\2") %>% type.convert, 
                  .y = gsub(x = a1, pattern = "(.*):(\\d+)-(.*)", replacement = "\\3") %>% type.convert,
                  .f = function(b1, b2) {
                    
                    max(b1, b2) %>% unlist
                    
                  } ) %>% return
      
    } )
  
  # check if chromosomes are all equal. if they're not, then throw an error and die
  if (purrr::map2(.x = vector_VSR_chr, .y = vector_alternative_exon_chr, .f = ~.x == .y) %>% unlist %>% any == FALSE) {
    stop("VSR chromosome and alternative exon chromosome mismatch detected. Double-check the data.")
  }
  
  tibble_master_alternative_exons_chr_start_end_strand <- tibble("chr" = vector_VSR_chr,
                                                                 "VSR_start" = vector_VSR_start,
                                                                 "VSR_end" = vector_VSR_end,
                                                                 "alternative_exon_starts" = list_alternative_exon_starts,
                                                                 "alternative_exon_ends" = list_alternative_exon_ends,
                                                                 "strand" = if (tibble_alternative_exons %>% dplyr::select(contains("strand")) %>% dim %>% prod == 0) {"*"} else {tibble_alternative_exons %>% dplyr::select(contains("strand")) %>% .[, 1] %>% unlist},
                                                                 "splicemode" = tibble_alternative_exons$splicemode,
                                                                 "gene_name" = if (tibble_alternative_exons$gene_name %>% is.null != TRUE) {tibble_alternative_exons$gene_name} else {NA},
                                                                 "organism" = if (tibble_alternative_exons$organism %>% is.null != TRUE) {tibble_alternative_exons$organism} else {NA},
                                                                 "custom_identifier" = if (tibble_alternative_exons$custom_identifier %>% is.null != TRUE) {tibble_alternative_exons$custom_identifier} else {NA})
  
} else if (flag_chr_start_end_strand_exon_tibble == TRUE) {
  
  # check if chromosomes are all equal if there are multiple chromosome columns. if they're not, then throw an error and die
  if (tibble_alternative_exons %>% dplyr::select(contains("chr")) %>% ncol > 1) {
    if (tibble_alternative_exons %>% dplyr::select(contains("chr")) %>% purrr::array_tree(margin = 1) %>% future_map(~.x %>% unique %>% length > 1, .progress = TRUE) %>% unlist %>% any) {
      stop("VSR chromosome and alternative exon chromosome mismatch detected. Double-check the data.")
    }
  }
  
  tibble_master_alternative_exons_chr_start_end_strand <- tibble("chr" = tibble_alternative_exons %>% dplyr::select(contains("chr")) %>% .[, 1] %>% unlist,
                                                                 "VSR_start" = tibble_alternative_exons$VSR_start,
                                                                 "VSR_end" = tibble_alternative_exons$VSR_end,
                                                                 "alternative_exon_starts" = tibble_alternative_exons$alternative_exon_start,
                                                                 "alternative_exon_ends" = tibble_alternative_exons$alternative_exon_end,
                                                                 "strand" = if (tibble_alternative_exons %>% dplyr::select(contains("strand")) %>% dim %>% prod == 0) {"*"} else {tibble_alternative_exons %>% dplyr::select(contains("strand")) %>% .[, 1] %>% unlist},
                                                                 "splicemode" = tibble_alternative_exons$splicemode,
                                                                 "gene_name" = if (tibble_alternative_exons$gene_name %>% is.null != TRUE) {tibble_alternative_exons$gene_name} else {NA},
                                                                 "organism" = if (tibble_alternative_exons$organism %>% is.null != TRUE) {tibble_alternative_exons$organism} else {NA},
                                                                 "custom_identifier" = if (tibble_alternative_exons$custom_identifier %>% is.null != TRUE) {tibble_alternative_exons$custom_identifier} else {NA})
  
}

cat("import reference transcriptome GTF\n")
# extract only protein_coding transcripts
# if user has specified to output the FASTA, then we consider all transcripts. 
# if on poison exon detection mode only, then we go for only protein_coding transcript_biotype.
tibble_recon_gtf <- rtracklayer::import(reconstructed_transcript_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character) %>% type_convert

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

cat("begin matching VSRs to reference GTF to identify the relevant transcripts.\n")
# for exon skipping events, match transcripts with the exact exon.
# for partial exon extension events, match the interval between the free VSR side and the free differential region side.
# for IR events, do inexact matching.

# Partial exon extension matching: transcripts containing the extension junction as well as harbouring an exon that completely encapsulates the differential region.
# Inexact matching: match the IR junction then index all transcript ids which overlap with it. Filter those with at least one coord in common with the initially matched transcript. Filter these again for transcripts containing an exon which completely encapsulates the IR region. (This actually works really well.)

# UPDATE: IR events can simply match to the "VSR"!!
# UPDATE2: no, they cant. simply because it can't catch nested IR exons.

# IR inexact matching
# match the IR to transcript junctions
# obtain transcripts overlapping the matched transcript AND lying on the same strand AND (if available) the same gene_name

# remove elements where the exonic regions are less than 3nt long 
tibble_master_alternative_exons_chr_start_end_strand <- tibble_master_alternative_exons_chr_start_end_strand[
  purrr::map2(
  .x = tibble_master_alternative_exons_chr_start_end_strand$alternative_exon_starts, 
  .y = tibble_master_alternative_exons_chr_start_end_strand$alternative_exon_ends, 
  .f = function(a1, a2) {
    
    # DEBUG ###
    # a1 <- tibble_master_alternative_exons_chr_start_end_strand$alternative_exon_starts %>% .[[2567]]
    # a2 <- tibble_master_alternative_exons_chr_start_end_strand$alternative_exon_ends %>% .[[2567]]
    ###########
    
    purrr::map2(
      .x = a1 %>% unlist,
      .y = a2 %>% unlist,
      .f = function(b1, b2) {
        
        return((b2 - b1) > 3)
        
      } ) %>% unlist %>% all
    
  } ) %>% unlist,
  ]

# subset the recon GTF and exon table by chromosomes in common so that we can map2 over them
## split exon table by chromosome
list_alternative_exons_by_chr <- tibble_master_alternative_exons_chr_start_end_strand %>% 
  dplyr::group_split(chr) %>%
  set_names(nm = purrr::map(.x = ., .f = ~.x$chr %>% unique) %>% unlist)
## split GTF table by chromosome
list_recon_gtf_subset_by_chr <- tibble_recon_gtf %>% 
  dplyr::group_split(seqnames) %>%
  set_names(nm = purrr::map(.x = ., .f = ~.x$seqnames %>% unique) %>% unlist)
## get and filter for only chromosomes in common
vector_chr_in_commmon <- intersect(names(list_recon_gtf_subset_by_chr), names(list_alternative_exons_by_chr))

cat("checking reference genome directory\n")
vector_ref_genome_paths_by_chr <- paste(reference_genome_fasta_dir, list.files(reference_genome_fasta_dir)[list.files(reference_genome_fasta_dir) %>% grep(., pattern = ".*.fa$")], sep = "")

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

list_alternative_exons_by_chr <- list_alternative_exons_by_chr[vector_chr_in_commmon]
list_recon_gtf_subset_by_chr <- list_recon_gtf_subset_by_chr[vector_chr_in_commmon]

plan(list(tweak(multiprocess, workers = 8),
          tweak(multiprocess, workers = 16)))

# further split each chromosome list into smaller sectors to reduce parallel overhead. we do this by chunking across the genome.
list_alt_exons_sectored0 <- furrr::future_map(
  .x = list_alternative_exons_by_chr,
  .f = function(a1) {
    
    # DEBUG ###
    # a1 <- list_alternative_exons_by_chr[[1]]
    ###########
    
    tibble_alt_exons <- a1 %>% dplyr::arrange(VSR_start)
    
    # get the maximum coordinate range of each row
    tibble_coord_range_of_each_row <- tibble(
      "min" = purrr::map(.x = tibble_alt_exons[, c("VSR_start", "VSR_end", "alternative_exon_starts", "alternative_exon_ends")] %>% purrr::array_tree(), .f = ~min(.x %>% unlist %>% type.convert)) %>% unlist,
      "max" = purrr::map(.x = tibble_alt_exons[, c("VSR_start", "VSR_end", "alternative_exon_starts", "alternative_exon_ends")] %>% purrr::array_tree(), .f = ~max(.x %>% unlist %>% type.convert)) %>% unlist
    )
    
    # group into islands using a 0 flag
    tibble_coords_flag_table <- dplyr::bind_rows(
      tibble("coord" = tibble_coord_range_of_each_row$min, "flag" = 1),
      tibble("coord" = tibble_coord_range_of_each_row$max, "flag" = -1)
    ) %>% dplyr::arrange(`coord`)
    
    vec_cumulative_sum <- tibble_coords_flag_table$flag %>% purrr::accumulate(sum, .dir = "forward")
    
    tibble_island_intervals <- purrr::map2(
      # starts
      .x = tibble_coords_flag_table[which(vec_cumulative_sum == 1 & c(0, vec_cumulative_sum[1:(length(vec_cumulative_sum) - 1)]) == 0), ] %>% .$coord,
      # ends
      .y = tibble_coords_flag_table[which(vec_cumulative_sum == 0 & c(0, vec_cumulative_sum[1:(length(vec_cumulative_sum) - 1)]) == 1), ] %>% .$coord,
      .f = ~tibble("chr" = tibble_alt_exons$chr %>% unique, "start" = .x, "end" = .y, "strand" = tibble_alt_exons$strand %>% unique)
    ) %>% rbindlist(fill = TRUE, use.names = TRUE) %>% as_tibble
    
    # join islands to prevent having too many list elements and having the same transcript appear many many times in the split GTF
    list_islands_joined_sectored <- purrr::map(
      .x = split(x = 1:nrow(tibble_island_intervals), f = ceiling(seq_along(1:nrow(tibble_island_intervals))/200)),
      .f = ~tibble_island_intervals[.x, ])
    ## new, amalgamated intervals
    tibble_island_intervals <- list_islands_joined_sectored %>% purrr::map(~tibble("chr" = .x$chr %>% unique %>% .[1], "start" = .x$start %>% .[1], "end" = .x$end %>% .[nrow(.x)], "strand" = .x$strand %>% unique %>% .[1])) %>% 
      data.table::rbindlist() %>%
      tibble::as_tibble()
    
    L1_list_alt_exons_sectored <- purrr::map(
      .x = tibble_island_intervals %>% purrr::array_tree(),
      .f = function(b1) {
        
        tibble_alt_exons[which(tibble_coord_range_of_each_row$min >= (b1$start %>% type.convert) & tibble_coord_range_of_each_row$max <= (b1$end %>% type.convert)), ] %>% 
          return
        
      } )
    
    return(
      list(
        "L1_list_alt_exons_sectored" = L1_list_alt_exons_sectored,
        "tibble_island_intervals" = tibble_island_intervals
      )
    )
    
  } )

list_alt_exons_sectored <- purrr::map(.x = list_alt_exons_sectored0, .f = ~.x$L1_list_alt_exons_sectored) %>% purrr::flatten()
list_island_intervals <- purrr::map(.x = list_alt_exons_sectored0, .f = ~.x$tibble_island_intervals)

list_recon_gtf_sectored0 <- furrr::future_map2(
  .x = list_island_intervals,
  .y = list_recon_gtf_subset_by_chr,
  .f = function(a1, a2) {
    
    # DEBUG ###
    # a1 <- list_island_intervals[[1]]
    # a2 <- list_recon_gtf_subset_by_chr[[1]]
    ###########
    
    L1_list_island_intervals <- a1 %>% purrr::array_tree()
    
    # we now proceed to get all GTF transcripts which overlap with the intervals
    L1_list_recon_gtf_sectored <- purrr::map(
      .x = L1_list_island_intervals,
      .f = function(b1) {
        
        tibble_overlappiing_transcripts <- a2[a2$type == "transcript", ] %>% .[.$start <= (b1$end %>% type.convert(as.is = TRUE)) & .$end >= (b1$start %>% type.convert(as.is = TRUE)), ]
        
        tibble_overlapping_transcript_features <- dplyr::bind_rows(tibble_overlappiing_transcripts, a2[a2$transcript_id %in% tibble_overlappiing_transcripts$transcript_id, ])
        
        return(tibble_overlapping_transcript_features)
        
      } )
    
    return(L1_list_recon_gtf_sectored)
    
  } )

list_recon_gtf_sectored <- list_recon_gtf_sectored0 %>% purrr::flatten()

print("Number of GTF entries safely omitted:")
print((list_recon_gtf_sectored %>% purrr::map(~.x %>% nrow) %>% unlist %>% sum) - nrow(tibble_recon_gtf))

plan(list(tweak(multiprocess, workers = 8),
          tweak(multiprocess, workers = 16)))

cat("match VSRs to reconstructed transcriptome\n")
# list_recon_entries_matched_to_exons <- 
future_pmap(
  .l = list(
    "a1" = list_alt_exons_sectored,
    "a2" = list_recon_gtf_sectored,
    "a3" = 1:length(list_alt_exons_sectored)
  ),
  .f = function(a1, a2, a3) {
    
    # DEBUG ###
    # a1 <- list_alt_exons_sectored[[1]]
    # a2 <- list_recon_gtf_sectored[[1]]
    ##########
    
    message(a1$chr %>% unique)
    
    list_per_chromosome <- furrr::future_map(
      .x = a1 %>% array_tree,
      .f = function(b1) {
        
        # DEBUG ###
        # b1 <- a1 %>% array_tree %>% .[[4]]
        # b1 <- a1[106, ]
        ###########
        
        # message(b2)
        
        # detect exon extension/skipping
        ## get vector of VSR and exon coords
        vector_VSR_exon_coords <- b1[c("VSR_start", "VSR_end", "alternative_exon_starts", "alternative_exon_ends")] %>% unlist %>% type.convert
        
        # since we use tolerance, we have to think of values as contiguous bands because sometimes the VSR doesn't line up exactly with the alternative exon coords
        tibble_diffs <- tibble("n" = vector_VSR_exon_coords %>% sort %>% .[2:(length(.))], "n_minus_1" = vector_VSR_exon_coords %>% sort %>% .[1:(length(.) - 1)]) %>% dplyr::mutate("diff" = `n` - `n_minus_1`)
        
        ## if no. of unique VSR + alt. exon coords = 4, then it's exon skipping. If 3, then partial exon extension.
        # IR EVENTS
        if (grepl(x = b1$splicemode, pattern = intron_retention_string) == TRUE) {
          
          # VSR MATCHING
          # match transcripts to the IR region
          list_matched_GTF_exon_entries <- extract_overlapping.exons(query_chr = b1$chr, query_start = b1$VSR_start %>% type.convert, query_end = b1$VSR_end %>% type.convert, query_strand = b1$strand, tibble_gtf_table = a2, tolerance_left = 1, tolerance_right = 1, tolerance_inside = 0, tolerance_outside = 0, return_type = "exon") %>% dplyr::group_split(transcript_id) %>% set_names(nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist)
          # get parent transcript entries
          list_parent_GTF_transcript_entries <- purrr::map(.x = list_matched_GTF_exon_entries %>% names, .f = ~a2[which(a2$transcript_id == .x), ]) %>% set_names(list_matched_GTF_exon_entries %>% names)
          
          # INEXACT MATCHING IF VSR MATCHING FAILS
          if (list_matched_GTF_exon_entries %>% length == 0) {
            list_GTF_entries_matched_to_IR_junction <- extract_junction.flanking.exons(query_chr = b1$chr, query_start = b1$alternative_exon_starts %>% unlist %>% type.convert, query_end = b1$alternative_exon_ends %>% unlist %>% type.convert, query_strand = b1$strand, tibble_gtf_table = a2, tolerance_left = 0, tolerance_right = 0, tolerance_inside = 0, tolerance_outside = 0, match_consecutive = TRUE, return_type = "exon")
            # extract overlapping, same strand, (and maybe) same gene_name/gene_id transcripts.
            tibble_IR_junction_parent_transcript_entries <- a2[which(a2$transcript_id %in% names(list_GTF_entries_matched_to_IR_junction) &
                                                               a2$type == "transcript"), ]
            # get overlapping/samestrand/same gene_id transcript entries
            ## determine if it's called "gene_name" or "gene_id"
            tibble_gene_name <- a2 %>% head %>% dplyr::select(contains("gene_id"), contains("gene_name"))
            gene_name_or_gene_id <- colnames(tibble_gene_name[, 1])
            
            if (a2 %>% dplyr::select(contains("gene_id"), contains("gene_name")) %>% ncol != 0) {
              tibble_transcript_entries_overlapping_IR_junction_transcripts <- 
                a2[which(a2$type == "transcript" &
                     a2$start <= (tibble_IR_junction_parent_transcript_entries$end %>% max) &
                     a2$end >= (tibble_IR_junction_parent_transcript_entries$start %>% min) & 
                     a2$strand %in% (tibble_IR_junction_parent_transcript_entries$strand %>% unique) &
                     (a2[, gene_name_or_gene_id] %>% unlist) %in% (tibble_IR_junction_parent_transcript_entries[, gene_name_or_gene_id] %>% unlist)), ]
            } else if (a2 %>% dplyr::select(contains("gene_id"), contains("gene_name")) %>% ncol == 0) {
              tibble_transcript_entries_overlapping_IR_junction_transcripts <- 
                a2[which(a2$type == "transcript" &
                     a2$start <= (tibble_IR_junction_parent_transcript_entries$end %>% max) &
                     a2$end >= (tibble_IR_junction_parent_transcript_entries$start %>% min) & 
                     a2$strand %in% (tibble_IR_junction_parent_transcript_entries$strand %>% unique)), ]
            }
            
            # retrieve the full entries of overlapping transcripts
            tibble_all_candidate_overlapping_full_entries <- a2[which(a2$transcript_id %in% tibble_transcript_entries_overlapping_IR_junction_transcripts$transcript_id), ]
            
            # extract only transcripts with encapsulating exons. These are the final exon matches.
            list_matched_GTF_exon_entries <- 
              tibble_all_candidate_overlapping_full_entries[tibble_all_candidate_overlapping_full_entries$type == "exon" &
                                                         tibble_all_candidate_overlapping_full_entries$start <= b1$alternative_exon_starts %>% unlist %>% type.convert &
                                                         tibble_all_candidate_overlapping_full_entries$end >= b1$alternative_exon_ends %>% unlist %>% type.convert, ] %>% 
              dplyr::group_split(transcript_id) %>% set_names(nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist)
            # get parent transcript entries
            list_parent_GTF_transcript_entries <- purrr::map(.x = list_matched_GTF_exon_entries %>% names, .f = ~tibble_all_candidate_overlapping_full_entries[which(tibble_all_candidate_overlapping_full_entries$transcript_id == .x), ]) %>% set_names(list_matched_GTF_exon_entries %>% names)
          }
          
          ### EXON SKIPPING
        } else if (((which(tibble_diffs$diff > 1) %>% length) + 1) >= 4) {
          
          # exact exon match
          list_matched_GTF_exon_entries <- purrr::map2(
            .x = b1$alternative_exon_starts %>% unlist, 
            .y = b1$alternative_exon_ends %>% unlist, 
            .f = function(c1, c2) {
            
            extract_overlapping.exons(query_chr = b1$chr, query_start = c1 %>% type.convert, query_end = c2 %>% type.convert, query_strand = b1$strand, tibble_gtf_table = a2, tolerance_left = 1, tolerance_right = 1, tolerance_inside = 0, tolerance_outside = 0, return_type = "exon") %>% dplyr::group_split(transcript_id) %>% set_names(nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist) %>% return
            
          } ) %>% purrr::discard(.p = ~length(.x) == 0) 
          
          if (length(list_matched_GTF_exon_entries) > 0) {
            
            list_matched_GTF_exon_entries <- list_matched_GTF_exon_entries %>% purrr::reduce(.f = function(c1, c2) {purrr::map2(.x = c1[intersect(names(c1), names(c2))], .y = c2[intersect(names(c1), names(c2))], .f = ~dplyr::bind_rows(.x, .y)) %>% return} )
            
          }
          
          # need to have consecutive exons 
          list_matched_GTF_exon_entries <- list_matched_GTF_exon_entries %>% purrr::discard(.p = ~nrow(.x) != length(min(.x$exon_number):max(.x$exon_number)))
          
          # get parent transcript entries
          list_parent_GTF_transcript_entries <- purrr::map(.x = list_matched_GTF_exon_entries %>% names, .f = ~a2[which(a2$transcript_id == .x), ]) %>% set_names(list_matched_GTF_exon_entries %>% names)
          
          ### PARTIAL EXON EXTENSIONS
        } else if (((which(tibble_diffs$diff > 1) %>% length) + 1) == 3) {
          
          # ENHANCED EXON EXTENSION MATCHING
          # Strategy: 
          # 1. Determine the common VSR and exon coord.
          # 2. Calculate the genomic coord of the DR + 1 nucleotide.
          # This is done using the equation: genome_coord_first_nucleotide_extending_into_common_region = common + (common_coord-VSR_only_coord)/(abs(common_coord-VSR_only_coord))
          # 3. Find the CONSTITUTIVE (native, non-extended) exon by requiring no overlap with the DR but has either the start or end coord = the DR + 1 nucleotide.
          # 4. The exonic matches are any exon which overlaps the region of the constitutive span AND contains the WHOLE DR on the same strand.
          common_coord <- intersect(c(b1$VSR_start %>% type.convert, b1$VSR_end %>% type.convert), 
                                    c(b1$alternative_exon_starts %>% unlist %>% type.convert, b1$alternative_exon_ends %>% unlist %>% type.convert))
          VSR_only_coord <- setdiff(c(b1$VSR_start %>% type.convert, b1$VSR_end %>% type.convert), common_coord)
          exon_only_coord <- setdiff(c(b1$alternative_exon_starts %>% unlist %>% type.convert, b1$alternative_exon_ends %>% unlist %>% type.convert), common_coord)
          genome_coord_first_nucleotide_extending_into_common_region = common_coord + ((common_coord - VSR_only_coord)/(abs(common_coord - VSR_only_coord)))
          
          tibble_GTF_exons_with_DR_plus_one_nucleotide <- 
            a2[(which(a2$type == "exon" &
                        (a2$start == genome_coord_first_nucleotide_extending_into_common_region |
                           a2$end == genome_coord_first_nucleotide_extending_into_common_region))), ]
          # require no overlap with DR
          tibble_constitutive_exons <- 
            tibble_GTF_exons_with_DR_plus_one_nucleotide[!(tibble_GTF_exons_with_DR_plus_one_nucleotide$start <= max(common_coord, exon_only_coord) &
                                                             tibble_GTF_exons_with_DR_plus_one_nucleotide$end >= min(common_coord, exon_only_coord)), ]
          
          # find exons overlapping with the constitutive and containing the whole DR regions.
          list_matched_GTF_exon_entries <- 
            a2[which(a2$type == "exon" &
                       a2$strand %in% (tibble_constitutive_exons$strand %>% unique) &
                       # constitutive region overlap
                       a2$start <= max(tibble_constitutive_exons$start %>% max, tibble_constitutive_exons$end %>% max) &
                       a2$end >= min(tibble_constitutive_exons$start %>% min, tibble_constitutive_exons$end %>% min) &
                       # DR encapsulation
                       a2$start <= min(common_coord, exon_only_coord) &
                       a2$end >= max(common_coord, exon_only_coord)), ] %>%
            dplyr::group_split(transcript_id) %>% set_names(nm = purrr::map(.x = ., .f = ~.x$transcript_id %>% unique) %>% unlist)
          
          # get parent transcript entries
          list_parent_GTF_transcript_entries <- purrr::map(.x = list_matched_GTF_exon_entries %>% names, .f = ~a2[which(a2$transcript_id == .x), ]) %>% set_names(list_matched_GTF_exon_entries %>% names)
          
        } # END CONDITIONS ###
        
        # GET STRAND ###
        list_matched_strand <- list_parent_GTF_transcript_entries %>% purrr::map(~.x$strand %>% unique)
        
        # START AND STOP CODONS ###
        # extract the start codon entries
        parent_transcript_start_codon <- list_parent_GTF_transcript_entries %>% purrr::map(~.x %>% dplyr::filter(type == "start_codon"))
        # generate a flag per matched transcript indicating the presence of start codon
        if (use_start_codon == "YES") {
          list_start_codon_present <- parent_transcript_start_codon %>% purrr::map(~.x %>% nrow > 0)
        } else if (use_start_codon == "NO") {
          list_start_codon_present <- parent_transcript_start_codon %>% purrr::map(~FALSE)
        }
        
        # extract the stop codon entries
        parent_transcript_stop_codon <- list_parent_GTF_transcript_entries %>% purrr::map(~.x %>% dplyr::filter(type == "stop_codon"))
        # generate a flag per matched transcript indicating the presence of start codon
        list_stop_codon_present <- parent_transcript_stop_codon %>% purrr::map(~.x %>% nrow > 0)
        ###
        
        return(purrr::splice(
          b1,
          "list_matched_GTF_exon_entries" = list_matched_GTF_exon_entries %>% list,
          "list_parent_GTF_transcript_entries" = list_parent_GTF_transcript_entries %>% list,
          "list_matched_strand" = list_matched_strand %>% list,
          "list_start_codon_present" = list_start_codon_present %>% list,
          "parent_transcript_start_codon" = parent_transcript_start_codon %>% list,
          "list_stop_codon_present" = list_stop_codon_present %>% list,
          "parent_transcript_stop_codon" = parent_transcript_stop_codon %>% list
        ))
        
      } )
    
    # keep elements that didn't have any matches to GTF
    list_per_chromosome_unmatched <- list_per_chromosome %>% purrr::keep(.p = ~.x$list_matched_GTF_exon_entries %>% length == 0)
    # prune elements that didn't have any matches to GTF
    list_per_chromosome_pruned <- list_per_chromosome %>% purrr::discard(.p = ~.x$list_matched_GTF_exon_entries %>% length == 0)
    
    # return(list("unmatched_list" = list_per_chromosome_unmatched,
    #             "pruned_list" = list_per_chromosome_pruned))
    
    save(list_per_chromosome_unmatched, file = paste(output_dir, "/temp/", output_name, "_unmatched_list_", a3, "_temp.RData", sep = ""))
    save(list_per_chromosome_pruned, file = paste(output_dir, "/temp/", output_name, "_pruned_list_", a3, "_temp.RData", sep = ""))
    
  }, .progress = TRUE)

if (save_workspace_when_done == "DEBUG") {
  cat("saving workspace...\n")
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# # write the entries that didn't have any matches to the reference GTF.
# tibble_recon_entries_matched_to_exons_unmatched <- list_recon_entries_matched_to_exons %>% purrr::map(~.x$unmatched_list) %>% flatten %>% purrr::map(~.x[c("chr", "VSR_start", "VSR_end", "alternative_exon_starts", "alternative_exon_ends", "strand", "splicemode", "gene_name", "organism", "custom_identifier")] %>% as_tibble) %>% rbindlist %>% as_tibble
# # write
# write.table(x = tibble_recon_entries_matched_to_exons_unmatched %>% tidyr::unnest(c(alternative_exon_starts, alternative_exon_ends)), file = paste(output_dir, "/", output_name, "_unmatched.exons.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# # extract the pruned entries to continue
# list_recon_entries_matched_to_exons_pruned <- list_recon_entries_matched_to_exons %>% purrr::map(~.x$pruned_list)

# load(file = "/mnt/LTS/projects/2020_RNA_atlas/results/results_proteome_validation/list_recon_entries_matched_to_exons_pruned.Rlist")

plan(list(tweak(multicore, workers = 8),
          tweak(multicore, workers = 16)))

# simple round robin callr
## true round robin
## one line per process
## global environment automatically exported via temp file
round_robin_pmap_callr <- function(.l, .f, .num_workers = 1, .temp_path = NULL, .re_export = TRUE, .status_messages_dir, .job_name, ...) {
  
  # DEBUG ###
  # .l <- list("a1" = 1:100)
  # .f <- function(a1) {Sys.sleep(30); return(a1 + 1)}
  # # .f <- function(a1, a2) {return(a1 + a2)}
  # .num_workers <- 8
  # .temp_path <- "/mnt/LTS/projects/2020_RNA_atlas/results/results_proteome_validation/temp/atlas_polya_psisigma_exons_3FT_unmatched_list_131_temp.RData"
  # .status_messages_dir <- "/mnt/LTS/projects/2020_RNA_atlas/results/results_proteome_validation/temp/"
  # .job_name <- output_name
  # .re_export <- FALSE
  
  # .l <- list(
  #   "a1" = 1:length(list_alt_exons_sectored),
  #   "a2" = purrr::map2(.x = list_island_intervals, .y = vector_ref_genome_fasta_path, .f = ~rep(.y, times = nrow(.x))) %>% flatten
  #   # setdiff(1:length(list_alt_exons_sectored), list.files(paste(output_dir, "/temp/", sep = ""), pattern = "list_three_frame_translation") %>% gsub(pattern = ".*_(\\d+)_temp.RData", replacement = "\\1", perl = TRUE) %>% type.convert(as.is = TRUE))
  # )
  # .f <- testfun
  # .num_workers <- 8
  # .temp_path <- paste(output_dir, "/temp/", output_name, "_callr_export_3FT_temp.RData", sep = "")
  # .status_messages_dir <- paste(output_dir, "temp/", sep = "")
  # .job_name <- output_name
  # .re_export <- FALSE
  ###########
  
  if (is.null(.job_name)) {
    .job_name <- as.numeric(Sys.time())
  }
  
  print(.temp_path)
  
  print(file.exists(.temp_path))
  print(.re_export)
  
  print((!file.exists(.temp_path)) | .re_export == TRUE)
  
  # save the current env into temp path
  if (!file.exists(.temp_path) | .re_export == TRUE) {
    save.image(file = .temp_path, ...)
  }
  
  # check if all the list elements are of equal length
  map_length <- unique(unlist(purrr::map(.x = .l, .f = ~length(.x))))
  
  if (length(map_length) != 1) {
    stop("ERROR: list args have incompatible lengths.")
  }
  
  # preallocate worker list
  list_workers <- list()
  
  # initial allocation of tasks to maximum number of workers
  for (..i in 1:.num_workers) {
    
    function_current_instance <- .f
    
    list_current_arguments <- purrr::map(.x = .l, .f = ~.x[[..i]])
    
    formals(fun = function_current_instance) <- list_current_arguments
    
    # essential to spike the exported function with a command to read the environment back in
    body(function_current_instance) <- as.call(purrr::prepend(as.list(body(function_current_instance)), expression(load(file = commandArgs()[2])), before = 2))
    
    assign(
      x = paste("worker_", ..i, sep = ""), 
      value = callr::r_bg(
        cmdargs	= .temp_path,
        func = function_current_instance
      )
    )
    
    list_workers <- purrr::splice(list_workers, get(paste("worker_", ..i, sep = "")))
    names(list_workers)[length(list_workers)] <- paste("worker_", ..i, sep = "")
      
  }
  
  # keep track of iteration progress
  current_map_index <- .num_workers
  
  # keep running while there are no errors
  ## NULL values dont count as nonzero - this is good for us
  while(all(unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) == 0)) {
    
    vector_logical_indices_workers_completed_reported <- unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) == 0
    
    # check on process status and write progress file
    purrr::map2(
      .x = list_workers, 
      .y = names(list_workers), 
      .f = function(a1, a2) {
        
        stdout_lines <- a1$read_output_lines()
        stderr_lines <- a1$read_error_lines()
        
        if (length(stdout_lines) != 0) {write.table(x = stdout_lines, file = paste(.status_messages_dir, .job_name, "_", a2, "_stdout.txt", sep = ""), append = TRUE, row.names = FALSE, quote = FALSE)}
        
        if (length(stderr_lines) != 0) {write.table(x = stderr_lines, file = paste(.status_messages_dir, .job_name, "_", a2, "_stderr.txt", sep = ""), append = TRUE, row.names = FALSE, quote = FALSE)}
        
      } )
    
    # round robin action until the list is fully mapped
    if (current_map_index < map_length) {
      
      number_of_processes_alive <- length(which(unlist(purrr::map(.x = list_workers, .f = ~.x$is_alive()))))
      
      if (number_of_processes_alive < .num_workers) {
        
        new_map_start <- current_map_index + 1
        new_map_end <- min(c(current_map_index + .num_workers - number_of_processes_alive, map_length))
        
        current_map_index <- new_map_end
        
        for (..i in (new_map_start):min(c(new_map_end, map_length))) {
          
          function_current_instance <- .f
          worker_1$read_output_lines()
          list_current_arguments <- purrr::map(.x = .l, .f = ~.x[[..i]])
          
          formals(fun = function_current_instance) <- list_current_arguments
          
          # essential to spike the exported function with a command to read the environment back in
          body(function_current_instance) <- as.call(prepend(as.list(body(function_current_instance)), expression(load(file = commandArgs()[2])), before = 2))
          
          assign(
            x = paste("worker_", ..i, sep = ""), 
            value = callr::r_bg(
              cmdargs	= .temp_path,
              func = function_current_instance
            )
          )
          
          list_workers <- purrr::splice(list_workers, get(paste("worker_", ..i, sep = "")))
          names(list_workers)[length(list_workers)] <- paste("worker_", ..i, sep = "")
          
        }
        
      }
      
      # an actually useful progress bar although rudimentary
      cat(paste("\rPercent map completion: ", current_map_index, "/", map_length, " (", round(x = 100*current_map_index/map_length, digits = 1), "%)", sep = ""))
      
      Sys.sleep(1)
      
    } else if (current_map_index == map_length) {
      
      if (length(which(vector_logical_indices_workers_completed_reported)) == map_length) {
        break()
      }
      
    }
    
  }
  
  if (any(unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) != 0)) {
    purrr::map(.x = list_workers, .f = ~.x$signal(9))
    stop(print("Exit status failure received on workers: ", paste(names(list_workers)[which(unlist(purrr::map(.x = list_workers, .f = ~.x$get_exit_status())) != 0)]), collapse = ","))
  }
  
  rm(list_workers)
  rm(list = ls(pattern = "worker_"))
  
  cat("\n")
  
}

# retrieve all forward nucleotides of each parent transcript
cat("do three frame translation\n")
round_robin_pmap_callr(
  .l = list(
    "a1" = 1:length(list_alt_exons_sectored),
    "a2" = purrr::map2(.x = list_island_intervals, .y = vector_ref_genome_fasta_path, .f = ~rep(.y, times = nrow(.x))) %>% flatten
    # setdiff(1:length(list_alt_exons_sectored), list.files(paste(output_dir, "/temp/", sep = ""), pattern = "list_three_frame_translation") %>% gsub(pattern = ".*_(\\d+)_temp.RData", replacement = "\\1", perl = TRUE) %>% type.convert(as.is = TRUE))
  ),
  .num_workers = 16,
  .temp_path = paste(output_dir, "/temp/", output_name, "_callr_export_3FT_temp.RData", sep = ""),
  .status_messages_dir <- paste(output_dir, "temp/", sep = ""),
  .job_name <- output_name,
  .re_export = FALSE,
  .f = function(a1, a2) {
    
        library(seqinr)
        library(tidyverse)
        library(furrr)
        library(rtracklayer)
        library(data.table)
        library(optparse)
        library(regioneR)
        
        library(tictoc)
        # start counting execution time of the whole script
        tictoc::tic("Overall execution time")
        
        options(future.globals.maxSize = 30000000000, future.fork.enable = TRUE)
        
        plan(list(tweak(multicore, workers = 7)))
        
    # DEBUG ###
    # a1 <- 1:length(list_alt_exons_sectored) %>% .[[12]]
    # a2 <- purrr::map2(.x = list_island_intervals, .y = vector_ref_genome_fasta_path, .f = ~rep(.y, times = nrow(.x))) %>% flatten %>% .[[12]]
    ###########
    
    # message(a1)
    
    load(file = paste(output_dir, "/temp/", output_name, "_pruned_list_", a1, "_temp.RData", sep = ""))
    
    list_L2_result <- furrr::future_map2(
      .x = list_per_chromosome_pruned,
      .y = (1:length(a1)),
      .f = function(b1, b2) {
        
        # DEBUG ###
        # b1 <- list_per_chromosome_pruned[[1]]
        # b2 <- 1:length(a1) %>% .[[1]]
        ###########
        
        # message(b2)
        
        # temporary allocation to ref genome fasta list
        reference_genome_fasta_chr_temp <- seqinr::read.fasta(file = a2, forceDNAtolower = FALSE)
        
        # generate all forward genome-relative coords for each matched transcript
        list_all_forward_genome_relative_coords_of_parent_transcript <- b1$list_parent_GTF_transcript_entries %>% purrr::map(~.x[.x$type == "exon", ] %>% purrr::map2(.x = .$start, .y = .$end, .f = ~.x:.y) %>% unlist %>% unique %>% sort)
        # generate stranded coords
        list_all_stranded_genome_relative_coords_of_parent_transcript <- purrr::map2(
          .x = b1$list_matched_strand,
          .y = list_all_forward_genome_relative_coords_of_parent_transcript,
          .f = function(c1, c2) {
            
            if (c1 == "+") {
              c2 %>% return  
            } else if (c1 == "-") {
              c2 %>% rev %>% return  
            }
            
          } )
        
        # hence generate the forward nucleotides - STRAND-DEPENDENT
        # fwd. nucleotides are reversed for reverse strand.
        # don't forget to complement on top of reversing
        list_all_stranded_nucleotides_of_parent_transcript <- purrr::map2(
          .x = list_all_stranded_genome_relative_coords_of_parent_transcript,
          .y = b1$list_matched_strand,
          .f = function(c1, c2) {
            
            # DEBUG ###
            # c1 <- list_all_stranded_genome_relative_coords_of_parent_transcript[[1]]
            # c2 <- b1$list_matched_strand %>% .[[1]]
            ##########
            
            if (c2 == "+") {
              reference_genome_fasta_chr_temp[[b1$chr %>% paste]][c1] %>% return  
            } else if (c2 == "-") {
              reference_genome_fasta_chr_temp[[b1$chr %>% paste]][c1] %>% seqinr::comp(forceToLower = FALSE) %>% return  
            }
            
          } )
        
        # to account for the presence of a start codon or stop codon, we have to find the transcript-relative position of the first nucleotide of the start codon/last nucleotide of the stop codon.
        # to do this,
        # 1. refer to the "start_codon_present" flag. IF FALSE, then transcript-relative position of the first nt. of start codon is 1. (translate whole transcript). 
        # 1b. Similarly, if there is no stop codon, then the "stop" position will be the end of the transcript.
        # 2. If TRUE, then we refer to the $parent_transcript_start_codon entry.
        # 3. Extract the genome-relative position of the first nt. (STRAND-DEPENDENT)
        # 4. Find the transcript-relative position using the forward coords generated (STRAND-DEPENDENT)
        # 5. Start translating from the transcript-relative first nt. position.
        list_genomic_coord_first_nt_of_start_codon <- purrr::map2(
          .x = b1$list_matched_strand,
          .y = b1$parent_transcript_start_codon,
          .f = function(c1, c2) {
            
            if (c1 == "+") {
              c2$start %>% min %>% return
            } else if (c1 == "-") {
              c2$end %>% max %>% return
            }
            
          } )
        
        list_transcript_relative_first_nt_of_start_codon <- purrr::pmap(
          .l = list(
            "c1" = b1$list_matched_strand,
            "c2" = list_genomic_coord_first_nt_of_start_codon,
            "c3" = list_all_stranded_genome_relative_coords_of_parent_transcript,
            "c4" = b1$list_start_codon_present
            ), 
          .f = function(c1, c2, c3, c4) {
            
            if (c4 == FALSE) {
              return(1)
            } else {
              which(c3 == c2) %>% return
            } 
            
          } )
        
        list_genomic_coord_last_nt_of_stop_codon <- purrr::map2(
          .x = b1$list_matched_strand,
          .y = b1$parent_transcript_stop_codon,
          .f = function(c1, c2) {
            
            if (c1 == "+") {
              c2$end %>% max %>% return
            } else if (c1 == "-") {
              c2$start %>% min %>% return
            }
            
          } )
        
        list_transcript_relative_last_nt_of_stop_codon <- purrr::pmap(
          .l = list(
            "c1" = b1$list_matched_strand,
            "c2" = list_genomic_coord_last_nt_of_stop_codon,
            "c3" = list_all_stranded_genome_relative_coords_of_parent_transcript,
            "c4" = b1$list_stop_codon_present
          ), 
          .f = function(c1, c2, c3, c4) {
            
            if (c4 == FALSE) {
              c3 %>% length %>% return
            } else {
              which(c3 == c2) %>% return
            }
            
          } )
        
        # get transcript-relative start and end positions
        list_transcript_relative_position_alternative_exon_start <- purrr::map(
          .x = list_all_stranded_genome_relative_coords_of_parent_transcript,
          .f = function(c1) {
            
            # DEBUG ###
            # c1 <- list_all_stranded_genome_relative_coords_of_parent_transcript[[1]]
            ###########
            
            which(c1 == max((b1$alternative_exon_starts %>% unlist %>% min %>% type.convert(as.is = TRUE)), min(c1[c1 > (b1$alternative_exon_starts %>% unlist %>% min %>% type.convert(as.is = TRUE))]))) %>% return
            
          } )
        
        list_transcript_relative_position_alternative_exon_end <- purrr::map(
          .x = list_all_stranded_genome_relative_coords_of_parent_transcript,
          .f = function(c1) {
            
            # DEBUG ###
            # c1 <- list_all_stranded_genome_relative_coords_of_parent_transcript[[44]]
            ###########
            
            which(c1 == min((b1$alternative_exon_ends %>% unlist %>% max %>% type.convert(as.is = TRUE)), max(c1[c1 < (b1$alternative_exon_ends %>% unlist %>% max %>% type.convert(as.is = TRUE))]))) %>% return
            
          } )
        
        ## create flag indicating whether the alternative exon is in the 5'/3' UTR. 
        list_flag_alternative_exon_location <- purrr::pmap(
          .l = list(
            "c1" = list_transcript_relative_position_alternative_exon_start,
            "c2" = list_transcript_relative_position_alternative_exon_end,
            "c3" = list_transcript_relative_first_nt_of_start_codon,
            "c4" = list_transcript_relative_last_nt_of_stop_codon
          ), 
          .f = function(c1, c2, c3, c4) {
            
            # 5' UTR check. Alternative exon must contain at least 5 nucleotides which were translated.
            if (c1 < (c3 + 5) & c2 < (c3 + 5)) {
              return("alternative_exon_is_in_five_prime_utr")
              # 3' UTR check.
            } else if (c1 > (c4 - 5) & c2 > (c4 - 5)) {
              return("alternative_exon_is_in_three_prime_utr")
              # 5' UTR overlap check
            } else if (c3 %in% c1:c2) {
              return("alternative_exon_overlaps_five_prime_utr")
              # 3' UTR overlap check
            } else if (c4 %in% c1:c2) {
              return("alternative_exon_overlaps_three_prime_utr")
            } else {
              return("alternative_exon_is_in_translated_region")
            }
            
          } )
        
        # three-frame translation
        # NOTE: if there are not at least 5 exons translatable, then we don't translate.
        list_raw_three_frame_translation <- purrr::pmap(
          .l = list(
            "c1" = list_all_stranded_nucleotides_of_parent_transcript,
            "c2" = list_transcript_relative_first_nt_of_start_codon,
            "c3" = list_transcript_relative_last_nt_of_stop_codon,
            "c4" = b1$list_start_codon_present,
            "c5" = list_flag_alternative_exon_location
          ),
          .f = function(c1, c2, c3, c4, c5) {
            
            # DEBUG ###
            # c1 <- list_all_stranded_nucleotides_of_parent_transcript[[4]]
            # c2 <- list_transcript_relative_first_nt_of_start_codon[[4]]
            # c3 <- list_transcript_relative_last_nt_of_stop_codon[[4]]
            # c4 <- b1$list_start_codon_present %>% .[[4]]
            # c5 <- list_flag_alternative_exon_location[[4]]
            ##########
            
            if (c5 %in% c("alternative_exon_is_in_five_prime_utr", "alternative_exon_is_in_three_prime_utr") == FALSE) {
              list_3FT_result <- nt.sequence_strand_threeframetranslate(vector_forward_nucleotides = c1[c2:c3], strand = "+")
            } else if (c5 %in% c("alternative_exon_is_in_five_prime_utr", "alternative_exon_is_in_three_prime_utr") == TRUE) {
              list_3FT_result <- c(0:2) %>% 
                purrr::map(.f = ~"*") %>% 
                set_names(c("translation_frame_0", "translation_frame_1", "translation_frame_2"))
            } 
            
            # keep only frame 0 if start codon was present
            if (c4 == TRUE) {
              list_3FT_result <- list_3FT_result %>% purrr::map(~list_3FT_result$translation_frame_0)
            }
            
            return(list_3FT_result)
            
          } )
        
        # filter 3FT result for valid ORF
        # Method:
        # 1. Get transcript-relative co-ords of the alternative exon start and end (moved to above)
        # 2. Determine whether the alternative exon was in the translated range (i.e. between the start and stop codons)
        # 3. Find valid ORF after adjusting for the translatable range.
        
        # determine valid ORFs
        ## 1. first determine the effective translated ES and EE.
        ## 2. Determine the valid ORFs:
        ## ONLY THE uORF is used if valid stop codon is present.
        ## If the stop codon is not present, then we are going to have to consider both uORF and dORFs.
        list_valid_ORFs <- purrr::pmap(
          .l = list(
            "c1" = list_flag_alternative_exon_location,
            "c2" = list_transcript_relative_position_alternative_exon_start,
            "c3" = list_transcript_relative_position_alternative_exon_end,
            "c4" = list_all_stranded_genome_relative_coords_of_parent_transcript,
            "c5" = list_raw_three_frame_translation,
            "c6" = b1$list_start_codon_present,
            "c7" = b1$list_stop_codon_present,
            "c8" = list_transcript_relative_first_nt_of_start_codon,
            "c9" = list_transcript_relative_last_nt_of_stop_codon,
            "c10" = 1:length(list_flag_alternative_exon_location)
          ),
          .f = function(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10) {
            
            # DEBUG ###
            # c1 <- list_flag_alternative_exon_location[[1]]
            # c2 <- list_transcript_relative_position_alternative_exon_start[[1]]
            # c3 <- list_transcript_relative_position_alternative_exon_end[[1]]
            # c4 <- list_all_stranded_genome_relative_coords_of_parent_transcript[[1]]
            # c5 <- list_raw_three_frame_translation[[1]]
            # c6 <- b1$list_start_codon_present %>% .[[1]]
            # c7 <- b1$list_stop_codon_present %>% .[[1]]
            # c8 <- list_transcript_relative_first_nt_of_start_codon[[1]]
            # c9 <- list_transcript_relative_last_nt_of_stop_codon[[1]]
            ##########
            
            # message(c10)
            
            if (c1 %in% c("alternative_exon_is_in_five_prime_utr", "alternative_exon_is_in_three_prime_utr") == FALSE) {
              
              transcript_relative_original_ES <- min(c2, c3)
              transcript_relative_original_EE <- max(c2, c3)
              
              transcript_relative_effective_ES <- max(c8, transcript_relative_original_ES)
              transcript_relative_effective_EE <- min(c9, transcript_relative_original_EE)
              
              list_translation_frame_relative_effective_ES <- c(0:2) %>% 
                purrr::map(.f = ~calculate_translation_frame_relative_start_end_position(TL = c4 %>% length, 
                                                                                        AUG = c8, 
                                                                                        ES = transcript_relative_effective_ES, 
                                                                                        EE = transcript_relative_effective_EE,
                                                                                        frame = .x,
                                                                                        greedy = FALSE) %>% .$exon_start_AA_position) %>% 
                set_names(c("translation_frame_0", "translation_frame_1", "translation_frame_2"))
              
              list_translation_frame_relative_effective_EE <- c(0:2) %>% 
                purrr::map(.f = ~calculate_translation_frame_relative_start_end_position(TL = c4 %>% length, 
                                                                                         AUG = c8, 
                                                                                         ES = transcript_relative_effective_ES, 
                                                                                         EE = transcript_relative_effective_EE,
                                                                                         frame = .x,
                                                                                         greedy = FALSE) %>% .$exon_end_AA_position) %>% 
                set_names(c("translation_frame_0", "translation_frame_1", "translation_frame_2"))
              
              # condition for stop codon unavailable - test for both valid uORF and dORF.
              if (c7 == FALSE) {
                
                list_uORF <- purrr::pmap(.l = list(
                  "d1" = c5,
                  "d2" = list_translation_frame_relative_effective_ES,
                  "d3" = list_translation_frame_relative_effective_EE
                ),
                .f = function(d1, d2, d3) {
                  find_valid_uORF(AA_sequence = d1, exon_start_AA_position = d2, exon_end_AA_position = d3) %>% return
                } )
                
                list_dORF <- purrr::pmap(.l = list(
                  "d1" = c5,
                  "d2" = list_translation_frame_relative_effective_ES,
                  "d3" = list_translation_frame_relative_effective_EE
                ),
                .f = function(d1, d2, d3) {
                  find_valid_dORF(AA_sequence = d1, exon_start_AA_position = d2, exon_end_AA_position = d3) %>% return
                } )
                  
              } else if (c.87 == TRUE) {
                
                list_uORF <- purrr::pmap(.l = list(
                  "d1" = c5,
                  "d2" = list_translation_frame_relative_effective_ES,
                  "d3" = list_translation_frame_relative_effective_EE
                ),
                .f = function(d1, d2, d3) {
                  find_valid_uORF(AA_sequence = d1, exon_start_AA_position = d2, exon_end_AA_position = d3) %>% return
                } )
                
                list_dORF <- c5 %>% purrr::map(~"NONE_VALID")
                
              }
              
            } else if (c1 %in% c("alternative_exon_is_in_five_prime_utr", "alternative_exon_is_in_three_prime_utr") == TRUE) {
              # if the alternative exon was not in the translated region at all, then "NONE_VALID" for all translation frames.
              # and effective ES/EE will be null.
              list_translation_frame_relative_effective_ES <- c(0:2) %>% 
                purrr::map(.f = ~NA) %>% 
                set_names(c("translation_frame_0", "translation_frame_1", "translation_frame_2"))
              list_translation_frame_relative_effective_EE <- c(0:2) %>% 
                purrr::map(.f = ~NA) %>% 
                set_names(c("translation_frame_0", "translation_frame_1", "translation_frame_2"))
              list_uORF <- c5 %>% purrr::map(~"NONE_VALID")
              list_dORF <- c5 %>% purrr::map(~"NONE_VALID")
            }
            
            return(list(
              "list_translation_frame_relative_effective_ES" = list_translation_frame_relative_effective_ES,
              "list_translation_frame_relative_effective_EE" = list_translation_frame_relative_effective_EE,
              "list_uORF" = list_uORF,
              "list_dORF" = list_dORF
            ))
            
          } )
        
        list_L2_return <- purrr::splice(b1,
                                        list(
                                          "list_all_stranded_genome_relative_coords_of_parent_transcript" = list_all_stranded_genome_relative_coords_of_parent_transcript %>% list,
                                          "list_genomic_coord_first_nt_of_start_codon" = list_genomic_coord_first_nt_of_start_codon %>% list,
                                          "list_transcript_relative_first_nt_of_start_codon" = list_transcript_relative_first_nt_of_start_codon %>% list,
                                          "list_genomic_coord_last_nt_of_stop_codon" = list_genomic_coord_last_nt_of_stop_codon %>% list,
                                          "list_transcript_relative_last_nt_of_stop_codon" = list_transcript_relative_last_nt_of_stop_codon %>% list,
                                          "list_raw_three_frame_translation" = list_raw_three_frame_translation %>% list,
                                          "list_transcript_relative_position_alternative_exon_start" = list_transcript_relative_position_alternative_exon_start %>% list,
                                          "list_transcript_relative_position_alternative_exon_end" = list_transcript_relative_position_alternative_exon_end %>% list,
                                          "list_flag_alternative_exon_location" = list_flag_alternative_exon_location %>% list,
                                          "list_translation_frame_relative_effective_ES" = list_valid_ORFs %>% purrr::map(~.x$list_translation_frame_relative_effective_ES) %>% list,
                                          "list_translation_frame_relative_effective_EE" = list_valid_ORFs %>% purrr::map(~.x$list_translation_frame_relative_effective_EE) %>% list,
                                          "list_uORF" = list_valid_ORFs %>% purrr::map(~.x$list_uORF) %>% list,
                                          "list_dORF" = list_valid_ORFs %>% purrr::map(~.x$list_dORF) %>% list
                                        )
        )
        
        rm(list = ls() %>% .[. != "list_L2_return"])
        
        # rm(reference_genome_fasta_chr_temp)
        # rm(list_all_forward_genome_relative_coords_of_parent_transcript)
        # rm(list_all_stranded_genome_relative_coords_of_parent_transcript)
        # rm(list_all_stranded_nucleotides_of_parent_transcript)
        # rm(list_genomic_coord_first_nt_of_start_codon)
        # rm(list_transcript_relative_first_nt_of_start_codon)
        # rm(list_genomic_coord_last_nt_of_stop_codon)
        # rm(list_transcript_relative_last_nt_of_stop_codon)
        # rm(list_transcript_relative_position_alternative_exon_start)
        # rm(list_transcript_relative_position_alternative_exon_end)
        # rm(list_flag_alternative_exon_location)
        # rm(list_raw_three_frame_translation)
        # rm(list_valid_ORFs)
        
        gc()
        
        return(list_L2_return)
        
      } ) # L2 - per exon
    
    # return(list_result)
    
    save(list_L2_result, file = paste(output_dir, "/temp/", output_name, "_list_three_frame_translation_full_", a1, "_temp.RData", sep = ""))
    
    return(NULL)
    
  }) # L1 - per chromosome

if (save_workspace_when_done == "YES" | save_workspace_when_done == "DEBUG") {
  cat("saving list...\n")
  save(list_three_frame_translation, file = paste(output_dir, "/", output_name, "_list_three_frame_translation.list", sep = ""))
}

if (save_workspace_when_done == "DEBUG") {
  cat("saving workspace...\n")
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

cat("percolate and tibblise\n")
list_summarised_results <- furrr::future_imap(
  .x = list_three_frame_translation %>% flatten,
  .f = function(a1, a2) {
    
    # DEBUG ###
    # a1 <- list_three_frame_translation %>% flatten %>% .[[40]]
    ###########
    
    cat(a2, "\n")
    
    # good practice: collapse from the inside-out
    # collapse 3FT AA data
    tibble_3FT_AA_data <- purrr::map2(.x = a1[c("list_raw_three_frame_translation", "list_uORF", "list_dORF")],
                                      .y = c("list_raw_three_frame_translation", "list_uORF", "list_dORF"),
                                      .f = function(b1, b2) {
                                        
                                        # DEBUG ###
                                        # b1 <- a1[c("list_raw_three_frame_translation", "list_uORF", "list_dORF")][[1]]
                                        # b2 <- c("list_raw_three_frame_translation", "list_uORF", "list_dORF")[[1]]
                                        ###########
                                        
                                        purrr::map2(.x = b1,
                                                    .y = names(b1),
                                                    .f = function(c1, c2) {
                                                      
                                                      # DEBUG ###
                                                      # c1 <- b1[[1]]
                                                      # c2 <- names(b1) %>% .[[1]]
                                                      ###########
                                                      
                                                      c1 %>% purrr::map(~.x %>% paste(collapse = "")) %>% as_tibble(.name_repair = "unique") %>% t %>% as_tibble(rownames = "translation_frame", .name_repair = "unique") %>% setNames(c("translation_frame", b2)) %>% add_column("matched_transcripts" = c2) %>% return
                                                      
                                                    } ) %>% 
                                          
                                          rbindlist(fill = TRUE) %>% as_tibble %>% return
                                        
                                      } ) %>%
      purrr::reduce(dplyr::full_join, by = c("translation_frame", "matched_transcripts")) %>% 
      dplyr::distinct(matched_transcripts, list_uORF, list_dORF, .keep_all = TRUE)
    
    # tibblise the other 3FT info
    tibble_3FT_other_data <- purrr::map2(.x = a1[c("list_translation_frame_relative_effective_ES", "list_translation_frame_relative_effective_EE")],
                                         .y = c("list_translation_frame_relative_effective_ES", "list_translation_frame_relative_effective_EE"),
                                         .f = function(b1, b2) {
                                           
                                           # DEBUG ###
                                           # b1 <- a1[c("list_translation_frame_relative_effective_ES", "list_translation_frame_relative_effective_EE")][[1]]
                                           # b2 <- c("list_translation_frame_relative_effective_ES", "list_translation_frame_relative_effective_EE")[[1]]
                                           ###########
                                           
                                           purrr::map2(.x = b1,
                                                       .y = names(b1),
                                                       .f = function(c1, c2) {
                                                         
                                                         # DEBUG ###
                                                         # c1 <- b1[[1]]
                                                         # c2 <- names(b1) %>% .[[1]]
                                                         ###########
                                                         
                                                         c1 %>% as_tibble(.name_repair = "unique") %>% t %>% as_tibble(rownames = "translation_frame", .name_repair = "unique") %>% setNames(c("translation_frame", b2)) %>% add_column("matched_transcripts" = c2) %>% return
                                                         
                                                       } ) %>% 
                                             
                                             rbindlist(fill = TRUE) %>% as_tibble %>% return
                                           
                                         } ) %>%
      purrr::reduce(dplyr::full_join, by = c("translation_frame", "matched_transcripts")) %>%
      unique
    
    # collapse the coordinate information with commas and tibblise
    tibble_coordinate_info_per_transcript <- a1$`list_all_stranded_genome_relative_coords_of_parent_transcript` %>% purrr::map(~.x %>% paste(collapse = ",")) %>% as_tibble(.name_repair = "unique") %>% t %>% as_tibble(rownames = "matched_transcripts", .name_repair = "unique") %>% setNames(c("matched_transcripts", "list_all_stranded_genome_relative_coords_of_parent_transcript"))
    
    # tibblise the information per matched transcript
    tibble_other_info_per_matched_transcript <- a1[c("list_matched_strand", "list_start_codon_present", "list_stop_codon_present", "list_genomic_coord_first_nt_of_start_codon", "list_transcript_relative_first_nt_of_start_codon", "list_genomic_coord_last_nt_of_stop_codon", "list_transcript_relative_last_nt_of_stop_codon", "list_transcript_relative_position_alternative_exon_start", "list_transcript_relative_position_alternative_exon_end", "list_flag_alternative_exon_location")] %>% 
      purrr::map2(.x = ., .y = names(.), ~.x %>% as_tibble(.name_repair = "unique") %>% t %>% as_tibble(rownames = "matched_transcripts", .name_repair = "unique") %>% setNames(c("matched_transcripts", .y))) %>%
      purrr::reduce(dplyr::full_join, by = c("matched_transcripts"))
    
    # CREATE FASTA HEADER ###
    # finalise identifiers
    final_identifier <- if (a1$custom_identifier %>% is.na != TRUE) {a1$custom_identifier
    } else {
      paste(source_tag, "_VSR_", 
            a1$chr, ":", a1$VSR_start %>% type.convert, "-", a1$VSR_end %>% type.convert, 
            "_exon_", 
            a1$chr, ":", a1$alternative_exon_starts %>% type.convert, "-", a1$alternative_exon_ends %>% type.convert, 
            if ((a1$strand == "+" | a1$strand == "-") & a1$strand %>% is.na != TRUE) {
              paste(":", a1$strand, sep = "")
            } else {
              ""
            }, sep = "")
    }
    
    # create FASTA header
    fasta_header <- paste(source_tag, 
                          "|", 
                          final_identifier, 
                          "|", 
                          a1$list_parent_GTF_transcript_entries %>% rbindlist %>% .$transcript_id %>% unique %>% paste(collapse = ","), 
                          " OS=",
                          a1$organism, 
                          " GN=",
                          a1$gene_name, sep = "")
    
    # create tibble of exon information including the fasta header
    tibble_exon_info <- purrr::splice(a1[c("chr", "VSR_start", "VSR_end", "alternative_exon_starts", "alternative_exon_ends", "strand", "splicemode", "gene_name")] %>% 
                                        (function(x) {x[c("alternative_exon_starts", "alternative_exon_ends")] <- purrr::map(.x = x[c("alternative_exon_starts", "alternative_exon_ends")], .f = ~.x %>% paste(collapse = ",")); return(x)} ),
                                      "fasta_header" = fasta_header,
                                      "final_identifier" = final_identifier) %>%
      as_tibble
    
    # bind all tables together by left_join
    tibble_3FT_combined <- dplyr::left_join(tibble_3FT_AA_data, tibble_3FT_other_data, by = c("translation_frame", "matched_transcripts"))
    tibble_info_per_transcript_combined <- dplyr::left_join(tibble_coordinate_info_per_transcript, tibble_other_info_per_matched_transcript, by = c("matched_transcripts"))
    
    final_tibble <- dplyr::left_join(tibble_3FT_combined, tibble_info_per_transcript_combined, by = c("matched_transcripts")) %>% 
      dplyr::bind_cols(tibble_exon_info, .)

    return(final_tibble)
    
  }, .progress = TRUE )

if (save_workspace_when_done == "DEBUG") {
  cat("saving workspace...\n")
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

cat("cleanup\n")

# rbind and tibblise
tibble_summarised_results0 <- list_summarised_results %>% rbindlist(use.names = TRUE, fill = TRUE) %>% as_tibble %>% dplyr::mutate_at(.vars = "translation_frame", .funs = function(x) {gsub(x = x, pattern = "translation_frame_", replacement = "")} )

# stack the uORFs and the dORFs, take unique virtual peptide sequences + fasta headers only
tibble_summarised_results1 <- dplyr::bind_rows(tibble_summarised_results0 %>% dplyr::select(-list_uORF) %>% dplyr::rename("virtual_peptide_sequence" = "list_dORF") %>% add_column("ORF_type" = "dORF"),
                                               tibble_summarised_results0 %>% dplyr::select(-list_dORF) %>% dplyr::rename("virtual_peptide_sequence" = "list_uORF") %>% add_column("ORF_type" = "uORF")) %>%
  dplyr::distinct(fasta_header, virtual_peptide_sequence, .keep_all = TRUE)

# filter for virtual peptides less than 7 AA
tibble_summarised_results2 <- tibble_summarised_results1[tibble_summarised_results1$virtual_peptide_sequence %>% purrr::map(~.x %>% nchar >= 7) %>% unlist %>% which, ]

# filter for untranslatable
tibble_summarised_results3 <- tibble_summarised_results2 %>% dplyr::filter(virtual_peptide_sequence != "NONE_VALID")

# get and write the chucked out items
tibble_chucked_out_results <- dplyr::anti_join(tibble_summarised_results1, tibble_summarised_results3)
# write
write.table(x = tibble_chucked_out_results, file = paste(output_dir, "/", output_name, "_discarded.entries.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# filter substrings
## add column to indicate if the virtual peptide per row is a substring of another row
# replace asterisk with \\\\\\\\* for grep
vector_virtual_peptide_sequence <- tibble_summarised_results3$virtual_peptide_sequence

substring_or_not <- future_map2(.x = tibble_summarised_results3$virtual_peptide_sequence, 
                                .y = 1:nrow(tibble_summarised_results3),
                                .f = function(a1, a2) {
                                  
                                  cat(a2, "\n")
                                  
                                  # DEBUG ###
                                  # a1 <- tibble_summarised_results3$virtual_peptide_sequence %>% .[[5]]
                                  # a2 <- 1:nrow(tibble_summarised_results3) %>% .[[5]]
                                  ###########
                                  
                                  (grepl(x = vector_virtual_peptide_sequence[-a2], pattern = a1) & vector_virtual_peptide_sequence[-a2] != a1) %>% any == TRUE
                                  
                                }, .progress = TRUE, .options = future_options(globals = c("vector_virtual_peptide_sequence"))) %>% unlist

## filter
tibble_summarised_results3 <- tibble_summarised_results3 %>% add_column("substring_or_not" = substring_or_not)

# filter out substrings
tibble_summarised_results_no_substring <- tibble_summarised_results3 %>% dplyr::filter(substring_or_not == FALSE)

# tally up the number of valid frames we ended up with
tibble_exons_frame_tally <- tibble_summarised_results_no_substring %>% dplyr::distinct(translation_frame, fasta_header) %>% dplyr::group_by(fasta_header) %>% dplyr::summarise("tally" = n())

cat("\nnumber of exons input: ", tibble_master_alternative_exons_chr_start_end_strand %>% dplyr::distinct(chr, VSR_start, VSR_end, alternative_exon_starts, alternative_exon_ends, strand) %>% nrow, "\n")
cat("\nnumber of exons translated: ", tibble_summarised_results_no_substring$final_identifier %>% unique %>% length, "\n")
cat("\naverage number of translation frames for UNIQUE VSRs + exons: ", mean(tibble_exons_frame_tally$tally), "\n")

# WE WRITE A TIBBLE CONTAINING THE GENOME COORD-PEPTIDE MAPPING
# write a table
write.table(x = tibble_summarised_results_no_substring, file = paste(output_dir, "/", output_name, "_3FT.summary.info.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# FINALLY! WE WRITE THE FASTA!
write.fasta(sequences = tibble_summarised_results_no_substring$virtual_peptide_sequence %>% array_tree %>% flatten, names = tibble_summarised_results_no_substring$fasta_header, file.out = paste(output_dir, "/", output_name, ".fasta", sep = ""), open = "w", nbchar = 40, as.string = TRUE)

# write final exon table as .bed file
exon_bed_table <- tibble_summarised_results_no_substring[, c("chr", "alternative_exon_starts", "alternative_exon_ends", "final_identifier", "strand")] %>% unique %>% setNames(c("chr", "start", "end", "name", "strand")) %>% add_column(., "score" = 1000, .after = "name") %>% type_convert

write.table(exon_bed_table, file = paste(output_dir, "/", output_name, "_exons.bed", sep = ""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = FALSE)

if (save_workspace_when_done == "YES" | save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# finish counting
tictoc::toc()

q()

