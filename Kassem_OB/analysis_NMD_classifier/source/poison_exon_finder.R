script_description <- "# Poison Exon Finder ######
Accepts an input of exon/VSR chr, start, end, strand. Each row describes an alternative exon. (Like from PSI-Sigma).

This program will use the reference GTF (like Ensembl) to deduce whether each alternative exon would contain a premature termination codon (PTC), hence becoming a poison exon.

To do this, the program needs the reference GTF as well as the reference genome assembly.

BEHAVIOUR:
1. Match VSRs (Event Region) to reference transcripts (protein_coding only!!!!)
2. Delete any reference exons overlapping with the alternative exon in question
3. Replace the reference transcript with the exon in question and the associated upstream transcript segment.
4. Three-frame translation, using start-codon frame whenever possible
5. Any PTCs which touch the alternative exon means that the alternative exon is a poison exon. 
6. In addition, we calculate whether the differential exon preserves or alters the original frame.
7. Annotate each Event Region/Target Exon pair with the column \"introduces_PTC\" along with the information associatd with the matched reference transcript."

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
  "1" = make_option(c("-E", "--exon_table_path"), type = "character", default = NULL, 
                    help = "Compulsory. path to the actual exon table file (e.g. from PSI-Sigma) that you want checked for poison exon candidates. NOT THE CONTAINING DIRECTORY. Must be tab-separated, with two options of column names. You can either have three columns with identifier coords like: 1:230124:230245 and splicemode or chr, start, end, splicemode and (optional) strand. If the former, the first column must be called \"VSR_coords\", and the second column \"alternative_exon_coords\". Strand info is placed in columns called \"VSR_strand\" or \"alternative_exon_strand\", however the program will only look at one column because it assumes the VSR is on the same strand as the alternative exon. As long as it contains the word \"strand\". If the latter, then VSR_chr, VSR_start, VSR_end etc..., and alternative_exon_chr, alternative_exon_start etc...
                    For the splicemode, finds transcripts by matching VSRs to flanking exons for non-IR events, and IR regions to flanking exons.
                    You have to manually specify the string which indicates each entry is an IR event, or the program will automatically assume that if the splicemode contains \"IR\" in any capacity, it means the entry is an IR event.
                    If a fasta is to be outputted (three-frame translation mode), then the columns necessary to create the fasta header can be optionally be supplied i.e. gene_name, organism,  custom_identifier.", metavar = "character"),
  "2" = make_option(c("-I", "--intron_retention_string"), type = "character", default = "IR", 
                    help = "Compulsory. A regular expression which matches to all characters in the splicemode column which are associated with IR events.", metavar = "character"),
  "3" = make_option(c("-T", "--source_tag"), type = "character", default = "three_frame_translation_junctions", 
                    help = "Compulsory. A character string that will be added to the FASTA headers to indicate the source. It is in the same position as \"sp\" for UniProt fasta files.", metavar = "character"),
  "4" = make_option(c("-R", "--reference_transcript_gtf"), type = "character", default = NULL, 
                    help = "Compulsory. path to the actual reference GTF file (e.g. from Cufflinks, Strawberry) that you want checked for NMD candidates. NOT THE CONTAINING DIRECTORY.", metavar = "character"),
  "5" = make_option(c("-F", "--reference_genome_fasta_dir"), type = "character", default = NULL, 
                    help = "Compulsory. path to the directory containing the genome FASTA files. Ideally from Ensembl... you need separate files by chromosomes, NOT the primary assembly. 
              FORMATTING IMPORTANT!!!! MAKE SURE THE REF. GENOME FASTA FILES ARE IN THE FORMAT: <_anything_><chr>.fa e.g. \"Homo_sapiens.GRCh38.dna.chromosome.MT.fa\" OR \"chr16.fa\" OR \"Y.fa\". What will not work: anything which does not have .fa extension e.g. \"chr16.fasta\", anything between the chromosome number and the .fa extension e.g. \"chromosome1.ensembl.fa\"", metavar = "character"),
  "6" = make_option(c("-D", "--output_dir"), type = "character", default = NULL, 
                    help = "Compulsory. output file directory. where do you want to save the annotated exon table? IMPORTANT: MUST BE A FULL DIRECTORY AND NOT A FILE PATH. e.g. correct: ~/outputdir/ correct: ~/outputdir incorrect: /outputdir/annotated_exons.txt", metavar = "character"),
  "7" = make_option(c("-O", "--output_name"), type = "character", default = NULL, 
                    help = "Compulsory. output file name, to be saved in the output directory a.k.a. what do you want to save the annotated exon table as? IMPORTANT: MUST BE A STRING WITHOUT THE EXTENSION AND NOT A DIRECTORY. THE .txt EXTENSION WILL AUTOMATICALLY BE ADDED FOR THE OUTPUT FILE. e.g. correct: annotated_sample incorrect: annotated_exons.txt incorrect: annotated_sample/", metavar = "character"),
  "8" = make_option(c("-C", "--ncores"), type = "integer", default = 0, 
                    help = "Optional. Number of cores to use. possible inputs: numbers 1 to any integer. By default, uses all cores (ncores = 0).", metavar = "integer"),
  "9" = make_option(c("-S", "--use_start_codon"), type = "character", default = "YES", 
                    help = "Optional but you should really choose the right option. This option tells the program whether or not there are start codons provided for each transcript (where available). It's useful if you want to re-annotate e.g. the reference Ensembl GTF for the presence of PTCs, so that the program will not consider all 3 translation frames but only consider the frame containing the annotated start codon. If for any reason a transcript doesn't have an associated start codon, the program will revert to considering all 3 frames. By default, start codons are used where applicable, and 3-frame translation is used whenever a start codon annotation is not available (YES). NO: always 3-frame translation. ALWAYS: discard transcripts without start codon annotation.", metavar = "character"),
  "10" = make_option(c("-A", "--output_fasta"), type = "logical", default = FALSE, 
                    help = "Optional. Either TRUE or FALSE. Use this option if you wish to generate a custom peptide database (.fasta format) using the three-frame translation information while we're at it.", metavar = "logical"),
  "11" = make_option(c("-H", "--chrmode"), type = "integer", default = 0, 
                    help = "Optional. Specifies which chromosomes to do: select what chromosomes you want considered. possible inputs: numbers 0-2. 0 (default): nuclear chromosomes only i,e, 1:22, X & Y. 1: nuclear + mitochondrial i.e. 1:22, X & Y, M. 2: everything including haplotype/fusion chromosomes etc... this is possible provided the chromosome names.", metavar = "integer"),
  "12" = make_option(c("-N", "--nonchrname"), type = "character", default = NULL, 
                    help = "Compulsory only if you have specified \"--chrmode 2\". nonchromosomal file name. if you want to consider haplotypes, please specify what the reference genome FASTA file for it is called or the script won't know. This single FASTA file must contain all the haplotype information. The program won't try to search for a second file. In ensembl, this file is called something like \"Homo_sapiens.GRCh38.dna.nonchromosomal.fa\" or more generally, \"*nonchromosomal.fa\". So if you want to use that file, then for this option, you would specify \"--nonchrname nonchromosomal\".", metavar = "character"),
  "13" = make_option(c("-V", "--save_workspace_when_done"), type = "character", default = FALSE,
                     help = "Turn this on if you want to save the R workspace in the same name as the --output_name. YES: saves at the end. DEBUG: saves at each critical step. NO: doesn't save.", metavar = "character")
)

input_arg_info <- OptionParser(option_list = list_input_arg_info, description = script_description)
input_args <- input_arg_info %>% parse_args

# check if the input arguments are O.K
if ((list(input_args$reference_transcript_gtf, input_args$reference_genome_fasta_dir, input_args$output_name) %>% lapply(is.null) %>% unlist %>% any == TRUE) | 
    (input_args$chrmode == 2 & is.null(input_args$nonchrname) == TRUE)) {
  
  print_help(input_arg_info)
  
  stop("Make sure you entered the arguments correctly", call. = FALSE)
  
}


# DEBUG #######

# reference_transcript_gtf_path <- "Z:/hg38_ensembl_reference/gtf/Homo_sapiens.GRCh38.98.gtf"
# reference_transcript_gtf_path <- "/media/Ubuntu/sharedfolder/hg38_ensembl_reference/gtf/Homo_sapiens.GRCh38.98.gtf"
# reference_genome_fasta_dir <- "Z:/hg38_ensembl_reference/raw_genome_fasta/genome_fasta_extract2/"
# reference_genome_fasta_dir <- "/media/Ubuntu/sharedfolder/hg38_ensembl_reference/raw_genome_fasta/genome_fasta_extract2/"
# output_dir <- "Z:/PGNEXUS_kassem_MSC/Kassem_OB/analysis_NMD_classifier/results/"
# output_dir <- "/media/Ubuntu/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_NMD_classifier/results/"
# output_name <- "long_tibble_of_psisigma_results_allcomparisons_differential_poison_exon_finder.txt"
# intron_retention_string <- "IR"
# reference_transcript_gtf_path <- "/media/Ubuntu/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/GRAND_OBseries_ref_denovo_reconstructed_stringtiemerged.gtf"
# reference_genome_fasta_dir <- "/media/Ubuntu/sharedfolder/hg38_ensembl_reference/raw_genome_fasta/genome_fasta_extract2/"
# output_name <- "Homo_sapiens.GRCh38.98_NMDflagger_qualitycheck"
# ncores <- 6
# use_start_codon <- "YES"
# output_fasta <- TRUE
# chrmode <- 1
# psisigma_result_path <- "Z:/PGNEXUS_kassem_MSC/Kassem_OB/analysis_NMD_classifier/results/long_tibble_of_psisigma_DEXSeq_results.txt"
# psisigma_result_path <- "/media/Ubuntu/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_NMD_classifier/results/long_tibble_of_psisigma_DEXSeq_results.txt"
# save_workspace_when_done <- "YES"

# edit our psi-sigma exon tables ###
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
                                                        # "source" = tibble_alternative_exons_identifier$source,
                                                        # "custom_identifier" = tibble_alternative_exons_identifier$custom_identifier,)

# write.table(x = tibble_alternative_exons_chr_start_end_strand, file = paste(output_dir, "tibble_alternative_exons_chr.start.end.strand_all.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
# 
# exon_table_path <- paste("/media/Ubuntu/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_PSIsigma/results/run_1.9g_in_parallel_with_denominator_sorted_GTF/R_processing_results/long_tibble_of_psisigma_allcomparisons_differential_exons1810_dpsi15_DEXSeq_padj0.01_anysig_with_na_poison_exon_finder_3FT.txt", sep = "")
# exon_table_path <- paste(output_dir, "tibble_alternative_exons_chr.start.end.strand_all.txt", sep = "")

###############

exon_table_path <- input_args$exon_table_path
intron_retention_string <- input_args$intron_retention_string
source_tag <- input_args$source_tag
reference_transcript_gtf_path <- input_args$reference_transcript_gtf
reference_genome_fasta_dir <- input_args$reference_genome_fasta_dir
output_dir <- input_args$output_dir
output_name <- input_args$output_name
ncores <- input_args$ncores
use_start_codon <- input_args$use_start_codon
output_fasta <- input_args$output_fasta
chrmode <- input_args$chrmode
nonchrname <- input_args$nonchrname
save_workspace_when_done <- input_args$save_workspace_when_done

cat("exon_table_path:", exon_table_path, "\n")
cat("intron_retention_string:", intron_retention_string, "\n")
cat("source_tag:", source_tag, "\n")
cat("reference_transcript_gtf_path:", reference_transcript_gtf_path, "\n")
cat("reference_genome_fasta_dir:", reference_genome_fasta_dir, "\n")
cat("output_dir:", output_dir, "\n")
cat("output_name:", output_name, "\n")
cat("use_start_codon:", use_start_codon, "\n")
cat("output_fasta:", output_fasta, "\n")
cat("chrmode:", chrmode, "\n")
cat("nonchrname:", nonchrname, "\n")
cat("save_workspace_when_done:", save_workspace_when_done, "\n")

# manage parrallellisation rrlllRll

if (ncores != 0) {
  number_of_workers <- ncores
} 

cat(number_of_workers, "cores will be used\n")
options(future.globals.maxSize = 30000000000, mc.cores = number_of_workers, future.fork.enable = TRUE)
future::plan(multiprocess)

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
  # tibble_gtf_table = tibble_reference_gtf
  # tolerance_left = 1
  # tolerance_right = 1
  # tolerance_inside = 1
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
  
  return(list_of_tibbles_flanking_exon_gtf.entries_per_transcript)
  
}

# END extract_junction.flanking.exons_JUM()

# FUNCTIONS to test if the left/right side of stop codons in an exon are translatable or not. (i.e. whether uORF or dORF exists or not)
find_valid_uORF <- function(AA_sequence, exon_start_AA_position, exon_end_AA_position) {
  
  exon_start_AA_position <- exon_start_AA_position %>% as.numeric
  exon_end_AA_position <- exon_end_AA_position %>% as.numeric
  
  validity_test <- stringr::str_detect(AA_sequence[1:(exon_start_AA_position - 1)] %>% rev %>% paste(collapse = ""), "^[^\\*]+M|^[^\\*]+$")
  
  exonic_uORF_sequence <- AA_sequence[exon_start_AA_position:exon_end_AA_position] %>% paste(collapse = "") %>% strsplit(., split = "\\*") %>% unlist %>% first
  
  # if there is indeed a valid uORF then translate the first part within the exon
  if (validity_test == TRUE & nchar(exonic_uORF_sequence) >= 7) {
    
    return(exonic_uORF_sequence)
    
  } else {
    
    return("NONE_VALID")
    
  }
  
}

find_valid_dORF <- function(AA_sequence, exon_start_AA_position, exon_end_AA_position) {
  
  exon_start_AA_position <- exon_start_AA_position %>% as.numeric
  exon_end_AA_position <- exon_end_AA_position %>% as.numeric
  
  validity_test <- stringr::str_detect(AA_sequence[exon_start_AA_position:exon_end_AA_position] %>% rev %>% paste(collapse = ""), "^[^\\*]+M|^[^\\*]+$")
  
  exonic_dORF_sequence <- AA_sequence[exon_start_AA_position:exon_end_AA_position] %>% paste(collapse = "") %>% strsplit(., split = "\\*") %>% unlist %>% last
  
  # if there is indeed a valid dORF then translate the first part within the exon
  if (validity_test == TRUE & nchar(exonic_dORF_sequence) >= 7) {
    
    return(exonic_dORF_sequence)
    
  } else {
    
    return("NONE_VALID")
    
  }
  
}

# END find_valid_uORF() and find_valid_dORF()

# BEGIN EXECUTION ###########################

cat("checking alternative exon table\n")
tibble_alternative_exons <- read.delim(file = exon_table_path, stringsAsFactors = FALSE, sep = "\t", header = TRUE, row.names = NULL, check.names = FALSE) %>% as_tibble %>% unique

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# check for whether we have all the complete VSR/alternative exon identifiers or whether we have complete chr start end info.
flag_identifier_exon_tibble <- all(c("VSR_coords", "alternative_exon_coords", "splicemode") %in% colnames(tibble_alternative_exons) == TRUE)
flag_chr_start_end_strand_exon_tibble <- all(c("VSR_start", "VSR_end", "alternative_exon_start", "alternative_exon_end", "splicemode") %in% colnames(tibble_alternative_exons) == TRUE) & any(grepl(x = colnames(tibble_alternative_exons), pattern = "chr") == TRUE)

if (flag_identifier_exon_tibble == TRUE) {cat("exon table format detected: genomic identifier coords\n")}
if (flag_chr_start_end_strand_exon_tibble == TRUE) {cat("exon table format detected: chr start end strand\n")}

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
  vector_alternative_exon_start <- purrr::map2(.x = gsub(x = tibble_alternative_exons$alternative_exon_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\2") %>% type.convert, 
                                               .y = gsub(x = tibble_alternative_exons$alternative_exon_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\3") %>% type.convert, 
                                               .f = ~min(.x, .y)) %>% unlist
  vector_alternative_exon_end <- purrr::map2(.x = gsub(x = tibble_alternative_exons$alternative_exon_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\2") %>% type.convert, 
                                             .y = gsub(x = tibble_alternative_exons$alternative_exon_coords, pattern = "(.*):(\\d+)-(.*)", replacement = "\\3") %>% type.convert, 
                                             .f = ~max(.x, .y)) %>% unlist

  # check if chromosomes are all equal. if they're not, then throw an error and die
  if (purrr::map2(.x = vector_VSR_chr, .y = vector_alternative_exon_chr, .f = ~.x == .y) %>% unlist %>% any == FALSE) {
    stop("VSR chromosome and alternative exon chromosome mismatch detected. Double-check the data.")
  }
  
  tibble_master_alternative_exons_chr_start_end_strand <- tibble("chr" = vector_VSR_chr,
                                                          "VSR_start" = vector_VSR_start,
                                                          "VSR_end" = vector_VSR_end,
                                                          "alternative_exon_start" = vector_alternative_exon_start,
                                                          "alternative_exon_end" = vector_alternative_exon_end,
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
         "alternative_exon_start" = tibble_alternative_exons$alternative_exon_start,
         "alternative_exon_end" = tibble_alternative_exons$alternative_exon_end,
         "strand" = if (tibble_alternative_exons %>% dplyr::select(contains("strand")) %>% dim %>% prod == 0) {"*"} else {tibble_alternative_exons %>% dplyr::select(contains("strand")) %>% .[, 1] %>% unlist},
         "splicemode" = tibble_alternative_exons$splicemode,
         "gene_name" = if (tibble_alternative_exons$gene_name %>% is.null != TRUE) {tibble_alternative_exons$gene_name} else {NA},
         "organism" = if (tibble_alternative_exons$organism %>% is.null != TRUE) {tibble_alternative_exons$organism} else {NA},
         "custom_identifier" = if (tibble_alternative_exons$custom_identifier %>% is.null != TRUE) {tibble_alternative_exons$custom_identifier} else {NA})
  
}

cat("checking reference genome directory\n")
vector_ref_genome_paths_by_chr <- paste(reference_genome_fasta_dir, list.files(reference_genome_fasta_dir)[list.files(reference_genome_fasta_dir) %>% grep(., pattern = ".*.fa$")], sep = "")
cat("import reference transcriptome GTF\n")
# extract only protein_coding transcripts
# if user has specified to output the FASTA, then we consider all transcripts. 
# if on poison exon detection mode only, then we go for only protein_coding transcript_biotype.
if (output_fasta == TRUE) {
  tibble_reference_gtf <- rtracklayer::import(reference_transcript_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
} else if (output_fasta == FALSE) {
  tibble_reference_gtf <- rtracklayer::import(reference_transcript_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(transcript_biotype == "protein_coding") 
}

cat("checking reference transcriptome GTF\n")
# automatically detect if exons are always numbered in increasing order regardless of strand (common for ref. transcripts)
## sample the first transcript on the negative strand with more than 1 exon
temp_number <- 1

first_transcript_id <- tibble_reference_gtf$transcript_id %>% unique %>% na.omit %>% .[temp_number]

while (tibble_reference_gtf[tibble_reference_gtf$transcript_id == first_transcript_id, "exon_number"] %>% nrow == 1 | 
       tibble_reference_gtf[tibble_reference_gtf$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "+" | 
       tibble_reference_gtf[tibble_reference_gtf$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "." | 
       tibble_reference_gtf[tibble_reference_gtf$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "*") {
  
  temp_number <- temp_number + 1
  
  first_transcript_id <- tibble_reference_gtf$transcript_id %>% unique %>% na.omit %>% .[temp_number]
  
}

# the test condition
max_test <- tibble_reference_gtf[tibble_reference_gtf$transcript_id == first_transcript_id, "exon_number"] %>% unlist %>% na.omit %>% max
max_exon_start_test <- tibble_reference_gtf[tibble_reference_gtf$transcript_id == first_transcript_id & tibble_reference_gtf$exon_number == max_test, "start"] %>% na.omit %>% paste
min_exon_start_test <- tibble_reference_gtf[tibble_reference_gtf$transcript_id == first_transcript_id & tibble_reference_gtf$exon_number == 1, "start"] %>% na.omit %>% paste

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

# list-ify the master exon table for looping
list_tibble_master_alternative_exons_chr_start_end_strand_array.tree <- tibble_master_alternative_exons_chr_start_end_strand %>% unique %>% purrr::array_tree()
# list_tibble_master_alternative_exons_chr_start_end_strand_array.tree <- tibble_master_alternative_exons_chr_start_end_strand %>% unique %>% .[1:1000,] %>% purrr::array_tree()

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# get the transcript_ids matched by VSR or IR regions
list_GTF_entries_matched_to_VSRs <- future_map(.x = list_tibble_master_alternative_exons_chr_start_end_strand_array.tree, 
                                               .f = function(a1) {
                                                 
                                                 # DEBUG ###
                                                 # a1 <- list_tibble_master_alternative_exons_chr_start_end_strand_array.tree[[4]]
                                                 ###########
                                                 
                                                 # if it's an IR event, we have to match the IR exon, not the VSR.
                                                 if (grepl(x = a1$splicemode, pattern = intron_retention_string) == TRUE) {
                                                   
                                                   list_tibbles_GTF_entries_matched_to_VSR <- extract_junction.flanking.exons(query_chr = a1$chr, query_start = a1$alternative_exon_start %>% type.convert, query_end = a1$alternative_exon_end %>% type.convert, query_strand = a1$strand, tibble_gtf_table = tibble_reference_gtf, tolerance_left = 0, tolerance_right = 0, tolerance_inside = 0, tolerance_outside = 0, match_consecutive = FALSE, return_type = "exon")
                                                   
                                                 } else {
                                                   
                                                   list_tibbles_GTF_entries_matched_to_VSR <- extract_junction.flanking.exons(query_chr = a1$chr, query_start = a1$VSR_start %>% type.convert, query_end = a1$VSR_end %>% type.convert, query_strand = a1$strand, tibble_gtf_table = tibble_reference_gtf, tolerance_left = 0, tolerance_right = 0, tolerance_inside = 0, tolerance_outside = 0, match_consecutive = FALSE, return_type = "exon")
                                                   
                                                 }
                                                 
                                                 # retrieve the entire entry for each transcript as well while we're at it.
                                                 list_tibbles_GTF_entries_parent_transcript <- purrr::map(.x = list_tibbles_GTF_entries_matched_to_VSR %>% names,
                                                                                                          .f = function(b1) {
                                                                                                            
                                                                                                            tibble_reference_gtf %>% dplyr::filter(transcript_id == b1) %>% return
                                                                                                            
                                                                                                          } ) %>% set_names(list_tibbles_GTF_entries_matched_to_VSR %>% names)
                                                 
                                                 # get rid of transcripts without start codon annotations if use_start_codon = ALWAYS
                                                 if (use_start_codon == "ALWAYS") {
                                                   list_tibbles_GTF_entries_parent_transcript <- list_tibbles_GTF_entries_parent_transcript %>% purrr::discard(.p = ~.x %>% dplyr::filter(type == "start_codon") %>% dim %>% prod == 0)
                                                 } else if (use_start_codon == "NO") {
                                                   list_tibbles_GTF_entries_parent_transcript <- list_tibbles_GTF_entries_parent_transcript %>% purrr::map(.f = ~.x %>% dplyr::filter(type != "start_codon"))
                                                 }
                                                 
                                                 # add a flag stating whether start codon exists for each matched transcript
                                                 # if TRUE, then start_codon present.
                                                 list_logical_start_codon_present <- list_tibbles_GTF_entries_parent_transcript %>% purrr::map(.f = ~.x %>% dplyr::filter(type == "start_codon") %>% dim %>% prod > 0)
                                                 
                                                 return(a1 %>% purrr::splice("list_tibbles_GTF_entries_matched_to_VSR" = list(list_tibbles_GTF_entries_matched_to_VSR),
                                                                             "list_tibbles_GTF_entries_parent_transcript" = list(list_tibbles_GTF_entries_parent_transcript),
                                                                             "list_logical_start_codon_present" = list(list_logical_start_codon_present)))
                                                 
                                               }, .progress = TRUE, .options = future_options(globals = c("tibble_reference_gtf", "extract_junction.flanking.exons", "intron_retention_string", "use_start_codon", "dplyr"))) %>% set_names(1:length(list_tibble_master_alternative_exons_chr_start_end_strand_array.tree))

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# prune L1 elements which don't have VSR matches.
list_GTF_entries_matched_to_VSRs_pruned <- list_GTF_entries_matched_to_VSRs %>% purrr::discard(.p = ~length(.x$list_tibbles_GTF_entries_matched_to_VSR) == 0)

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

cat("modify each matched transcript by tacking on the alternative exon.\n")
# for reference entries overlapping the alternative exon, remove them.
list_modified_reference_transcript_segments <- future_imap(.x = list_GTF_entries_matched_to_VSRs_pruned,
                                               .f = function(a1, a2) { 
                                                 
                                                 # DEBUG ###
                                                 # a1 <- list_GTF_entries_matched_to_VSRs[[1]]
                                                 ###########
                                                 
                                                 if (save_workspace_when_done == "DEBUG") {
                                                   cat(a2, "\n")
                                                 }
                                                 
                                                 # get the exon numbers of the reference exons overlapping the alternative exon
                                                 list_ref_exon_numbers_overlapping_alternative_exon <- purrr::map(
                                                   .x = a1$list_tibbles_GTF_entries_parent_transcript, 
                                                   .f = function(b1) {
                                                     
                                                     # DEBUG ###
                                                     # b1 <- a1$list_tibbles_GTF_entries_parent_transcript %>% .[[1]]
                                                     ###########
                                                     
                                                     # extract exonic entries of the parent transcript only
                                                     tibble_exonic_parent_transcript_entries <- b1 %>% dplyr::filter(type == "exon")
                                                     
                                                     # test for reference entry overlap with alternative exon
                                                     row.indices_reference_exons_overlapping_with_alternative_exon <- 
                                                       which(tibble_exonic_parent_transcript_entries$start <= a1$alternative_exon_end %>% type.convert & 
                                                               tibble_exonic_parent_transcript_entries$end >= a1$alternative_exon_start %>% type.convert)
                                                     
                                                     if (row.indices_reference_exons_overlapping_with_alternative_exon %>% length == 0) {
                                                       return(NULL)
                                                     } else {
                                                       return(b1[row.indices_reference_exons_overlapping_with_alternative_exon, "exon_number"] %>% na.omit %>% unique %>% unlist)
                                                     }
                                                     
                                                   } )
                                                 
                                                 # calculate whether the presence of the differential exon instead of the original exons changes the frame.
                                                 
                                                 list_frameshift_info <- purrr::map2(.x = a1$list_tibbles_GTF_entries_parent_transcript,
                                                                                     .y = list_ref_exon_numbers_overlapping_alternative_exon,
                                                                                     .f = function(b1, b2) {
                                                                                       
                                                                                       # DEBUG ###
                                                                                       # b1 <- a1$list_tibbles_GTF_entries_parent_transcript %>% .[[1]]
                                                                                       # b2 <- list_ref_exon_numbers_overlapping_alternative_exon[[1]]
                                                                                       ########### 
                                                                                       
                                                                                       width_overlapped_exon <- b1 %>% dplyr::filter(type == "exon" & exon_number %in% b2) %>% .$width %>% sum
                                                                                       width_alternative_exon <- (a1$alternative_exon_end %>% type.convert) - (a1$alternative_exon_start %>% type.convert) + 1
                                                                                       
                                                                                       width_difference <- abs(width_overlapped_exon - width_alternative_exon)
                                                                                       
                                                                                       # if modulus 3 is 0, then it's frame-preserving, if not, then it's not frame preserving.
                                                                                       has_frameshift <- width_difference %% 3 != 0
                                                                                       
                                                                                       return(has_frameshift)
                                                                                       
                                                                                     } )
                                                 
                                                 # knock the differential exon into the reference transcript.
                                                 # BEHAVIOUR:
                                                 # 1. if there were no overlapping exons found, then leave alternative exon as-is.
                                                 # 2. if overlapping exons found, delete all overlapping entries.
                                                 # 3. Depending on stranded or increasing exon numbering of the GTF, we will then extract all entries with exon numbers from the start to just before the alternative exon, then the alternative exon.
                                                 # 4. now if the annotated stop codon is before the alternative exon or start codon is after the alternative exon, then we immediately know that it does not generate a poison exon.
                                                 
                                                 list_transcript_segments <- purrr::map2(.x = a1$list_tibbles_GTF_entries_parent_transcript,
                                                                                       .y = list_ref_exon_numbers_overlapping_alternative_exon,
                                                                                       .f = function(b1, b2) {
                                                                                         
                                                                                         # DEBUG ###
                                                                                         # b1 <- a1$list_tibbles_GTF_entries_parent_transcript %>% .[[1]]
                                                                                         # b2 <- list_ref_exon_numbers_overlapping_alternative_exon[[1]]
                                                                                         ###########
                                                                                         
                                                                                         tibble_exonic_parent_transcript_entries <- b1 %>% dplyr::filter(type == "exon")
                                                                                         vector_overlapped_exon_numbers <- b2
                                                                                         
                                                                                         # CHECK FOR STOP CODON POSITION
                                                                                         if ("stop_codon" %in% b1$type == TRUE) {
                                                                                           if (b1$strand %>% unique == "+" &
                                                                                               a1$alternative_exon_start %>% type.convert >= b1 %>% dplyr::filter(type == "stop_codon") %>% .$end %>% min - 2) {
                                                                                             return("alt_exon_is_in_threeprime_utr")
                                                                                           } else if (b1$strand %>% unique == "-" &
                                                                                                      a1$alternative_exon_end %>% type.convert <= b1 %>% dplyr::filter(type == "stop_codon") %>% .$start %>% max + 2) {
                                                                                             return("alt_exon_is_in_threeprime_utr")
                                                                                           }
                                                                                         }
                                                                                         
                                                                                         # CHECK FOR START CODON POSITION
                                                                                         if ("start_codon" %in% b1$type == TRUE) {
                                                                                           if (b1$strand %>% unique == "+" &
                                                                                               a1$alternative_exon_end %>% type.convert <= b1 %>% dplyr::filter(type == "start_codon") %>% .$start %>% max + 2) {
                                                                                             return("alt_exon_is_in_fiveprime_utr")
                                                                                           } else if (b1$strand %>% unique == "-" &
                                                                                                      a1$alternative_exon_start %>% type.convert >= b1 %>% dplyr::filter(type == "start_codon") %>% .$end %>% min - 2) {
                                                                                             return("alt_exon_is_in_fiveprime_utr")
                                                                                           }
                                                                                         }
                                                                                         
                                                                                         # CHECK FOR UNKNOWN REFERENCE STRAND
                                                                                         if (b1$strand %>% unique %in% c("+", "-") == FALSE) {
                                                                                           return("strand_unknown_in_reference")
                                                                                         }
                                                                                         
                                                                                         # create the transcript segment going from exon 1 to the alternative exon.
                                                                                         ## if positive strand: choose all exons which ALL have both start/end less than the alternative exon start
                                                                                         ## if negative strand: choose all exons which ALL have both start/end greater than the alternative exon end
                                                                                         if (b1$strand %>% unique == "+") {
                                                                                           # fetch row indices of upstream exons
                                                                                           row.indices_exons_upstream_of_alternative_exon <- which(tibble_exonic_parent_transcript_entries$start < a1$alternative_exon_start %>% type.convert & tibble_exonic_parent_transcript_entries$end < a1$alternative_exon_start %>% type.convert)
                                                                                           
                                                                                           # subset tibble based on row indices
                                                                                           tibble_upstream_exons <- tibble_exonic_parent_transcript_entries[row.indices_exons_upstream_of_alternative_exon, ]
                                                                                           # add the alternative exon as the entry
                                                                                           tibble_alternative_exon_entry <- tibble_exonic_parent_transcript_entries[1, ] %>% 
                                                                                             dplyr::mutate_at(.vars = "start", .funs = function(x) {return(a1$alternative_exon_start %>% type.convert)}) %>% 
                                                                                             dplyr::mutate_at(.vars = "end", .funs = function(x) {return(a1$alternative_exon_end %>% type.convert)}) %>% 
                                                                                             dplyr::mutate_at(.vars = "width", .funs = function(x) {.$end - .$start + 1}) %>% 
                                                                                             dplyr::mutate_at(.vars = c("source"), .funs = function(x) {return("poison_exon_finder")}) %>% 
                                                                                             dplyr::mutate_at(.vars = "type", .funs = function(x) {return("exon")}) %>% 
                                                                                             dplyr::mutate_at(.vars = colnames(tibble_exonic_parent_transcript_entries)[!colnames(tibble_exonic_parent_transcript_entries) %in% c("seqnames", "start", "end", "width", "strand", "source", "type", "score", "gene_id")], .funs = function(x) {return(NA)})
                                                                                           
                                                                                           tibble_transcript_segment <- dplyr::bind_rows(tibble_upstream_exons, tibble_alternative_exon_entry) %>% 
                                                                                                              # sort by increasing start positions
                                                                                                              dplyr::arrange(start) %>%
                                                                                                              # assign STRANDED exon numbers
                                                                                                              dplyr::mutate_at(.vars = "exon_number", .funs = function(x) {return(1:nrow(.))} )
                                                                                                            
                                                                                           } else if (b1$strand %>% unique == "-") {
                                                                                             # fetch row indices of upstream exons
                                                                                             row.indices_exons_upstream_of_alternative_exon <- which(tibble_exonic_parent_transcript_entries$start > a1$alternative_exon_end %>% type.convert & tibble_exonic_parent_transcript_entries$end > a1$alternative_exon_end %>% type.convert)
                                                                                             
                                                                                             # subset tibble based on row indices
                                                                                             tibble_upstream_exons <- tibble_exonic_parent_transcript_entries[row.indices_exons_upstream_of_alternative_exon, ]
                                                                                             # add the alternative exon as the entry
                                                                                             tibble_alternative_exon_entry <- tibble_exonic_parent_transcript_entries[1, ] %>% 
                                                                                               dplyr::mutate_at(.vars = "start", .funs = function(x) {return(a1$alternative_exon_start %>% type.convert)}) %>% 
                                                                                               dplyr::mutate_at(.vars = "end", .funs = function(x) {return(a1$alternative_exon_end %>% type.convert)}) %>% 
                                                                                               dplyr::mutate_at(.vars = "width", .funs = function(x) {.$end - .$start + 1}) %>% 
                                                                                               dplyr::mutate_at(.vars = c("source"), .funs = function(x) {return("poison_exon_finder")}) %>% 
                                                                                               dplyr::mutate_at(.vars = colnames(tibble_exonic_parent_transcript_entries)[!colnames(tibble_exonic_parent_transcript_entries) %in% c("seqnames", "start", "end", "width", "strand", "source", "type", "score", "gene_id")], .funs = function(x) {return(NA)})
                                                                                             
                                                                                             tibble_transcript_segment <- dplyr::bind_rows(tibble_upstream_exons, tibble_alternative_exon_entry) %>% 
                                                                                               # sort by DECREASING start positions
                                                                                               dplyr::arrange(desc(start)) %>%
                                                                                               # assign STRANDED exon numbers
                                                                                               dplyr::mutate_at(.vars = "exon_number", .funs = function(x) {return(1:nrow(.))} )
                                                                                             
                                                                                           } 
                                                                                         
                                                                                         return(tibble_transcript_segment)
                                                                                         
                                                                                       } )
                                                 
                                                 return(splice(a1,
                                                               "list_ref_exon_numbers_overlapping_alternative_exon" = list(list_ref_exon_numbers_overlapping_alternative_exon),
                                                               "has_frameshift" = list(list_frameshift_info),
                                                               "list_transcript_segments" = list(list_transcript_segments)))
                                                 
                                               }, .progress = TRUE, .options = future_options(globals = c("type.convert", "dplyr")))

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# drop elements which dont have any more valid transcript segments
list_modified_reference_transcript_segments_pruned <- list_modified_reference_transcript_segments %>% purrr::discard(.p = ~.x$list_transcript_segments %>% purrr::map(~length(.x) == 0 | length(.x) == 1) %>% unlist %>% all == TRUE)

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

cat("turn each transcript segment into vectors of genome-relative nucleotide positions.\n")
# also calculate the transcript-relative start codon position.
list_genome_relative_positions_of_transcript_segments <- future_imap(.x = list_modified_reference_transcript_segments_pruned,
                                                                     .f = function(a1, a2) {
                                                                       
                                                                       # DEBUG ###
                                                                       # a1 <- list_modified_reference_transcript_segments_pruned[[73]]
                                                                       ###########
                                                                       
                                                                       if (save_workspace_when_done == "DEBUG") {
                                                                         cat(a2, "\n")
                                                                       }
                                                                       
                                                                       # create genome-relative nucleotide positions of each transcript
                                                                       vector_transcript_segment_forward_genome_relative_coords <- purrr::map(
                                                                         .x = a1$list_transcript_segments, 
                                                                         .f = function(b1) {
                                                                           
                                                                           # DEBUG ###
                                                                           # b1 <- a1$list_transcript_segments %>% .[[1]]
                                                                           ###########
                                                                           
                                                                           if (b1 == "strand_unknown_in_reference") {
                                                                             return("strand_unknown_in_reference")
                                                                           } else if (b1 == "alt_exon_is_in_threeprime_utr") {
                                                                             return("alt_exon_is_in_threeprime_utr")
                                                                           } else if (b1 == "alt_exon_is_in_fiveprime_utr") {
                                                                             return("alt_exon_is_in_fiveprime_utr")
                                                                           }
                                                                           
                                                                           # include 2 more nucleotides past the VSR for translation, in order to catch reconstituted stop codons.
                                                                           vector_forward_coords_up_to_including_alt_exon <- purrr::map2(.x = b1$start, .y = b1$end, .f = ~.x:.y) %>% unlist %>% type.convert %>% sort
                                                                           
                                                                           if (b1$strand %>% unique == "+") {
                                                                             vector_forward_coords_for_translation <- c(vector_forward_coords_up_to_including_alt_exon, a1$VSR_end %>% type.convert + 1, a1$VSR_end %>% type.convert + 2)
                                                                           } else if (b1$strand %>% unique == "-") {
                                                                             vector_forward_coords_for_translation <- c(a1$VSR_start %>% type.convert - 2, a1$VSR_start %>% type.convert - 1, vector_forward_coords_up_to_including_alt_exon)
                                                                           }
                                                                           
                                                                           return(vector_forward_coords_for_translation)
                                                                           
                                                                         } )
                                                                       
                                                                       transcript_segment_relative_first_nt_of_start_codon_position <- purrr::map2(
                                                                         .x = vector_transcript_segment_forward_genome_relative_coords,
                                                                         .y = a1$list_tibbles_GTF_entries_parent_transcript,
                                                                         .f = function(b1, b2) {
                                                                           
                                                                           # DEBUG ###
                                                                           # b1 <- vector_transcript_segment_forward_genome_relative_coords[[3]]
                                                                           # b2 <- a1$list_tibbles_GTF_entries_parent_transcript %>% .[[3]]
                                                                           ###########
                                                                           
                                                                           if (b1 == "strand_unknown_in_reference") {
                                                                             return("strand_unknown_in_reference")
                                                                           } else if (b1 == "alt_exon_is_in_threeprime_utr") {
                                                                             return("alt_exon_is_in_threeprime_utr")
                                                                           } else if (b1 == "alt_exon_is_in_fiveprime_utr") {
                                                                             return("alt_exon_is_in_fiveprime_utr")
                                                                           }
                                                                           
                                                                           # get strand 
                                                                           current_strand <- b2$strand %>% unique
                                                                           
                                                                           # take cases.
                                                                           ## if + strand, then 1st nt of start codon is the min
                                                                           ## if - strand, then 1st nt of start codon is the max
                                                                           if (b2 %>% dplyr::filter(type == "start_codon") %>% dim %>% prod == 0) {
                                                                           } else if (current_strand == "+") {
                                                                             genome_relative_position_first_nt_of_start_codon <- b2 %>% dplyr::filter(type == "start_codon") %>% purrr::map2(.x = .$start, .y = .$end, .f = ~.x:.y) %>% unlist %>% type.convert %>% min
                                                                           } else if (current_strand == "-") {
                                                                             genome_relative_position_first_nt_of_start_codon <- b2 %>% dplyr::filter(type == "start_codon") %>% purrr::map2(.x = .$start, .y = .$end, .f = ~.x:.y) %>% unlist %>% type.convert %>% max
                                                                           }
                                                                           
                                                                           # get transcript-segment-relative (NOT FORWARD) positions of the first nt. of start codon
                                                                           ## if start codon didn't exist, then the transcript_segment_relative_position_first_nt_of_start_codon = 1
                                                                           if (b2 %>% dplyr::filter(type == "start_codon") %>% dim %>% prod == 0) {
                                                                             transcript_segment_relative_position_first_nt_of_start_codon <- 1 
                                                                           } else if (current_strand == "+") {
                                                                             transcript_segment_relative_position_first_nt_of_start_codon <- which(b1 == genome_relative_position_first_nt_of_start_codon)
                                                                           } else if (current_strand == "-") {
                                                                             transcript_segment_relative_position_first_nt_of_start_codon <- which((b1 %>% rev) == genome_relative_position_first_nt_of_start_codon)
                                                                           }
                                                                           
                                                                           if (length(transcript_segment_relative_position_first_nt_of_start_codon) == 0) {
                                                                             if (use_start_codon == "ALWAYS") {
                                                                               return("start_codon_removed_by_knock_in")
                                                                             } else {
                                                                               return(1)
                                                                             }
                                                                           }
                                                                           
                                                                           return(transcript_segment_relative_position_first_nt_of_start_codon)
                                                                           
                                                                         } )
                                                                       
                                                                       # if start codon was removed by knock-in, then we have to prune and update the entire list structure to reflect this.
                                                                       if (length(transcript_segment_relative_first_nt_of_start_codon_position) != 0) {
                                                                         element.indices_start_codon_removed <- integer(0) 
                                                                       }
                                                                       
                                                                       element.indices_start_codon_removed <- transcript_segment_relative_first_nt_of_start_codon_position %>% purrr::map(~.x == "start_codon_removed_by_knock_in") %>% unlist %>% which
                                                                       
                                                                       if (length(element.indices_start_codon_removed) > 0) {
                                                                         
                                                                         a1[c("list_tibbles_GTF_entries_matched_to_VSR", "list_tibbles_GTF_entries_parent_transcript", "list_logical_start_codon_present", "list_ref_exon_numbers_overlapping_alternative_exon", "has_frameshift", "list_transcript_segments")] <- a1[c("list_tibbles_GTF_entries_matched_to_VSR", "list_tibbles_GTF_entries_parent_transcript", "list_logical_start_codon_present", "list_ref_exon_numbers_overlapping_alternative_exon", "has_frameshift", "list_transcript_segments")] %>% purrr::map(~.x[-element.indices_start_codon_removed])
                                                                         
                                                                         vector_transcript_segment_forward_genome_relative_coords <- vector_transcript_segment_forward_genome_relative_coords[-element.indices_start_codon_removed]
                                                                         transcript_segment_relative_first_nt_of_start_codon_position <- transcript_segment_relative_first_nt_of_start_codon_position[-element.indices_start_codon_removed]
                                                                         
                                                                       }
                                                                       
                                                                       # if there are no more matched GTF entries remaining after prune, then exit the loop immediately.
                                                                       if (length(a1$list_tibbles_GTF_entries_matched_to_VSR) == 0) {
                                                                         return(a1)
                                                                       } else {
                                                                         
                                                                         transcript_segment_relative_first_nt_of_alternative_exon <- purrr::map2(
                                                                           .x = vector_transcript_segment_forward_genome_relative_coords,
                                                                           .y = a1$list_tibbles_GTF_entries_parent_transcript,
                                                                           .f = function(b1, b2) {
                                                                             
                                                                             # DEBUG ###
                                                                             # b1 <- vector_transcript_segment_forward_genome_relative_coords[[1]]
                                                                             # b2 <- a1$list_tibbles_GTF_entries_parent_transcript %>% .[[1]]
                                                                             ###########
                                                                             
                                                                             if (b1 == "strand_unknown_in_reference") {
                                                                               return("strand_unknown_in_reference")
                                                                             } else if (b1 == "alt_exon_is_in_threeprime_utr") {
                                                                               return("alt_exon_is_in_threeprime_utr")
                                                                             } else if (b1 == "alt_exon_is_in_fiveprime_utr") {
                                                                               return("alt_exon_is_in_fiveprime_utr")
                                                                             }
                                                                             
                                                                             # get strand 
                                                                             current_strand <- b2$strand %>% unique
                                                                             
                                                                             # get transcript-segment-relative (NOT FORWARD) positions of the first nt. of start codon
                                                                             ## if start codon didn't exist, then the transcript_segment_relative_position_first_nt_of_start_codon = 1
                                                                             if (current_strand == "+") {
                                                                               transcript_segment_relative_position_first_nt_of_alternative_exon <- which(b1 == a1$alternative_exon_start %>% type.convert)
                                                                             } else if (current_strand == "-") {
                                                                               transcript_segment_relative_position_first_nt_of_alternative_exon <- which((b1 %>% rev) == a1$alternative_exon_end %>% type.convert)
                                                                             }
                                                                             
                                                                             return(transcript_segment_relative_position_first_nt_of_alternative_exon)
                                                                             
                                                                           } )
                                                                         
                                                                         return(splice(a1,
                                                                                       "vector_transcript_segment_forward_genome_relative_coords" = list(vector_transcript_segment_forward_genome_relative_coords),
                                                                                       "transcript_segment_relative_first_nt_of_start_codon_position" = list(transcript_segment_relative_first_nt_of_start_codon_position),
                                                                                       "transcript_segment_relative_first_nt_of_alternative_exon" = list(transcript_segment_relative_first_nt_of_alternative_exon)))
                                                                         
                                                                       }
                                                                       
                                                                     }, .progress = TRUE, .options = future_options(globals = c("type.convert", "dplyr", "use_start_codon")) )

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# prune L1 elements which don't have VSR matches remaining in the event that the start codon was knocked out.
list_genome_relative_positions_of_transcript_segments_pruned <- list_genome_relative_positions_of_transcript_segments %>% purrr::discard(.p = ~length(.x$list_tibbles_GTF_entries_matched_to_VSR) == 0)

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

cat("retrieve the forward nucleotides of the transcript segment\n")
## loop thru each chromosome. first, we split the list by chromosome
list_genome_relative_positions_of_transcript_segments_split_by_chr <- purrr::map(.x = chr_to_run,
                                                                                 .f = function(a1) {
                                                                                   
                                                                                   list_genome_relative_positions_of_transcript_segments_pruned %>% 
                                                                                     purrr::keep(.p = ~.x$chr == a1) %>% 
                                                                                     return
                                                                                   
                                                                                 } ) %>% set_names(chr_to_run) %>% purrr::compact()

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

cat("loop thru the split list\n")
list_genome_forward_genomic_nucleotides <- future_map2(.x = list_genome_relative_positions_of_transcript_segments_split_by_chr,
                                                                     .y = names(list_genome_relative_positions_of_transcript_segments_split_by_chr),
                                                                     .f = function(a1, a2) {
                                                                       
                                                                       # DEBUG ###
                                                                       # a1 <- list_genome_relative_positions_of_transcript_segments_split_by_chr[[2]]
                                                                       # a2 <- list_genome_relative_positions_of_transcript_segments_split_by_chr[2] %>% names
                                                                       ###########
                                                                       
                                                                       # get positions of the vector where the path of the ref. genome fasta
                                                                       vector_ref_genome_paths_by_chr_position <- grep(x = vector_ref_genome_paths_by_chr, pattern = paste("(\\D|^)", a2, ".fa$", sep = ""))
                                                                       
                                                                       if (length(vector_ref_genome_paths_by_chr_position) != 1) {
                                                                         
                                                                         stop("Something is wrong with the contents of the fasta file directory. Please check that it's structured in the desired format.")
                                                                         }
                                                                       
                                                                       # temporary allocation to ref genome fasta list
                                                                       reference_genome_fasta_chr_temp <- seqinr::read.fasta(file = paste(vector_ref_genome_paths_by_chr[vector_ref_genome_paths_by_chr_position]), forceDNAtolower = FALSE)
                                                                       # map over each exon, splice in the fwd. nucleotides
                                                                       output_exon_list <- purrr::map(
                                                                         .x = a1,
                                                                         .f = function(b1) {
                                                                           
                                                                           # DEBUG ###
                                                                           # b1 <- a1[[1]]
                                                                           ###########
                                                                           
                                                                           list_forward_nucleotides <- purrr::map(
                                                                             .x = b1$vector_transcript_segment_forward_genome_relative_coords,
                                                                             .f = function(c1) {
                                                                               
                                                                               # DEBUG ###
                                                                               # c1 <- b1$vector_transcript_segment_forward_genome_relative_coords %>% .[[1]]
                                                                               ###########
                                                                               
                                                                               if (c1 %in% c("strand_unknown_in_reference", "alt_exon_is_in_threeprime_utr", "alt_exon_is_in_fiveprime_utr")) {return(NA)}                                                          
                                                                               reference_genome_fasta_chr_temp[[b1$chr]][c1] %>% return
                                                                               
                                                                             } )
                                                                           
                                                                           return(splice(b1,
                                                                                         "list_forward_nucleotides" = list(list_forward_nucleotides)))
                                                                           
                                                                         } ) %>% return 
                                                                       
                                                                     }, .progress = TRUE, .options = future_options(globals = c("vector_ref_genome_paths_by_chr", "seqinr::read.fasta"))) %>% flatten

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

cat("three-frame translation\n")
# after translation, find location of stop codon
# if there are stop codons in the alternative exon, then it's a poison exon (introduces PTC)
list_three_frame_translation <- future_imap(.x = list_genome_forward_genomic_nucleotides, 
                                           .f = function(a1, a2) {
                                             
                                             # DEBUG ###
                                             # a1 <- list_genome_relative_positions_of_transcript_segments[[10]]
                                             # a1 <- test[[1]]
                                             ##########
                                             
                                             if (save_workspace_when_done == "DEBUG") {
                                               cat(a2, "\n")
                                             }
                                             
                                             list_three_frame_translation <- purrr::pmap(
                                               .l = list("b1" = a1$list_forward_nucleotides,
                                                         "b2" = a1$list_tibbles_GTF_entries_parent_transcript,
                                                         "b3" = a1$transcript_segment_relative_first_nt_of_start_codon_position,
                                                         "b4" = a1$list_logical_start_codon_present),
                                               .f = function(b1, b2, b3, b4) {
                                                 
                                                 # DEBUG ###
                                                 # b1 <- a1$list_forward_nucleotides %>% .[[7]]
                                                 # b2 <- a1$list_tibbles_GTF_entries_parent_transcript %>% .[[7]]
                                                 # b3 <- a1$transcript_segment_relative_first_nt_of_start_codon_position %>% .[[7]]
                                                 # b4 <- a1$list_logical_start_codon_present %>% .[[7]]
                                                 ###########
                                                 
                                                 transcript_segment_forward_nucleotides <- b1
                                                 if (is.na(transcript_segment_forward_nucleotides) == TRUE) {return(NA)}
                                                 
                                                 transcript_segment_relative_first_nt_of_start_codon_position <- b3
                                                   
                                                 # get the current strand
                                                 current_strand <- b2$strand %>% unique
                                                 
                                                 # subset the forward nucleotides from the start of the start codon
                                                 if (is.na(transcript_segment_forward_nucleotides) == TRUE) {
                                                   return(NA)
                                                   # list("translation_frame_0" = "*",
                                                   #      "translation_frame_1" = "*",
                                                   #      "translation_frame_2" = "*")
                                                 } else if (current_strand == "+") {
                                                   list_raw_3FT_results <- nt.sequence_strand_threeframetranslate(vector_forward_nucleotides = transcript_segment_forward_nucleotides[transcript_segment_relative_first_nt_of_start_codon_position:length(transcript_segment_forward_nucleotides)], strand = "+")
                                                   # if start codon was present, then all frames = 0th frame
                                                   if (b4 == TRUE) {
                                                     list_raw_3FT_results <- list_raw_3FT_results %>% purrr::map(.f = ~list_raw_3FT_results[[1]])
                                                   }
                                                 } else if (current_strand == "-") {
                                                   list_raw_3FT_results <- nt.sequence_strand_threeframetranslate(vector_forward_nucleotides = transcript_segment_forward_nucleotides %>% rev %>% .[transcript_segment_relative_first_nt_of_start_codon_position:length(transcript_segment_forward_nucleotides)] %>% seqinr::comp(seq = ., forceToLower = FALSE), strand = "+")
                                                   # if start codon was present, then all frames = 0th frame
                                                   if (b4 == TRUE) {
                                                     list_raw_3FT_results <- list_raw_3FT_results %>% purrr::map(.f = ~list_raw_3FT_results[[1]])
                                                   }
                                                 }
                                                 
                                                 return(list_raw_3FT_results)
                                                 
                                               } ) # L2
                                             
                                             return(splice(a1,
                                                           "list_three_frame_translation" = list(list_three_frame_translation)))
                                             
                                           }, .progress = TRUE, .options = future_options(globals = c("vector_ref_genome_paths_by_chr", "seqinr::read.fasta", "nt.sequence_strand_threeframetranslate")) )

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

cat("finally test for the presence of a stop codon in the alternative exon.\n")
list_test_for_PTC <- future_imap(.x = list_three_frame_translation,
                                .f = function(a1, a2) {
                                  
                                  # DEBUG ###
                                  # a1 <- list_three_frame_translation[[4]]
                                  ##########
                                  
                                  cat(a2, "\n")
                                  
                                  # create a list of PTC testing information:
                                  ## coords of stop 
                                  list_PTC_testing_info <- purrr::pmap(
                                    .l = list("b1" = a1$list_three_frame_translation,
                                              "b2" = a1$transcript_segment_relative_first_nt_of_start_codon_position,
                                              "b3" = a1$transcript_segment_relative_first_nt_of_alternative_exon,
                                              "b4" = a1$vector_transcript_segment_forward_genome_relative_coords,
                                              "b5" = a1$list_tibbles_GTF_entries_parent_transcript),
                                    .f = function(b1, b2, b3, b4, b5) {
                                      
                                      # DEBUG ###
                                      # b1 <- a1$list_three_frame_translation %>% .[[1]]
                                      # b2 <- a1$transcript_segment_relative_first_nt_of_start_codon_position %>% .[[1]]
                                      # b3 <- a1$transcript_segment_relative_first_nt_of_alternative_exon %>% .[[1]]
                                      # b4 <- a1$vector_transcript_segment_forward_genome_relative_coords %>% .[[1]]
                                      # b5 <- a1$list_tibbles_GTF_entries_parent_transcript %>% .[[1]]
                                      ###########
                                      
                                      current_strand <- b5$strand %>% unique
                                      
                                      # retrieve all translation-relative coords of all stop codons.
                                      list_translation_relative_stop_codon_positions <- b1 %>% purrr::map(~which(.x == "*"))
                                      
                                      # convert translation-relative coords of all stop codons to transcript-segment-relative coords.
                                      ## only take the last nt. for simplicity of deteremining whether it lies in the alternative exon.
                                      ## 3x-2 and 3x
                                      list_transcript_segment_relative_last_nt_of_stop_codon_positions <- purrr::map(
                                        .x = list_translation_relative_stop_codon_positions,
                                        .f = function(c1) {
                                          if (c1 %>% length == 0) {
                                            return(NA)
                                          } else {
                                            return(c1*3 + b2 - 1)
                                          }
                                        } )
                                      
                                      # finally test for PTC in the alternative exon
                                      list_logical_PTC_exists_in_alternative_exon <- purrr::map(
                                        .x = list_transcript_segment_relative_last_nt_of_stop_codon_positions,
                                        .f = function(c1) {
                                          if (is.na(c1)) {
                                            return(FALSE)
                                          } else {
                                            return(c1 >= b3)
                                          }
                                        } )
                                      
                                      # return list of genomic coords of stop codons
                                      list_genome_relative_coords_of_stop_codons <- purrr::map2(
                                        .x = list_transcript_segment_relative_last_nt_of_stop_codon_positions,
                                        .y = list_logical_PTC_exists_in_alternative_exon,
                                        .f = function(c1, c2) {
                                          if (all(c2 == FALSE) == TRUE) {
                                            return(NA)
                                          } else {
                                           # depending on strand, look-up the genome-relative forward coords of each transcript segment.
                                             if (current_strand == "+") {
                                               vector_all_genome_relative_stop_codon_positions <- b4[c(c1, c1 - 1, c1 - 2)] %>% sort
                                               return(vector_all_genome_relative_stop_codon_positions)
                                             } else if (current_strand == "-") {
                                               vector_all_genome_relative_stop_codon_positions <- b4 %>% rev %>% .[c(c1, c1 - 1, c1 - 2)] %>% sort(decreasing = TRUE)
                                               return(vector_all_genome_relative_stop_codon_positions)
                                             }
                                          }
                                        } ) # L3
                                      
                                      return(list(
                                        "list_transcript_segment_relative_last_nt_of_stop_codon_positions" = list_transcript_segment_relative_last_nt_of_stop_codon_positions,
                                        "list_logical_PTC_exists_in_alternative_exon" = list_logical_PTC_exists_in_alternative_exon,
                                        "list_genome_relative_coords_of_stop_codons" = list_genome_relative_coords_of_stop_codons,
                                        "matched_strand" = current_strand
                                        ))
                                      
                                    } ) # L2
                                  
                                  return(splice(a1,
                                                "list_transcript_segment_relative_last_nt_of_stop_codon_positions" = list_PTC_testing_info %>% purrr::map(~.x$list_transcript_segment_relative_last_nt_of_stop_codon_positions) %>% list,
                                                "list_logical_PTC_exists_in_alternative_exon" = list_PTC_testing_info %>% purrr::map(~.x$list_logical_PTC_exists_in_alternative_exon) %>% list,
                                                "list_genome_relative_coords_of_stop_codons" = list_PTC_testing_info %>% purrr::map(~.x$list_genome_relative_coords_of_stop_codons) %>% list))
                                  
                                }, .progress = TRUE, .options = future_options(globals = FALSE )) # L1

if (save_workspace_when_done == "DEBUG") {
  save(list_test_for_PTC, file = paste(output_dir, "/", output_name, ".Rlist", sep = ""))
}

# last prune for elements without any valid forward genome coords.
list_test_for_PTC_pruned <- list_test_for_PTC %>% purrr::discard(.p = ~.x$list_forward_nucleotides %>% purrr::map(~is.na(.x)) %>% unlist %>% all == TRUE)

cat("percolate and tibblise the results.\n")
list_results <- future_map2(.x = list_test_for_PTC_pruned,
                            .y = 1:length(list_test_for_PTC_pruned),
                            .f = function(a1, a2) {
                              
                              # DEBUG ###
                              # a1 <- list_test_for_PTC[[1]]
                              ##########
                              
                              cat(a2, "\n")
                              
                              if (save_workspace_when_done == "DEBUG") {
                                cat(a2, "\n")
                              }
                              
                              # test line
                              # a1$list_logical_PTC_exists_in_alternative_exon %>% purrr::map_depth(.depth = 2, .f = ~.x %>% paste(collapse = ",")) %>% purrr::map(~ .x %>% as_tibble %>% t %>% (function(x) {tibble <- x; colnames(tibble) <- c("PTC_exists_in_alternative_exon"); return(tibble)} ) %>% as_tibble(rownames = "translation_frame"))
                              
                              # tibblise...
                              ## it seems the easiest to go from the roots-up, since we collect the loosest ends first and work our way back up to a common format.
                              
                              # concatenate the three-frame translation results
                              a1 <- a1 %>% purrr::map_at(.at = "list_three_frame_translation", .f = ~purrr::map_depth(.x = .x, .depth = 2, .f = ~.x %>% paste(collapse = "")))
                              # concatenate the three-frame translation-related elements
                              a1 <- a1 %>% purrr::map_at(.at = c("list_transcript_segment_relative_last_nt_of_stop_codon_positions", "list_logical_PTC_exists_in_alternative_exon", "list_genome_relative_coords_of_stop_codons"), .f = ~purrr::map_depth(.x = .x, .depth = 2, .f = ~.x %>% paste(collapse = ",")))
                              
                              # tibblise the three-frame translation elements
                              list_three_frame_translation_elements0 <- purrr::map2(.x = a1[c("list_three_frame_translation", "list_transcript_segment_relative_last_nt_of_stop_codon_positions", "list_logical_PTC_exists_in_alternative_exon", "list_genome_relative_coords_of_stop_codons")],
                                                                                    .y = c("list_three_frame_translation", "transcript_segment_relative_last_nt_of_stop_codon_positions", "PTC_exists_in_alternative_exon", "genome_relative_coords_of_stop_codons"),
                                                                                    .f = function(b1, b2) {
                                                                                      
                                                                                      # DEBUG ###
                                                                                      # b1 <- a1[c("list_three_frame_translation", "list_transcript_segment_relative_last_nt_of_stop_codon_positions", "list_logical_PTC_exists_in_alternative_exon", "list_genome_relative_coords_of_stop_codons")] %>% .[[1]]
                                                                                      # b2 <- c("three_frame_translation_result", "transcript_segment_relative_last_nt_of_stop_codon_positions", "PTC_exists_in_alternative_exon", "genome_relative_coords_of_stop_codons") %>% .[[1]]
                                                                                      ##########
                                                                                      
                                                                                      list_tibblised_elements <- b1 %>% purrr::discard(.p = ~length(.x) == 1) %>% purrr::map(~.x %>% as_tibble %>% t %>% (function(x) {tibble <- x; colnames(tibble) <- b2; return(tibble)} ) %>% as_tibble(rownames = "translation_frame"))
                                                                                      
                                                                                      return(list_tibblised_elements)
                                                                                      
                                                                                    } )
                              
                              # table join
                              list_three_frame_translation_elements1 <- list_three_frame_translation_elements0 %>% purrr::pmap(.l = ., .f = ~list(...) %>% purrr::reduce(dplyr::left_join))
                              # add matched transcript_id as a column for each element, then rbind, tibblise
                              tibble_three_frame_translation_elements <- purrr::map2(.x = list_three_frame_translation_elements1,
                                                                                     .y = names(list_three_frame_translation_elements1), 
                                                                                     .f = ~.x %>% add_column("matched_transcript_id" = .y)) %>% rbindlist(use.names = TRUE) %>% as_tibble
                              
                              # now tibblise each matched transcript element
                              list_matched_transcript_elements0 <- purrr::map2(.x = a1[c("list_logical_start_codon_present", "list_ref_exon_numbers_overlapping_alternative_exon", "has_frameshift", "vector_transcript_segment_forward_genome_relative_coords", "transcript_segment_relative_first_nt_of_start_codon_position", "transcript_segment_relative_first_nt_of_alternative_exon", "list_forward_nucleotides")],
                                                                               .y = c("start_codon_present", "ref_exon_numbers_overlapping_alternative_exon", "has_frameshift", "vector_transcript_segment_forward_genome_relative_coords", "transcript_segment_relative_first_nt_of_start_codon_position", "transcript_segment_relative_first_nt_of_alternative_exon", "list_forward_nucleotides"),
                                                                               .f = function(b1, b2) {
                                                                                 
                                                                                 # DEBUG ###
                                                                                 # b1 <- a1[c("list_logical_start_codon_present", "list_ref_exon_numbers_overlapping_alternative_exon", "has_frameshift", "vector_transcript_segment_forward_genome_relative_coords", "transcript_segment_relative_first_nt_of_start_codon_position", "transcript_segment_relative_first_nt_of_alternative_exon", "list_forward_nucleotides")] %>% .[[4]]
                                                                                 # b2 <- c("start_codon_present", "ref_exon_numbers_overlapping_alternative_exon", "has_frameshift", "vector_transcript_segment_forward_genome_relative_coords", "transcript_segment_relative_first_nt_of_start_codon_position", "transcript_segment_relative_first_nt_of_alternative_exon", "list_forward_nucleotides") %>% .[[4]]
                                                                                 ##########
                                                                                 
                                                                                 list_tibblised_elements <- b1 %>% purrr::map(~.x %>% paste(collapse = ",")) %>% as_tibble %>% t %>% (function(x) {tibble <- x; colnames(tibble) <- b2; return(tibble)} ) %>% as_tibble(rownames = "matched_transcript_id")
                                                                                 
                                                                               } )
                              
                              # join each element onto tibble_three_frame_translation_elements by matched_transcript_id
                              tibble_3FT_and_transcript_elements_joined <- purrr::splice(tibble_three_frame_translation_elements, list_matched_transcript_elements0) %>% purrr::reduce(dplyr::left_join)
                              
                              # finally tibblise with the exon/VSR information
                              tibble_all_elements_combined <- purrr::splice(a1[c("chr", "VSR_start", "VSR_end", "alternative_exon_start", "alternative_exon_end", "strand", "splicemode")], tibble_3FT_and_transcript_elements_joined) %>% flatten %>% as_tibble
                              
                              return(tibble_all_elements_combined)
                              
                            }, .progress = TRUE, .options = future_options(globals = c("add_column", "as_tibble", "rbindlist")) )

if (save_workspace_when_done == "YES" | save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# rbind and tibblise
tibble_results <- list_results %>% rbindlist(use.names = TRUE) %>% as_tibble

# write tibble
write.table(tibble_results, paste(output_dir, "/", output_name, ".txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)

# if user has specified that FASTA needs to be written (i.e. this script is to be repurposed for three-frame translation of exons), then continue working.

if (output_fasta == TRUE) {
  
  cat("FASTA output mode specified. now running 3-frame translation\n")
  
  cat("get exon-only virtual AA sequence\n")
  # annotate each artificial transcript entry with width preceding exon, and alternative exon width.
  # we want this because we're going to do ES = ceiling(x/3) + 1, EE = floor(y/3). 
  # x = width_preceding_alternative_exon - width before start codon - translation frame
  # y = width_preceding_alternative_exon + width_of_alternative_exon - width before start codon - translation frame
  list_3FT_result_exon_only <- future_map2(.x = list_test_for_PTC_pruned, 
                                           .y = 1:length(list_test_for_PTC_pruned),
                                           .f = function(a1, a2) {
                                             
                                             # DEBUG ###
                                             # a1 <- list_test_for_PTC_pruned[[33]]
                                             ###########
                                             
                                             cat(a2, "\n")
                                             
                                             list_filtered_3FT_results_alternative_exon_only <- purrr::pmap(.l = list(
                                               "b1" = a1$list_transcript_segments,
                                               "b2" = a1$transcript_segment_relative_first_nt_of_start_codon_position,
                                               "b3" = a1$list_three_frame_translation,
                                               "b4" = a1$list_logical_start_codon_present),
                                               .f = function(b1, b2, b3, b4) {
                                                 
                                                 # DEBUG ###
                                                 # b1 <- a1$list_transcript_segments %>% .[[1]]
                                                 # b2 <- a1$transcript_segment_relative_first_nt_of_start_codon_position %>% .[[1]]
                                                 # b3 <- a1$list_three_frame_translation %>% .[[1]]
                                                 # b4 <- a1$list_logical_start_codon_present %>% .[[1]]
                                                 ###########
                                                 
                                                 if (b1 %>% data.class == "character") {
                                                   return(b1)
                                                 } else {
                                                   
                                                   # retrieve widths of all exons barring the knocked-in exon.
                                                   width_preceding_alternative_exon <- b1 %>% 
                                                     dplyr::filter(source != "poison_exon_finder" & type == "exon") %>%
                                                     .$width %>% sum
                                                   
                                                   # retrieve width of the query exon
                                                   width_of_alternative_exon <- b1 %>% 
                                                     dplyr::filter(source == "poison_exon_finder" & type == "exon") %>%
                                                     .$width %>% sum
                                                   
                                                   # calculate ES and EE for each translation frame
                                                   # if start codon was present, we just have to take the 0th frame to calculate ES and EE. 
                                                   # if start codon was not present, we consider all frames.
                                                   # if the start codon lies in the alternative exon, then taking the max between the start codon first nt position and the alternative exon position (1) means we protect ourselves from negative indices.
                                                   if (b4 == FALSE) {
                                                     list_ES <- c(0:2) %>% purrr::map(.f = ~max(ceiling((width_preceding_alternative_exon - b2 + 1 - .x) / 3) + 1, 1))
                                                     list_EE <- c(0:2) %>% purrr::map(.f = ~floor((width_preceding_alternative_exon + width_of_alternative_exon - b2 + 1 - .x) / 3))
                                                   } else if (b4 == TRUE) {
                                                     list_ES <- c(0, 0, 0) %>% purrr::map(.f = ~max(ceiling((width_preceding_alternative_exon - b2 + 1 - .x) / 3) + 1, 1))
                                                     list_EE <- c(0, 0, 0) %>% purrr::map(.f = ~floor((width_preceding_alternative_exon + width_of_alternative_exon - b2 + 1 - .x) / 3))
                                                   }
                                                   
                                                   # filter the raw three-frame translation results according to ES and EE
                                                   three_frame_translation_result_alternative_exon_only <- purrr::pmap(.l = list(
                                                     "c1" = list_ES,
                                                     "c2" = list_EE,
                                                     "c3" = b3),
                                                     .f = function(c1, c2, c3) {
                                                       
                                                       # DEBUG ###
                                                       # c1 <- list_ES %>% .[[1]]
                                                       # c2 <- list_EE %>% .[[1]]
                                                       # c3 <- b3 %>% .[[1]]
                                                       ###########
                                                       
                                                       c3[c1:c2] %>% return
                                                       
                                                       # if there is start codon present, then only return the uORF, and cancel the dORF.
                                                       if (b4 == TRUE) {
                                                         valid_uORF <- find_valid_uORF(AA_sequence = c3, exon_start_AA_position = c1, exon_end_AA_position = c2)
                                                         valid_dORF <- "NONE_VALID"
                                                       } else {
                                                         valid_uORF <- find_valid_uORF(AA_sequence = c3, exon_start_AA_position = c1, exon_end_AA_position = c2)
                                                         valid_dORF <- find_valid_dORF(AA_sequence = c3, exon_start_AA_position = c1, exon_end_AA_position = c2)
                                                       }
                                                      
                                                       return(list("valid_uORF" = valid_uORF,
                                                              "valid_dORF" = valid_dORF)
                                                              )
                                                        
                                                     } 
                                                   ) %>% set_names(names(b3)) # L3
                                                   
                                                   # extract the uORF and dORF results
                                                   list_uORF_results <- three_frame_translation_result_alternative_exon_only %>% purrr::map(~.x$valid_uORF)
                                                   list_dORF_results <- three_frame_translation_result_alternative_exon_only %>% purrr::map(~.x$valid_dORF)
                                                   
                                                   return(list(
                                                     "width_preceding_alternative_exon" = width_preceding_alternative_exon,
                                                     "width_of_alternative_exon" = width_of_alternative_exon,
                                                     "list_uORF_results" = list_uORF_results,
                                                     "list_dORF_results" = list_dORF_results
                                                   ))
                                                   
                                                 }
                                                 
                                               } ) # L2
                                             
                                             # concatenate the three-frame translation results
                                             a1 <- a1 %>% purrr::map_at(.at = "list_three_frame_translation", .f = ~purrr::map_depth(.x = .x, .depth = 2, .f = ~.x %>% paste(collapse = "")))
                                             # concatenate the three-frame translation-related elements
                                             a1 <- a1 %>% purrr::map_at(.at = c("list_transcript_segment_relative_last_nt_of_stop_codon_positions", "list_logical_PTC_exists_in_alternative_exon", "list_genome_relative_coords_of_stop_codons"), .f = ~purrr::map_depth(.x = .x, .depth = 2, .f = ~.x %>% paste(collapse = ",")))
                                             # concatenate the individual coords (per transcript)
                                             a1 <- a1 %>% purrr::map_at(.at = c("vector_transcript_segment_forward_genome_relative_coords", "list_forward_nucleotides"), .f = ~purrr::map_depth(.x = .x, .depth = 1, .f = ~.x %>% paste(collapse = ",")))
                                             
                                             return(purrr::splice(a1,
                                                                  "width_preceding_alternative_exon" = list_filtered_3FT_results_alternative_exon_only %>% purrr::map(.f = function(.x) {if (.x %>% data.class == "character") {.x %>% return} else {.x$width_preceding_alternative_exon %>% return}} ) %>% list,
                                                                  "width_of_alternative_exon" = list_filtered_3FT_results_alternative_exon_only %>% purrr::map(.f = function(.x) {if (.x %>% data.class == "character") {.x %>% return} else {.x$width_of_alternative_exon %>% return}} ) %>% list,
                                                                  "list_uORF_results" = list_filtered_3FT_results_alternative_exon_only %>% purrr::map(.f = function(.x) {if (.x %>% data.class == "character") {.x %>% return} else {.x$list_uORF_results %>% return}} ) %>% list,
                                                                  "list_dORF_results" = list_filtered_3FT_results_alternative_exon_only %>% purrr::map(.f = function(.x) {if (.x %>% data.class == "character") {.x %>% return} else {.x$list_dORF_results %>% return}} ) %>% list))
                                             
                                          }, .progress = TRUE, .options = future_options(globals = c("dplyr::filter"))) # L1
  
  if (save_workspace_when_done == "DEBUG") {
    save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
  }
  
  cat("percolate, rbind and tibblise\n")
  list_summarised_3FT_result_exon_only <- future_map2(.x = list_3FT_result_exon_only,
                                                      .y = 1:length(list_3FT_result_exon_only),
                                                      .f = function(a1, a2) {
                                                        
                                                        # DEBUG ###
                                                        # a1 <- list_3FT_result_exon_only[[157]]
                                                        ###########
                                                        
                                                        cat(a2, "\n")
                                                        
                                                        # collapse 3FT results into a string
                                                        list_3FT_alternative_exon_only_collapsed <- a1[c("list_uORF_results", "list_dORF_results")] %>% purrr::map_depth(.depth = 3, .f = ~.x %>% paste(collapse = ""))
                                                        # tibblise
                                                        tibble_3FT_alternative_exon_only <- purrr::map(
                                                          .x = list_3FT_alternative_exon_only_collapsed,
                                                          .f = function(b1) {
                                                            
                                                            # DEBUG ###
                                                            # b1 <- list_3FT_alternative_exon_only_collapsed[[1]]
                                                            ###########
                                                            
                                                            tibble_3FT_result <- purrr::map2(.x = b1,
                                                                                           .y = names(b1),
                                                                                           .f = ~.x %>% as_tibble(.name_repair = "unique") %>% t %>% as_tibble(rownames = "translation_frame", .name_repair = "unique") %>% setNames(c("translation_frame", "virtual_peptide_sequence")) %>% add_column("matched_transcripts" = .y)) %>% 
                                                                                             rbindlist %>% as_tibble %>% 
                                                                                             dplyr::distinct(virtual_peptide_sequence, matched_transcripts, .keep_all = TRUE)
                                                            
                                                            # splice in the coordinate information
                                                            tibble_coordinate_info <- purrr::map2(.x = a1[c("vector_transcript_segment_forward_genome_relative_coords", "list_forward_nucleotides")], .y = c("vector_transcript_segment_forward_genome_relative_coords", "list_forward_nucleotides"), .f = ~.x %>% as_tibble(.name_repair = "unique") %>% t %>% as_tibble(.name_repair = "unique", rownames = "matched_transcripts") %>% setNames(c("matched_transcripts", .y))) %>% purrr::reduce(dplyr::full_join, by = "matched_transcripts")
                                                            
                                                            # add the coordinate data by joining using the matched transcript names
                                                            tibble_3FT_result_with_coordinate_data <- dplyr::left_join(tibble_3FT_result, tibble_coordinate_info, by = "matched_transcripts")
                                                            
                                                            # collapse the matched transcripts
                                                            # tibble_3FT_result_with_coordinate_data <- tibble_3FT_result_with_coordinate_data %>%
                                                            #   dplyr::mutate("matched_transcripts" = matched_transcripts %>% paste(collapse = ",")) %>%
                                                            #   dplyr::distinct(virtual_peptide_sequence, .keep_all = TRUE)
                                                            
                                                            return(tibble_3FT_result_with_coordinate_data)
                                                            
                                                          } ) %>% rbindlist %>% as_tibble
                                                          
                                                          # if three frame translation was unsuccessful, then we deal with it
                                                        tibble_3FT_alternative_exon_only[tibble_3FT_alternative_exon_only$virtual_peptide_sequence %in% c("alt_exon_is_in_threeprime_utr", "alt_exon_is_in_fiveprime_utr", "strand_unknown_in_reference", "NONE_VALID"), "translation_frame"] <- "untranslatable"
                                                        
                                                        final_identifier <- if (a1$custom_identifier %>% is.na != TRUE) {a1$custom_identifier
                                                        } else {
                                                          paste(a1$chr, ":", a1$VSR_start %>% type.convert, "-", a1$VSR_end %>% type.convert, 
                                                                ";", 
                                                                a1$chr, ":", a1$alternative_exon_start %>% type.convert, "-", a1$alternative_exon_end%>% type.convert, 
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
                                                                              tibble_3FT_alternative_exon_only$matched_transcripts %>% unique %>% paste(collapse = ","), 
                                                                              " OS=",
                                                                              a1$organism, 
                                                                              " GN=",
                                                                              a1$gene_name, sep = "")
                                                        
                                                        # add fasta header to the alternative exon-only 3FT results
                                                        tibble_3FT_alternative_exon_only_with_fasta_header <- tibble_3FT_alternative_exon_only %>% 
                                                          add_column("fasta_header" = fasta_header,
                                                                     "final_identifier" = final_identifier)
                                                        
                                                        # finalise by adding in the input data
                                                        final_tibble <- purrr::splice(
                                                          a1[c("chr", "VSR_start", "VSR_end", "alternative_exon_start", "alternative_exon_end", "strand", "splicemode")] %>% as_tibble,
                                                          tibble_3FT_alternative_exon_only_with_fasta_header) %>% flatten %>% as_tibble
                                                        
                                                        return(final_tibble)
                                                        
                                                      }, .progress = TRUE, .options = future_options(globals = c("rbindlist", "as_tibble", "dplyr::mutate", "dplyr::distinct", "add_column")) )
  
  if (save_workspace_when_done == "DEBUG") {
    save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
  }
  
  # rbind and tibblise
  tibble_summarised_3FT_result_exon_only <- list_summarised_3FT_result_exon_only %>% rbindlist(use.names = TRUE) %>% as_tibble %>% type_convert
  
  # filter for virtual peptides less than 7 AA
  tibble_summarised_3FT_result_exon_only <- tibble_summarised_3FT_result_exon_only[tibble_summarised_3FT_result_exon_only$virtual_peptide_sequence %>% purrr::map(~.x %>% nchar >= 7) %>% unlist %>% which, ]
  
  # write the chucked out items
  tibble_chucked_out_results <- tibble_summarised_3FT_result_exon_only %>% dplyr::filter(translation_frame == "untranslatable")
  
  # filter for untranslatable
  tibble_summarised_3FT_result_exon_only <- tibble_summarised_3FT_result_exon_only %>% dplyr::filter(translation_frame != "untranslatable") %>% dplyr::distinct(virtual_peptide_sequence, fasta_header, .keep_all = TRUE)
  
  # filter substrings
  ## add column to indicate if the virtual peptide per row is a substring of another row
  vector_virtual_peptide_sequence <- tibble_summarised_3FT_result_exon_only$virtual_peptide_sequence
  
  substring_or_not <- future_map2(.x = tibble_summarised_3FT_result_exon_only$virtual_peptide_sequence, 
                                  .y = 1:nrow(tibble_summarised_3FT_result_exon_only),
                                  .f = function(a1, a2) {
                                    
                                    cat(a2, "\n")
                                    
                                    # DEBUG ###
                                    # a1 <- tibble_summarised_3FT_result_exon_only$virtual_peptide_sequence %>% .[[5]]
                                    # a2 <- 1:nrow(tibble_summarised_3FT_result_exon_only) %>% .[[5]]
                                    ###########
                                    
                                    (grepl(x = vector_virtual_peptide_sequence[-a2], pattern = a1) & vector_virtual_peptide_sequence[-a2] != a1) %>% any == TRUE
                                    
                                  }, .progress = TRUE, .options = future_options(globals = c("vector_virtual_peptide_sequence"))) %>% unlist
  
  ## filter
  tibble_summarised_3FT_result_exon_only <- tibble_summarised_3FT_result_exon_only %>% add_column("substring_or_not" = substring_or_not)
  
  # filter out substrings
  tibble_summarised_3FT_result_exon_only_no_substring <- tibble_summarised_3FT_result_exon_only %>% dplyr::filter(substring_or_not == FALSE)
  
  # tally up the number of valid frames we ended up with
  tibble_exons_frame_tally <- tibble_summarised_3FT_result_exon_only_no_substring %>% dplyr::distinct(translation_frame, fasta_header) %>% dplyr::group_by(fasta_header) %>% dplyr::summarise("tally" = n())
  
  cat("\nnumber of exons input: ", list_tibble_master_alternative_exons_chr_start_end_strand_array.tree %>% length, "\n")
  cat("\nnumber of exons translated: ", tibble_summarised_3FT_result_exon_only_no_substring$final_identifier %>% unique %>% length, "\n")
  cat("\naverage number of translation frames for UNIQUE VSRs + exons: ", mean(tibble_exons_frame_tally$tally), "\n")
  
  # WE WRITE A TIBBLE CONTAINING THE GENOME COORD-PEPTIDE MAPPING
  # write a table
  write.table(x = tibble_summarised_3FT_result_exon_only_no_substring, file = paste(output_dir, "/", output_name, "_3FT.summary.info.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # discarded entries
  write.table(x = tibble_chucked_out_results, file = paste(output_dir, "/", output_name, "_discarded.entries.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # FINALLY! WE WRITE THE FASTA!
  write.fasta(sequences = tibble_summarised_3FT_result_exon_only_no_substring$virtual_peptide_sequence %>% array_tree %>% flatten, names = tibble_summarised_3FT_result_exon_only_no_substring$fasta_header, file.out = paste(output_dir, "/", output_name, ".fasta", sep = ""), open = "w", nbchar = 40, as.string = TRUE)
  
  # write final exon table as .bed file
  exon_bed_table <- tibble_summarised_3FT_result_exon_only_no_substring[, c("chr", "alternative_exon_start", "alternative_exon_end", "final_identifier", "strand")] %>% unique %>% setNames(c("chr", "start", "end", "name", "strand")) %>% add_column(., "score" = 1000, .after = "name") %>% type_convert
  
  write.table(exon_bed_table, file = paste(output_dir, "/", output_name, "_exons.bed", sep = ""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = FALSE)
  
}

if (save_workspace_when_done == "YES" | save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# finish counting
tictoc::toc()

q()
