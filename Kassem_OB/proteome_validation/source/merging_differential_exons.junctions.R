# MERGING DIFFERENTIAL JUNCTIONS FROM JUM AND CREATING FASTA HEADERS FOR 3-FRAME TRANSLATION #####
# This script will do the following: 
# 1. Accepts an input file of differential junction results with AS_event_ID and splicemode as columns (either one of: "A3SS_events", "A5SS_events", "cassette_exon_events", "composite_events", "intron_retention", "MXE_events" as per JUM nomenclature. You need to have done this manually beforehand.)
# 2. Matches the AS_event_ID with the junc_coords used by JUM. This requires all the detailed tables for every comparison to be imported.
# 3. Fetches the chr, start, end and strand params of each junc_coord ID using the UNION_junc_coor_with_junction_ID_more_than_xx_read_in_at_least_xx_samples.txt
# 4. Bumps up all co-ordinate values by 1 to account for JUM co-ordinates starting from 0 instead of 1.
# 5. Then matches with a reconstructed (e.g. cufflinks, stringtie, strawberry) AND reference transcriptome annotation (GTF) file (e.g. ensembl) and:
# a) determines which junctions are perfectly flanked by observed transcripts, then,
# b) for those junctions, build a FASTA header based on external_gene_name or transcript_id (ENST)
# c) junctions with no match are discarded 
# 6. Export final table with FASTA header appended to it, ready for 3-frame translation.

# arg1: path to differential junction table
# arg2: path to the actual reconstructed GTF file (e.g. from Cufflinks, Strawberry). NOT THE CONTAINING DIRECTORY. tip: for better junction matching, combine the reconstructed GTF with reference GTF beforehand e.g. using StringTie
# arg3: path to the actual genome FASTA file. Ideally from Ensembl... you need the whole primary assembly. Not separate files by chromosome.
# arg4: how many nucleotides to translate upstream W.R.T. the transcript. for junction-based mode, the total nt. length to be translated will be arg4 + arg5
# arg5: how many nt to translate downstream W.R.T. the transcript. for exon-based mode, this will be 0 by default. the total arg4 + arg5 will indicate how many nt to translate beyond the exon.
# arg6: output directory
# arg7: MODE: currently only accepts "exon" or "junction".

print(commandArgs())

args = commandArgs(trailingOnly=TRUE)

print(args)

print(paste("number of arguments: ", length(args), sep = ""))

# test if there are 7 arguments: if not, return an error
if (length(args) != 7) {
  stop("check the amount of arguments you have. should have specified 7", call.=FALSE)
}


# SET ENVIRONMENT ##########

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("tidyverse", "purrr", "data.table", "dplyr", "rtracklayer", "biomaRt"))

library(tidyverse)
library(purrr)
library(data.table)
library(dplyr)
library(rtracklayer)
library(biomaRt)

## Directories

differential_junction_table_path <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/input_files/5x5_SOM_1649_isoforms_PSI_OB_diff_any_qvalue0.01_any_deltaPSI_greaterthan_0.2.txt"

JUM_results_dir <- "/home/z3463471/PGNEXUS_kassem_MSC/OB/analysis_JUM/run_1_PGNEXUS_OBseries_pvalue1/results/"

UNION_junc_coor_path <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/bam_bai_sj_sam_files/JUM_diff/UNION_junc_coor_with_junction_ID_more_than_5_read_in_at_least_3_samples.txt"

ref_gtf_path <- "/srv/scratch/z3463471/hg38_ensembl_reference/gtf/Homo_sapiens.GRCh38.98.gtf"
  
reconstructed_gtf_path <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/alltimepoints_reconstructed_stringtiemerged.gtf"

output_file_name <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/results_database_generation/UNION_junc_coor_differential_allSOM.txt"

# Pooling AS junction co-ordinates from JUM

# accepts:
# 1. table of differential junctions (result from JUM), must contain AS_event_ID and splicemode columns
# 2. all the detailed tables used in analysis, as is (untouched)
# 3. UNION_junc_coor file from JUM, which contains all the junction ID to genome co-ordinate mapping.

## Import the detailed tables into the environment

list_of_timepoint_comparisons_final <- read.delim(paste(JUM_results_dir, "list_of_timepoint_comparisons_final.txt", sep = ""), stringsAsFactors = FALSE, sep = "\t", header = FALSE, row.names = NULL) %>% array_tree %>% flatten

list_of_AS_events <- c("A3SS_events", "A5SS_events", "cassette_exon_events", "composite_events", "intron_retention", "MXE_events") %>% array_tree %>% flatten

# read tables into environment

list_of_detailed_tables <- purrr::map(.x = list_of_timepoint_comparisons_final, .f = ~purrr::map2(.x = .x, .y = list_of_AS_events, .f = ~read_delim(file = paste(JUM_results_dir, "final_JUM_output_", .x, "/", list.files(path = paste(JUM_results_dir, "final_JUM_output_", .x, "/", sep = ""), pattern = paste("(.)", .y, "(.*)detailed.txt", sep = "")), sep = ""), delim = "\t") %>% add_column(., "splicemode" = .y)) %>% set_names(list_of_AS_events)) %>% set_names(list_of_timepoint_comparisons_final)

# bind rows for each comparison

list_of_detailed_tables_2 <- list_of_detailed_tables %>% purrr::map(.x = ., .f = ~rbindlist(.x))

# rename the colnames to be not comparison-specific - we will add a comparison column later on

list_of_detailed_tables_3 <- purrr::map2(.x = list_of_detailed_tables_2, .y = list_of_timepoint_comparisons_final, .f = ~setNames(., c("Gene", "AS_event_ID", "AS_structure_ID", "sub_junction_ID", "sub-junction_dispersion_estimate", "LRT_statistic-full_vs_reduced", "LRT_p_value-full_vs_reduced", "BH_adjusted_p-values", "logCPM_1", "logCPM_2", "fitting_parameter_log2fold_change_2_from_1", "sub_junction_chr", "sub_junction_start_coor", "sub_junction_end_coor", "sub_junction_size", "sub_junction_strand", "raw_count.1", "raw_count.2", "raw_count.3", "raw_count.4", "raw_count.5", "raw_count.6", "percentage_usage.1", "percentage_usage.2", "percentage_usage.3", "percentage_usage.4", "percentage_usage.5", "percentage_usage.6", "deltaPSI", "splicemode")) %>% add_column(., "comparison" = .y)) 

wide_table_of_all_detailed_tables <- list_of_detailed_tables_3 %>% rbindlist %>% as_tibble

### Import table of differential junctions

differential_junction_table <- read_delim(paste(differential_junction_table_path), delim = "\t")

### Filter the detailed tables and extract junction IDs

wide_table_of_all_detailed_tables_diffonly <- dplyr::semi_join(wide_table_of_all_detailed_tables, differential_junction_table, by = c("AS_event_ID", "splicemode"))

AS_junction.IDs <- tibble("junction_ID" = wide_table_of_all_detailed_tables_diffonly$AS_structure_ID %>% gsub(x = ., pattern = "(5|3)_(.*)", replacement = "\\2"))

AS_junction_IDs_split <- AS_junction.IDs$junction_ID %>% gsub(x = ., pattern = "Junction_", replacement = "") %>% strsplit(., split = "_")

AS_junction_IDs_split_2 <- tibble("junction_ID" = paste("Junction_", AS_junction_IDs_split %>% unlist, sep = "")) %>% unique

print(paste("there are", nrow(AS_junction_IDs_split_2), "unique junctions which have been plausibly differentially spliced"))


### import the UNION junc coordinate table and filter it for differential junctions we have just found

UNION_junc_coor_table <- read_delim(paste(UNION_junc_coor_path), delim = "\t", col_names = c("chr", "start", "end", "strand", "junction_ID"))

UNION_junc_coor_table_2 <- UNION_junc_coor_table

UNION_junc_coor_table_2[, "strand"] <- gsub(x = UNION_junc_coor_table_2[, "strand"] %>% as.data.frame %>% unlist, pattern = "0", replacement = ".")

UNION_junc_coor_table_3 <- dplyr::semi_join(UNION_junc_coor_table_2, AS_junction_IDs_split_2, by = "junction_ID") %>% na.omit

# shift back the coord positions because JUM co-ordinates start from 0

UNION_junc_coor_table_4 <- UNION_junc_coor_table_3
UNION_junc_coor_table_4[, "start"] <- UNION_junc_coor_table_4[, "start"] + 1
UNION_junc_coor_table_4[, "end"] <- UNION_junc_coor_table_4[, "end"] + 1

### import reference GTF, filter it for transcripts with JUNCTION-FLANKING EXONS ONLY and append to UNION junc coord table

# this annotation will be placed as FASTA headers in the custom database.
# 
# the plan is to grab reference annotation where possible. If not possible, refer to the reconstructed transcriptome merged from all the timepoints and extract the transcript id.
# 
# case 1: flanked by reference transcript: header contains ONLY REF transcript and REF gene name regardless of whether the junction is found in the reconstructed GTF or not
# case 2: junction not in ref GTF BUT IS in an annotated gene AND junction IS in reconstructed GTF: header contains REF external gene names, comma separated (no space), transcript id comes from reconstructed GTF, separated by comma.
# case 3: junction is not in ref GTF and DOES NOT lie in any annotated gene region: header contains no external gene name (NA) and transcript id(s) come(s) from reconstructed GTF.
# case 4: junction is not in either ref nor reconstructed GTF: junction is discarded because there's no way to translate it.

# FUNCTION TO EXTRACT TRANSCRIPTS WITH JUNCTION-FLANKING EXONS.
# NOTE: to be used with purrr
# input: spliceregion_list: a list containing details of ONE junction: $chr, $start, $end, $strand and $junction_ID
# tibble_gtf_table: rtracklayer::import(...) %>% as_tibble. can be any GTF. reconstructed or reference.
# index: loop progress marker to be used with imap

extract_junction.flanking.exons <- function(spliceregion_list, tibble_gtf_table, index) {
  
  # DEBUG ###################

  # spliceregion_list <- UNION_junc_coor_table_4_array.tree[[1]]
  # gtf_table <- tibble_ref_gtf

  ###########################

  print(paste("now processing junction number", index))
  
if (spliceregion_list$strand == ".") {
    
    tibble_gtf_subset_flanking_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == spliceregion_list$chr, ] %>% .[.$start <= ((spliceregion_list$end %>% as.numeric) + 1) & .$end >= ((spliceregion_list$start %>% as.numeric) - 1), ] %>% .[!(.$start <= (spliceregion_list$end %>% as.numeric) & .$end >= (spliceregion_list$start %>% as.numeric)), ] %>% .[.$type == "exon", ]
    
  } else if (spliceregion_list$strand == "+" | spliceregion_list$strand == "-") {
    
    tibble_gtf_subset_flanking_exons <- tibble_gtf_table[tibble_gtf_table$seqnames == spliceregion_list$chr, ] %>% .[.$strand == spliceregion_list$strand, ] %>% .[.$start <= ((spliceregion_list$end %>% as.numeric) + 1) & .$end >= ((spliceregion_list$start %>% as.numeric) - 1), ] %>% .[!(.$start <= (spliceregion_list$end %>% as.numeric) & .$end >= (spliceregion_list$start %>% as.numeric)), ] %>% .[.$type == "exon", ]
    
  } else {
    
    stop("Could not match the strand information in the transposed differential-only UNION_junc_coor_table. Make sure that the \"strand\" column in the UNION_junc_coor_table contains only +, - or .")
    
  }
  
  list_of_junction_associated_transcripts <- tibble_gtf_subset_flanking_exons$transcript_id %>% unique %>% array_tree %>% flatten
  
  # make a list for each transcript that directly flanks a junction.
  # then filter so that there are only a) exon PAIRS which b) are directly connected in the mature (spliced) transcript
  
  list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_gtf_subset_flanking_exons[tibble_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts) %>% keep(.x = ., .p = ~nrow(.x) == 2) %>% keep(.x = ., .p = ~abs((.x[2, "exon_number"] %>% paste %>% as.numeric) - (.x[1, "exon_number"] %>% paste %>% as.numeric)) == 1)
  
  return(list_of_tibbles_flanking_exon_gtf.entries_per_transcript)
  
}

# END extract_junction.flanking.exons #######################

# function to extract gene name from reference GTF 
# only FULL overlap considered, where the query junction interval is a subset of any transcript

extract_external.gene.name_from_reference_GTF <- function(spliceregion_list, ref_gtf_table, index) {
  
  # DEBUG ###################

  # spliceregion_list <- UNION_junc_coor_table_4_array.tree[[1]]
  # ref_gtf_table <- tibble_ref_gtf
  # index <- 1

  ###########################

  print(paste("now processing junction number", index))
  
if (spliceregion_list$strand == ".") {
    
    ref_gtf_subset <- ref_gtf_table[ref_gtf_table$seqnames == spliceregion_list$chr, ] %>% .[.$start <= (spliceregion_list$start %>% as.numeric) & .$end >= (spliceregion_list$end %>% as.numeric), ]
    
  } else if (spliceregion_list$strand == "+" | spliceregion_list$strand == "-") {
    
    ref_gtf_subset <- ref_gtf_table[ref_gtf_table$seqnames == spliceregion_list$chr, ] %>% .[.$start <= (spliceregion_list$start %>% as.numeric) & .$end >= (spliceregion_list$end %>% as.numeric), ]
    
  } else {
    
    stop("Could not match the strand information in the transposed differential-only UNION_junc_coor_table. Make sure that the \"strand\" column in the UNION_junc_coor_table contains only +, - or .")
    
  }
  
  column_matched_genes <- ref_gtf_subset$gene_name %>% unique
  
  return(column_matched_genes)
  
}

# END ########################

# FUNCTION TO RETURN FASTA HEADER FROM MATCHED REF/RECONSTRUCTED GTF ROWS.
# NOTE: to be used with purrr
# input: a list containing ref/reconstructed GTF matches for ONE junction, totalling 4 elements: 1. subsets of reference GTF (by transcript) each containing EXACTLY 2 entries of the reference GTF tibble - these describe the exons flanking the junction. 2. subsets of reconstructed GTF giving similar information. 3. amount of transcripts matched to ref. GTF 4. amount of transcripts matched to reconstructed GTF.

# SUMMARY OF BEHAVIOUR:
# the plan is to grab reference annotation where possible. If not possible, refer to the reconstructed transcriptome merged from all the timepoints and extract the transcript id.
# 
# case 1: flanked by reference transcript: header contains ONLY REF transcript and REF gene name regardless of whether the junction is found in the reconstructed GTF or not
# case 2: junction not in ref GTF BUT IS in an annotated gene AND junction IS in reconstructed GTF: header contains REF external gene names, comma separated (no space), transcript id comes from reconstructed GTF, separated by comma.
# case 3: junction is not in ref GTF and DOES NOT lie in any annotated gene region: header gene name is called "intergenic" and transcript id(s) come(s) from reconstructed GTF.
# case 4: junction is not in either ref nor reconstructed GTF: junction is discarded because there's no way to translate it. header is NA, and we shall execute na.omit to ultimately remove it from the junction co-ordinate table
# 
# FASTA HEADER FORMAT:
# format: >diff.splice.tool_transcript.reconstruction.tool|junction_ID|ENST.id|external_gene_name
# example: >JUM_strawberry|Junction_100084|ENST00000536720.1|ALKBH2

generate_fasta_header_from_subset_gtf <- function(element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length, diff.splice.tool_name, transcript.reconstruction.tool_name) {
  
  # DEBUG ###################
  
  # diff.splice.tool_name <- "JUM"
  # transcript.reconstruction.tool_name <- "strawberry"
  
  ###########################
  
  ref_gtf_combined <- element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$reference_gtf_flanking_exon_entries_per_transcript %>% rbindlist %>% as_tibble
  reconstructed_gtf_combined <- element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$reconstructed_gtf_flanking_exon_entries_per_transcript %>% rbindlist %>% as_tibble
  
  # case 1
  if (element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$reference.gtf_number.of.transcripts.matched %>% paste %>% as.numeric != 0) {
    
    gene_ids <- ref_gtf_combined$gene_name %>% unique
    
    transcript_ids <- ref_gtf_combined$transcript_id %>% unique
    # case 2
  } else if (element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$reference.gtf_number.of.transcripts.matched %>% paste %>% as.numeric == 0 & element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$reconstructed.gtf_number.of.transcripts.matched %>% paste %>% as.numeric != 0 &
             length(element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$external_gene_names) != 0) {
    
    gene_ids <- element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$external_gene_names
    transcript_ids <- reconstructed_gtf_combined$transcript_id %>% unique
    # case 3
  } else if (element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$reference.gtf_number.of.transcripts.matched %>% paste %>% as.numeric == 0 & element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$reconstructed.gtf_number.of.transcripts.matched %>% paste %>% as.numeric != 0 &
             length(element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$external_gene_names) == 0) {
    
    gene_ids <- "intergenic"
    transcript_ids <- reconstructed_gtf_combined$transcript_id %>% unique
    
  }
  
  # case 4
  
  if (element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$reference.gtf_number.of.transcripts.matched %>% paste %>% as.numeric == 0 & element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$reconstructed.gtf_number.of.transcripts.matched %>% paste %>% as.numeric == 0) {
    
    fasta_header <- NA
    
  } else {
    
    # format: >diff.splice.tool_transcript.reconstruction.tool|junction_ID|ENST.id|external_gene_name
    # example: >JUM_strawberry|Junction_100084|ENST00000536720.1|ALKBH2
    
    fasta_header <- paste(">", 
                          diff.splice.tool_name, 
                          "_", 
                          transcript.reconstruction.tool_name, 
                          "|", 
                          element_of_list_ref.and.reconstructed_gtf_matching.entries_with.length$junction_ID, 
                          "|", 
                          paste(transcript_ids, collapse = ","), 
                          "|",
                          paste(gene_ids, collapse = ","), sep = "")
    
  }
  
  return(fasta_header)
  
}

# END generate_fasta_header_from_subset_gtf ####################################

# DEBUG ###############

# set.seed(7)
# UNION_junc_coor_table_4_array.tree <- UNION_junc_coor_table_4 %>% .[sample(1:100, 20), ] %>% array_tree

#######################

UNION_junc_coor_table_4_array.tree <- UNION_junc_coor_table_4 %>% array_tree

tibble_ref_gtf <- rtracklayer::import(ref_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)

tibble_reconstructed_gtf <- rtracklayer::import(reconstructed_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)

# FUNCTION TO ANNOTATE EACH JUNCTION WITH FASTA HEADERS from reference GTF

# format: >diff.splice.tool_transcript.reconstruction.tool|junction_ID|ENST.id|external_gene_name
# example: >JUM_strawberry|Junction_100084|ENST00000536720.1|ALKBH2

list_ref.and.reconstructed_gtf_matching.entries <- purrr::imap(.x = UNION_junc_coor_table_4_array.tree, .f = ~list("reference_gtf_flanking_exon_entries_per_transcript" = extract_junction.flanking.exons(.x, tibble_ref_gtf, .y), "reconstructed_gtf_flanking_exon_entries_per_transcript" = extract_junction.flanking.exons(.x, tibble_reconstructed_gtf, .y), "external_gene_names" = extract_external.gene.name_from_reference_GTF(.x, tibble_ref_gtf, .y), "junction_ID" = .x$junction_ID))

list_ref.and.reconstructed_gtf_matching.entries_with.length <- purrr::map(.x = list_ref.and.reconstructed_gtf_matching.entries, .f = ~purrr::splice(.x, "reference.gtf_number.of.transcripts.matched" = length(.x$reference_gtf_flanking_exon_entries_per_transcript), "reconstructed.gtf_number.of.transcripts.matched" = length(.x$reconstructed_gtf_flanking_exon_entries_per_transcript)))

list_fasta_headers <- purrr::map(.x = list_ref.and.reconstructed_gtf_matching.entries_with.length, .f = ~generate_fasta_header_from_subset_gtf(.x, "JUM", "strawberry"))

UNION_junc_coor_table_5 <- add_column(UNION_junc_coor_table_4, "fasta_header" = list_fasta_headers %>% unlist)

UNION_junc_coor_table_5 <- UNION_junc_coor_table_5 %>% na.omit

### export final UNION_junc_coor table

write.table(UNION_junc_coor_table_5, file = paste(output_file_name), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


