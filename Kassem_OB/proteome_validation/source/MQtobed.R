script_description <- "# GENERATE COVERAGE .BED FILE FROM MAXQUANT OUTPUT TABLE ######
This script will automatically generate a .bed file showing the positions in the genome which codes for identified peptides, phosphosites etc... 

BEHAVIOUR:
1. Imports the Maxquant .txt file
2. Extracts all the protein-relative start & end, along with their parent protein IDs e.g. uniprotkb_entry
3. Converts to CDS-relative start & end by multiplying by 3.
4. Converts to transcript-relative coords by accounting for start codon position. NOTE: to do this, first import a reference GTF as well as a table of parent protein ID to ENSP matching. if there is no start_codon, then the entry will be discarded.
5. Converts to genome-relative coords.
6. Writes the blocked bedfile suitable for use in IGV.
"

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

library(tidyverse)
library(furrr)
options(future.globals.maxSize = 30000000000, future.fork.enable = TRUE)
library(rtracklayer)
library(data.table)
library(optparse)

library(tictoc)
# start counting execution time of the whole script
tictoc::tic("Overall execution time")

# manage arguments
# manage arguments
list_input_arg_info <- list(
  "1" = make_option(c("-M", "--maxquant_table_path"), type = "character", default = NULL, 
                    help = "Compulsory. Path to the actual table in the txt/ folder outputted by maxquant.", metavar = "character"),
  "2" = make_option(c("-T", "--maxquant_table_type"), type = "character", default = "NULL", 
                    help = "Optional. You can tell the script what kind of table it is, i.e. peptides, phosphosites etc...", metavar = "character"),
  "3" = make_option(c("-S", "--score_FDR_localisation_cutoffs"), type = "character", default = NULL, 
                    help = "Compulsory. Specify <score;FDR(or PEP);localisation_probability>. IN THAT ORDER. Separated by a semicolon.", metavar = "character"),
  "4" = make_option(c("-P", "--protein_id_mapping_path"), type = "character", default = NULL, 
                    help = "Compulsory. Path to a TAB-SEPARATED table which maps maxquant protein name to GTF protein ID", metavar = "character"),
  "5" = make_option(c("-Q", "--protein_id_mapping_column_name_info"), type = "character", default = NULL, 
                    help = "Compulsory. Gotta specify 1. which column contains your Maxquant protein name and 2. which column contains your GTF protein_id. IN THAT ORDER. Separated by a semicolon, \";\". E.g. \"uniprotkb_entry;ensembl_peptide_id\"", metavar = "character"),
  "6" = make_option(c("-G", "--reference_gtf_path"), type = "character", default = NULL, 
                    help = "Compulsory. Path to the actual reference GTF file (e.g. from Cufflinks, Strawberry) that you want to use to locate the position where proteins are coded. NOT THE CONTAINING DIRECTORY.", metavar = "character"),
  "7" = make_option(c("-H", "--reference_gtf_field_name_info"), type = "character", default = NULL, 
                    help = "Compulsory. Gotta specify which field in the reference GTF contains protein_ids such as ENSP identifiers. These would have been in your protein ID mapping file. If not specified, the program will try to guess which column has ENSP identifiers.", metavar = "character"),
  
  "8" = make_option(c("-R", "--track_colour"), type = "character", default = "255,0,0", 
                    help = "Optional. Specify the colour you want the feature tracks to be. Default is \"0,0,255\", which is blue. See the RGB colour specification for more options. Use a an online colour picker web tool to get the RGB code of the colour you like.", metavar = "character"),
  "9" = make_option(c("-D", "--output_dir"), type = "character", default = NULL, 
                    help = "Compulsory. output file directory. where do you want to save the annotated exon table? IMPORTANT: MUST BE A FULL DIRECTORY AND NOT A FILE PATH. e.g. correct: ~/outputdir/ correct: ~/outputdir incorrect: /outputdir/bedfile.bed", metavar = "character"),
  "10" = make_option(c("-O", "--output_name"), type = "character", default = NULL, 
                    help = "Compulsory. output file name, to be saved in the output directory a.k.a. what do you want to save the annotated exon table as? IMPORTANT: MUST BE A STRING WITHOUT THE EXTENSION AND NOT A DIRECTORY. THE .bed EXTENSION WILL AUTOMATICALLY BE ADDED FOR THE OUTPUT FILE. e.g. correct: \"bedfile\" incorrect: \"bedfile.txt\" incorrect: bedfile/", metavar = "character"),
  
  "11" = make_option(c("-C", "--ncores"), type = "character", default = 0, 
                    help = "Optional. Number of cores to use. possible inputs: numbers 1 to any integer. By default, uses all cores (ncores = 0). If a single number is specified, it will just tell future to loop thru chromosomes in parallel using the specified core count. If numberxnumber for example 7x4 then 28 cores will be used. 7 for chromosomes and 4 for inside each chromosome.", metavar = "character"), 
  "12" = make_option(c("-V", "--save_workspace_when_done"), type = "character", default = "NO",
                     help = "Turn this on if you want to save the R workspace in the same name as the --output_name. YES: saves at the end. DEBUG: saves at each critical step. NO: doesn't save.", metavar = "character")
)
input_arg_info <- OptionParser(option_list = list_input_arg_info, description = script_description)
input_args <- input_arg_info %>% parse_args

# check if the input arguments are O.K
# if ((list(input_args$reconstructed_transcript_gtf_path, input_args$reference_genome_fasta_dir, input_args$output_name) %>% lapply(is.null) %>% unlist %>% any == TRUE) | 
#     (input_args$chrmode == 2 & is.null(input_args$nonchrname) == TRUE)) {
#   
#   print_help(input_arg_info)
#   
#   stop("Make sure you entered the arguments correctly", call. = FALSE)
#   
# }

# DEBUG #######

# maxquant_table_path <- "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/phosphoproteomic_analysis/analysis_maxquant/results/2020_phosphoproteome_OBseries_con_sp.hsa.canonical.isoforms/txt/Phospho (STY)Sites.txt"
# maxquant_table_path <- "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/phosphoproteomic_analysis/analysis_maxquant/results/2020_phosphoproteome_OBseries_con_sp.hsa.canonical.isoforms/txt/peptides.txt"
# 
# maxquant_table_type <- "NULL"
# score_FDR_localisation_cutoffs <- "30;0.01;0.75"
# protein_id_mapping_path <- "/mnt/Tertiary/sharedfolder/table_ENSP_to_uniprot_entry_mapping.txt"
# protein_id_mapping_column_name_info <- "uniprotkb_entry;ensembl_peptide_id"
# reference_gtf_path <- "/mnt/Tertiary/sharedfolder/hg38_ensembl_reference/gtf/Homo_sapiens.GRCh38.98.gtf"
# reference_gtf_field_name_info <- "protein_id"
# track_colour <- "255,0,0"
# output_dir <- "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/phosphoproteomic_analysis/analysis_maxquant/results/2020_phosphoproteome_OBseries_con_sp.hsa.canonical.isoforms/"
# output_name <- "debug_swissprot_bedfile_generator_peptides"
# ncores <- "16x2"
# save_workspace_when_done <- "NO"

###############

maxquant_table_path <- input_args$maxquant_table_path
maxquant_table_type <- input_args$maxquant_table_type
score_FDR_localisation_cutoffs <- input_args$score_FDR_localisation_cutoffs
protein_id_mapping_path <- input_args$protein_id_mapping_path
protein_id_mapping_column_name_info <- input_args$protein_id_mapping_column_name_info
reference_gtf_path <- input_args$reference_gtf_path
reference_gtf_field_name_info <- input_args$reference_gtf_field_name_info
track_colour <- input_args$track_colour
output_dir <- input_args$output_dir
output_name <- input_args$output_name
ncores <- input_args$ncores
save_workspace_when_done <- input_args$save_workspace_when_done

cat("maxquant_table_path:", maxquant_table_path, "\n")
cat("maxquant_table_type:", maxquant_table_type, "\n")
cat("score_FDR_localisation_cutoffs:", score_FDR_localisation_cutoffs, "\n")
cat("protein_id_mapping_path:", protein_id_mapping_path, "\n")
cat("protein_id_mapping_column_name_info:", protein_id_mapping_column_name_info, "\n")
colname_mapping_info_table_containing_uniprot_id <- protein_id_mapping_column_name_info %>% strsplit(split = ";") %>% unlist %>% .[1]
colname_mapping_info_table_containing_ensembl_peptide_id <- protein_id_mapping_column_name_info %>% strsplit(split = ";") %>% unlist %>% .[2]
cat("this means:\n")
cat("the column name of the mapping table with maxquant protein IDs:", colname_mapping_info_table_containing_uniprot_id, "\n")
cat("the column name of the mapping table with GTF protein IDs:", colname_mapping_info_table_containing_ensembl_peptide_id, "\n")
cat("reference_gtf_path:", reference_gtf_path, "\n")
cat("reference_gtf_field_name_info:", reference_gtf_field_name_info, "\n")
cat("track_colour:", track_colour, "\n")
cat("output_dir:", output_dir, "\n")
cat("output_name:", output_name, "\n")
cat("ncores:", ncores, "\n")
cat("save_workspace_when_done:", save_workspace_when_done, "\n")

if(!dir.exists(output_dir) ) {
  dir.create(output_dir, recursive = TRUE)}

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

## function to split a whole table by a delimiter in a column
split_delimited_columns_in_table <- function(input_table, target_colname, split, columns_to_deduplicate = NULL) {
  
  # DEBUG ###
  # input_table <- test
  # target_colname <- c("Proteins", "Positions.within.proteins", "Fasta.headers", "Leading.proteins")
  # split = "\\;"
  # columns_to_deduplicate <- c("id")
  ###########
  
  # list-ify the target column
  list_target_column_strsplit_per_element <- input_table[, target_colname] %>% array_tree(margin = 2) %>% purrr::map(~.x %>% unlist %>% strsplit(., split = split))
  
  # map length
  vector_sum_lengths <- list_target_column_strsplit_per_element %>% purrr::map(~length(.x %>% unlist)) %>% unlist(use.names = FALSE) %>% unique
  
  # check. if the split lengths are different, then die.
  if (length(vector_sum_lengths) != 1) {
    
    stop("split lengths are uneven across specified columns.")
    
  }
  
  # repeat table according to the split lengths
  vector_split_lengths <- purrr::map_depth(.x = list_target_column_strsplit_per_element, .depth = 2, .f = ~length(.x)) %>% 
    purrr::map(~unlist(.x)) %>%
    purrr::pmap(.f = ~max(...)) %>%
    unlist
  
  row_indices_of_table_repeated_by_split <- purrr::map2(.x = 1:nrow(input_table), .y = vector_split_lengths, .f = ~rep(x = .x, times = .y)) %>% unlist
  
  input_table_repeated_by_split <- input_table[row_indices_of_table_repeated_by_split, ]
  
  # replace target column with split values
  input_table_repeated_by_split[, target_colname] <- list_target_column_strsplit_per_element %>% purrr::map(~unlist(.x)) %>% purrr::reduce(cbind)
  
  split_table <- input_table_repeated_by_split
  
  # if specified, append an index to a particular column
  if (is.null(columns_to_deduplicate) == FALSE) {
    
    # get the duplicated row indices where split lengths > 1
    indices_of_duplicates <- which(vector_split_lengths > 1)
    
    # get the repetition number where split lengths > 1
    repetition_numbers_of_duplicates <- vector_split_lengths[which(vector_split_lengths > 1)]
    
    # list-ify the columns to be appended
    list_deduplicated_columns <- input_table[, columns_to_deduplicate] %>% array_tree(margin = 2) %>% purrr::map(~array_tree(.x))
    
    # map over each column, split the target element and add _[0-9]+
    list_deduplicated_columns_split <- purrr::map(.x = list_deduplicated_columns, .f = function(a1) {
      
      # map a subset each of the L2 (elements of a column)
      a1[indices_of_duplicates] <- purrr::map2(.x = a1[indices_of_duplicates], .y = repetition_numbers_of_duplicates, 
                                               .f = ~rep(.x, times = .y) %>% unlist %>% paste(., 1:.y, sep = "_"))
      
      return(a1 %>% unlist)
      
    } )
    
    # tibblise
    tibble_deduplicated_columns_split <- list_deduplicated_columns_split %>% as_tibble
    
    # add back every row onto the split table
    for (dedupe_colname in columns_to_deduplicate) {
      
      split_table[, dedupe_colname] <- tibble_deduplicated_columns_split[, dedupe_colname]
      
    }
    
  }
  
  return(split_table)
  
}
# END split_delimited_columns_in_table() ###

# BEGIN EXECUTION #################################

cat("importing maxquant table\n")
tibble_maxquant_table <- read.delim(file = maxquant_table_path, stringsAsFactors = FALSE, sep = "\t", header = TRUE, row.names = NULL, check.names = FALSE) %>% as_tibble

cat("checking cutoffs\n")
score_cutoff <- score_FDR_localisation_cutoffs %>% strsplit(split = ";") %>% unlist %>% .[1]
FDR_cutoff <- score_FDR_localisation_cutoffs %>% strsplit(split = ";") %>% unlist %>% .[2]
localisation_probability_cutoff <- score_FDR_localisation_cutoffs %>% strsplit(split = ";") %>% unlist %>% .[3]

cat("checking the peptide ID mapping table\n")
tibble_uniprot_to_ENSP_mapping_table <- read.delim(file = protein_id_mapping_path, stringsAsFactors = FALSE, sep = "\t", header = TRUE, row.names = NULL, check.names = FALSE) %>% as_tibble

if (any(is.null(c(colname_mapping_info_table_containing_uniprot_id, colname_mapping_info_table_containing_ensembl_peptide_id)))) {
  stop("Missing field. Check the specified --protein_id_mapping_column_name_info option. User specified: ", protein_id_mapping_column_name_info, "\n")
}

if (any(c(colname_mapping_info_table_containing_uniprot_id, colname_mapping_info_table_containing_ensembl_peptide_id) %in% colnames(tibble_uniprot_to_ENSP_mapping_table) == FALSE)) {
  stop("Specified mapping column names are not in the protein ID mapping table. Check the specified --protein_id_mapping_column_name_info option. User specified: ", protein_id_mapping_column_name_info, "\nThe column names of the mapping table: ", colnames(tibble_uniprot_to_ENSP_mapping_table) %>% paste(collapse = "; "), 
       "\nThe problem column name specifications: ", c(colname_mapping_info_table_containing_uniprot_id, colname_mapping_info_table_containing_ensembl_peptide_id)[c(colname_mapping_info_table_containing_uniprot_id, colname_mapping_info_table_containing_ensembl_peptide_id) %in% colnames(tibble_uniprot_to_ENSP_mapping_table) == FALSE] %>% paste(collapse = "; "), 
       "\n")
}

cat("checking maxquant table structure and filter the table\n")
# check if the expected columns are present given the table type
## retrieve file name from full path
maxquant_file_name <- gsub(x = maxquant_table_path, pattern = "(.*)\\/(.*)$", replacement = "\\2")

# phosphosites
if ((grepl(x = maxquant_file_name, pattern = "Phospho.*Sites.txt", ignore.case = TRUE) & maxquant_table_type == "NULL") | 
    maxquant_table_type == "phosphosites") {
  
  detected_maxquant_table_type <- "phosphosites"
  
  tibble_maxquant_table_filtered_and_split <- tibble_maxquant_table %>%
    # filter for class I phosphosites
    # remove the reverse hits or table split wont work
    dplyr::filter(`Localization prob` >= 0.75 & `PEP` <= 0.01 & Score >= 30 & `Reverse` != "+") %>%
    # dedupe in preparation for the next step
    split_delimited_columns_in_table(input_table = ., target_colname = c("Proteins", "Positions within proteins"), split = ";", columns_to_deduplicate = "id") %>%
    dplyr::select(Proteins, `Positions within proteins`, contains("probabilities",ignore.case = TRUE), Score) %>%
    dplyr::rename(!!colname_mapping_info_table_containing_uniprot_id := "Proteins",
                  "start_position_within_protein" = "Positions within proteins",
                  "second_id" = colnames(.)[grep(x = colnames(.), pattern = "probabilities",ignore.case = TRUE)]) %>% 
    dplyr::mutate("end_position_within_protein" = `start_position_within_protein`)
  
  if (any(c("Proteins", "Positions within proteins") %in% colnames(tibble_maxquant_table) == FALSE)) {
    stop("We couldn't find the correct protein columns of the maxquant table. This could mean that the maxquant table name is incorrect/misleading, or the user specified option --maxquant_table_type is incorrect. It could also be that the version of MaxQuant is too new. \nThe column that could not be found is/are: ", c("Proteins", "Positions within proteins")[c("Proteins", "Positions within proteins") %in% colnames(tibble_maxquant_table) == FALSE] %>% paste(collapse = "; "), "\n")
  }
  
  # peptides
} else if ((grepl(x = maxquant_file_name, pattern = "peptides.txt", ignore.case = TRUE) & maxquant_table_type == "NULL") | 
           maxquant_table_type == "peptides") {
  
  detected_maxquant_table_type <- "peptides"
  
  if (any(c("Leading razor protein", "Start position", "End position") %in% colnames(tibble_maxquant_table) == FALSE)) {
    stop("We couldn't find the correct protein columns of the maxquant table. This could mean that the maxquant table name is incorrect/misleading, or the user specified option --maxquant_table_type is incorrect. It could also be that the version of MaxQuant is too new. 
         \nThe column that could not be found is/are: ", c("Leading razor protein", "Start position", "End position")[c("Leading razor protein", "Start position", "End position") %in% colnames(tibble_maxquant_table) == FALSE] %>% paste(collapse = "; "), 
         "\n")
  }
  
  tibble_maxquant_table_filtered_and_split <- tibble_maxquant_table %>%
    # filter for class I phosphosites
    # remove the reverse hits or table split wont work
    dplyr::filter(`PEP` <= 0.01 & Score >= 30 & `Reverse` != "+") %>%
    # dedupe in preparation for the next step
    split_delimited_columns_in_table(input_table = ., target_colname = c("Leading razor protein", "Start position", "End position"), split = ";", columns_to_deduplicate = "id") %>%
    dplyr::select(`Leading razor protein`, `Start position`, `End position`, Sequence, Score) %>%
    dplyr::rename(!!colname_mapping_info_table_containing_uniprot_id := "Leading razor protein",
                  "start_position_within_protein" = "Start position",
                  "end_position_within_protein" = "End position",
                  "second_id" = "Sequence")
  
}

cat("detected MaxQuant table type: ", detected_maxquant_table_type, "\n")

if (detected_maxquant_table_type == "phosphosites" & any(is.na(c(score_cutoff, FDR_cutoff, localisation_probability_cutoff)))) {
  stop("Check the specified --score_FDR_localisation_cutoffs option. User specified: ", score_FDR_localisation_cutoffs)
} else if (detected_maxquant_table_type == "peptides" & any(is.na(c(score_cutoff, FDR_cutoff)))) {
  stop("Check the specified --score_FDR_localisation_cutoffs option. User specified: ", score_FDR_localisation_cutoffs)
}

cat("import reference transcriptome GTF\n")
# extract only protein_coding transcripts
# if user has specified to output the FASTA, then we consider all transcripts. 
# if on poison exon detection mode only, then we go for only protein_coding transcript_biotype.
tibble_ref_gtf <- rtracklayer::import(reference_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)

cat("checking reference transcriptome GTF\n")
# automatically detect if exons are always numbered in increasing order regardless of strand (common for ref. transcripts)
## sample the first transcript on the negative strand with more than 1 exon
temp_number <- 1

first_transcript_id <- tibble_ref_gtf$transcript_id %>% unique %>% na.omit %>% .[temp_number]

while (tibble_ref_gtf[tibble_ref_gtf$transcript_id == first_transcript_id, "exon_number"] %>% nrow == 1 | 
       tibble_ref_gtf[tibble_ref_gtf$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "+" | 
       tibble_ref_gtf[tibble_ref_gtf$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "." | 
       tibble_ref_gtf[tibble_ref_gtf$transcript_id == first_transcript_id, "strand"] %>% unlist %>% unique %>% na.omit %>% paste == "*") {
  
  temp_number <- temp_number + 1
  
  first_transcript_id <- tibble_ref_gtf$transcript_id %>% unique %>% na.omit %>% .[temp_number]
  
}

# the test condition
max_test <- tibble_ref_gtf[tibble_ref_gtf$transcript_id == first_transcript_id, "exon_number"] %>% unlist %>% na.omit %>% max
max_exon_start_test <- tibble_ref_gtf[tibble_ref_gtf$transcript_id == first_transcript_id & tibble_ref_gtf$exon_number == max_test, "start"] %>% na.omit %>% paste
min_exon_start_test <- tibble_ref_gtf[tibble_ref_gtf$transcript_id == first_transcript_id & tibble_ref_gtf$exon_number == 1, "start"] %>% na.omit %>% paste

exon_order <- NULL

# if the exon 1 comes before the max exon, then the exon order is always increasing
if (min_exon_start_test < max_exon_start_test) {
  
  exon_order <- "increasing"
  
  # if the exon 1 cones after the max exon, then the exon order is stranded.
} else if (min_exon_start_test > max_exon_start_test) {
  
  exon_order <- "stranded"
  
}

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

cat("START PROCESSING\n")
# bind the mapping table onto the filtered split table
tibble_maxquant_table_filtered_and_split_with_mapping <- dplyr::left_join(
  tibble_maxquant_table_filtered_and_split,
  tibble_uniprot_to_ENSP_mapping_table,
  by = colname_mapping_info_table_containing_uniprot_id) %>%
  .[!is.na(.[, colname_mapping_info_table_containing_ensembl_peptide_id] %>% unlist), ] %>% unique

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# also write the tibble of entries that did not match 
tibble_maxquant_table_filtered_and_split_with_mapping_unmatched <- dplyr::left_join(
  tibble_maxquant_table_filtered_and_split,
  tibble_uniprot_to_ENSP_mapping_table,
  by = colname_mapping_info_table_containing_uniprot_id) %>%
  .[is.na(.[, colname_mapping_info_table_containing_ensembl_peptide_id] %>% unlist), ] %>% unique

# calculate CDS-relative start/end
tibble_added_CDS_relative_start <- tibble_maxquant_table_filtered_and_split_with_mapping %>%
  dplyr::mutate("CDS_relative_start" = 3*(start_position_within_protein %>% type.convert - 1) + 1,
                "CDS_relative_end" = 3*end_position_within_protein %>% type.convert)

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# calculate transcript-relative start/end
## to do this, we NEED the start codon position in the transcript.
## subset the GTF for only CDS entries with the ensembl peptide ids that we have
tibble_subset_GTF_for_ENSP <- tibble_ref_gtf[which((tibble_ref_gtf[, reference_gtf_field_name_info] %>% unlist) %in% tibble_added_CDS_relative_start$ensembl_peptide_id & 
                                                     tibble_ref_gtf$type == "CDS"), ]
## filtering join to obtain ENST ids.
tibble_added_transcript_ids <- dplyr::left_join(tibble_added_CDS_relative_start, 
                                          tibble_subset_GTF_for_ENSP[, c("transcript_id", reference_gtf_field_name_info)] %>% dplyr::rename(!!colname_mapping_info_table_containing_ensembl_peptide_id := reference_gtf_field_name_info) %>% unique,
                                          by = colname_mapping_info_table_containing_ensembl_peptide_id) %>%
  .[!is.na(.$transcript_id), ] %>% unique
  
if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# subset GTF for matched transcript_ids only
tibble_subset_GTF_for_ENST <- tibble_ref_gtf[which(tibble_ref_gtf$transcript_id %in% tibble_added_transcript_ids$transcript_id & 
                                                     tibble_ref_gtf$type %in% c("exon", "start_codon")), ]

# retrieve transcript-relative start codon positions as well as all the genome-relative positions of the transcript.
list_query_transcript_id_start_codon_and_genome_coord_info <- future_map(
  .x = tibble_added_transcript_ids$transcript_id,
  .f = function(a1) {
    
    # DEBUG ###
    # a1 <- tibble_added_transcript_ids$transcript_id %>% .[1]
    ###########
    
    # get GTF entry of transcript
    tibble_matched_transcript_entry <- tibble_subset_GTF_for_ENST[tibble_subset_GTF_for_ENST$transcript_id == a1, ]
    
    # get matched chr and strand
    current_chr <- tibble_matched_transcript_entry$seqnames %>% unique
    current_strand <- tibble_matched_transcript_entry$strand %>% unique
    
    # list out all the stranded genome-relative coords
    vector_stranded_genome_relative_nt_positions_of_transcript <- purrr::map2(.x = tibble_matched_transcript_entry %>% dplyr::filter(type == "exon") %>% .$start, 
                                                                              .y = tibble_matched_transcript_entry %>% dplyr::filter(type == "exon") %>% .$end, 
                                                                              .f = ~.x:.y) %>% unlist %>% unique %>% sort
      
    if (current_strand == "-") {
      vector_stranded_genome_relative_nt_positions_of_transcript <- vector_stranded_genome_relative_nt_positions_of_transcript %>% rev
    }
    
    # get stop codon position
    logical_start_codon_exists <- tibble_matched_transcript_entry %>% dplyr::filter(type == "start_codon") %>% nrow > 0
    
    if (current_strand == "+") {
      genome_relative_first_nt_of_start_codon <- tibble_matched_transcript_entry %>% dplyr::filter(type == "start_codon") %>% .$start %>% min
    } else if (current_strand == "-") {
      genome_relative_first_nt_of_start_codon <- tibble_matched_transcript_entry %>% dplyr::filter(type == "start_codon") %>% .$end %>% max
    }
    
    transcript_relative_first_nt_of_start_codon <- which(vector_stranded_genome_relative_nt_positions_of_transcript == genome_relative_first_nt_of_start_codon)
    
    # fill in start codon position with NA if there was no start codon entry in the GTF
    if (logical_start_codon_exists == FALSE) {
      transcript_relative_first_nt_of_start_codon <- NA
    }
    
    return(list(
      "chr" = current_chr,
      "strand" = current_strand,
      "vector_stranded_genome_relative_nt_positions_of_transcript" = vector_stranded_genome_relative_nt_positions_of_transcript %>% paste(collapse = ","),
      "transcript_relative_first_nt_of_start_codon" = transcript_relative_first_nt_of_start_codon
    ))
    
  }, .progress = TRUE)

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# now add the start codon info as new columns.
tibble_added_start_codon_info <- tibble_added_transcript_ids %>% 
  dplyr::mutate("chr" = list_query_transcript_id_start_codon_and_genome_coord_info %>% purrr::map(~.x$chr) %>% unlist,
                "strand" = list_query_transcript_id_start_codon_and_genome_coord_info %>% purrr::map(~.x$strand) %>% unlist,
                "vector_stranded_genome_relative_nt_positions_of_transcript" = list_query_transcript_id_start_codon_and_genome_coord_info %>% purrr::map(~.x$vector_stranded_genome_relative_nt_positions_of_transcript) %>% unlist,
                "transcript_relative_first_nt_of_start_codon" = list_query_transcript_id_start_codon_and_genome_coord_info %>% purrr::map(~.x$transcript_relative_first_nt_of_start_codon) %>% unlist) %>% 
  .[!is.na(.$transcript_relative_first_nt_of_start_codon), ] %>% unique

# now get transcript-relative start and end positions
tibble_added_transcript_relative_start_end <- tibble_added_start_codon_info %>%
  dplyr::mutate("transcript_relative_start" = `CDS_relative_start` + `transcript_relative_first_nt_of_start_codon` - 1,
                "transcript_relative_end" = `CDS_relative_end` + `transcript_relative_first_nt_of_start_codon` - 1)

# finally get all the genome-relative positions of the phosphosite/peptide etc...
tibble_added_genome_relative_coords_of_feature <- tibble_added_transcript_relative_start_end %>% 
  dplyr::mutate("all_genomic_positions_of_feature" = furrr::future_pmap(
    .l = list(
      "a1" = `transcript_relative_start`, 
      "a2" = `transcript_relative_end`, 
      "a3" = `vector_stranded_genome_relative_nt_positions_of_transcript`),
      .f = function(a1, a2, a3) {
        a3 %>% strsplit(split = ",") %>% unlist %>% type.convert %>% 
          .[a1:a2] %>% 
          paste(collapse = ",") %>% 
          return}, .progress = TRUE ) %>% unlist)

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# since there is a many-to-one uniprotkb to ENSP id mapping, we have to prune elements which have transcript-relative positions that lie outside of the mapped transcript length
list_added_genome_relative_coords_of_feature <- tibble_added_genome_relative_coords_of_feature %>% dplyr::group_split(chr) %>% purrr::reduce(bind_rows) %>% array_tree

# test for validity
vector_logical_invalid_transcript_start_end <- furrr::future_map(
  .x = tibble_added_genome_relative_coords_of_feature %>% dplyr::group_split(chr) %>% purrr::map(~.x %>% purrr::array_tree()), 
  .f = function(a1) {
    
    furrr::future_map(
      .x = a1, 
      .f = function(b1) {
        
        (b1$transcript_relative_start %>% type.convert > b1$vector_stranded_genome_relative_nt_positions_of_transcript %>% strsplit(split = ",") %>% unlist %>% length) | (b1$transcript_relative_end %>% type.convert > b1$vector_stranded_genome_relative_nt_positions_of_transcript %>% strsplit(split = ",") %>% unlist %>% length) %>% return
        
      } ) %>% unlist %>% return
    
  } , .progress = TRUE) %>% unlist

# prune 
list_added_genome_relative_coords_of_feature_pruned <- list_added_genome_relative_coords_of_feature[vector_logical_invalid_transcript_start_end == FALSE]

# convert to bedfile
list_feature_bedfile <- future_imap(.x = list_added_genome_relative_coords_of_feature_pruned, .f = function(a1, a2) {
  
  # DEBUG ###
  # a1 <- list_added_genome_relative_coords_of_feature_pruned[[108]]
  ###########
  
  # cat("\nnow processing:", a2, "/", length(tibble_added_genome_relative_coords_of_feature %>% array_tree))
  
  # vectorise the genome-relative feature positions
  vector_genome_relative_feature_positions <- a1$all_genomic_positions_of_feature %>% strsplit(split = ",") %>% unlist %>% type.convert %>% sort
  
  # convert the individual nt positions to a set of ranges
  ## achieve this by comparing n to n + 1
  ## for rows where diff > 1, the n represents the end of an exon. n.plus.1 represents the start of the exon right after the gap.
  tibble_n_n.plus.1 <- tibble("n" = vector_genome_relative_feature_positions %>% head(n = length(vector_genome_relative_feature_positions) - 1),
                              "n.plus.1" = vector_genome_relative_feature_positions[-1]) %>% 
    add_column("difference" = .$n.plus.1 - .$n)
  
  ## if there are no gaps, then just take the genomic range as the start:end
  if (tibble_n_n.plus.1$difference %>% unique %>% length == 1) {
    
    vec_blockCount <- 1
    vec_blockSizes <- length(vector_genome_relative_feature_positions)
    vec_blockStarts <- first(vector_genome_relative_feature_positions)
    
  } else if (tibble_n_n.plus.1$difference %>% unique %>% length > 1) {
    
    vec_blockCount <- tibble_n_n.plus.1$difference %>% unique %>% length
    vec_blockStarts <- c(first(vector_genome_relative_feature_positions), 
                         tibble_n_n.plus.1[tibble_n_n.plus.1$difference > 1, "n.plus.1"] %>% unlist(use.names = FALSE))
    vec_blockEnds <- c(tibble_n_n.plus.1[tibble_n_n.plus.1$difference > 1, "n"] %>% unlist(use.names = FALSE), 
                       last(vector_genome_relative_feature_positions))
    vec_blockSizes <- vec_blockEnds - vec_blockStarts + 1
    
  }
  
  # the block details are still vector at this point. 
  # must paste collapse in order to use in the BED file.
  tibble_block_info <- tibble(
    "chrom" = a1$chr,
    "start" = first(vector_genome_relative_feature_positions),
    "end" = last(vector_genome_relative_feature_positions),
    "name" = paste(a1$uniprotkb_entry, a1$`second_id`, sep = "|"),
    "score" = a1$Score %>% type.convert,
    "strand"= a1$strand,
    "thickStart" = first(vector_genome_relative_feature_positions),
    "thickEnd" = last(vector_genome_relative_feature_positions),
    "itemRgb" = track_colour,
    "blockCount" = vec_blockCount %>% paste(collapse = ","),
    "blockSizes" = vec_blockSizes %>% paste(collapse = ","),
    "blockStarts" = (vec_blockStarts - first(vector_genome_relative_feature_positions)) %>% paste(collapse = ","),
    "qName" = a1$`second_id` %>% paste %>% trimws)
  
}, .progress = TRUE)

if (save_workspace_when_done == "DEBUG") {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# rbind and tibblise
tibble_feature_bedfile <- list_feature_bedfile %>% rbindlist(use.names = TRUE, fill = TRUE) %>% as_tibble %>% dplyr::select(-qName)

# track line
# cat(paste("track name=\"", output_file_name_JUM, "\" description=\"", output_file_name_JUM, "\" graphType=junctions\n", sep = ""), file = paste(R_processing_results_dir, "proteomic_alignment_", a2, "_", b2, ".bed", sep = ""))
# 
# append = TRUE
write.table(x = tibble_feature_bedfile, file = paste(output_dir, "/", output_name, ".bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# write supplementary info 
write.table(x = tibble_added_genome_relative_coords_of_feature, file = paste(output_dir, "/", output_name, "_supp_info.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

if (save_workspace_when_done %in% c("DEBUG", "YES")) {
  save.image(file = paste(output_dir, "/", output_name, "_workspace.RData", sep = ""))
}

# finish counting
tictoc::toc()

q()

