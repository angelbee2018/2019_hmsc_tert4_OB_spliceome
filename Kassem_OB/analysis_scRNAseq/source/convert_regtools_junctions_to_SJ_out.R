script_description <- "A tiny script to convert regtools junctions to SJ.out.tab format like in STAR: http://www.trii.org/courses/rnaseq_course_materials/STARmanual.pdf"

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

# library(seqinr)
## can replicate the intronic motif from STAR, but not today.
library(tidyverse)
# library(rtracklayer)
library(data.table)
library(optparse)

library(tictoc)

# manage arguments
list_input_arg_info = list(
  "1" = make_option(c("-F", "--file"), type = "character", default = NULL, 
                    help = "Compulsory. Path to the raw output by Regtools, the junctions extracted from a BAM file.", metavar = "character"),
  "2" = make_option(c("-O", "--output_dir"), type = "character", default = NULL, 
                    help = "Compulsory. Directory where the output file should be placed.", metavar = "character"),
  "3" = make_option(c("-N", "--name"), type = "character", default = NULL, 
                    help = "Compulsory. Name of the file to be outputted.", metavar = "character")
)

input_arg_info <- OptionParser(option_list = list_input_arg_info, description = script_description)
input_args <- input_arg_info %>% parse_args

# check if the input arguments are O.
if (list(input_args$file, input_args$output_dir, input_args$name) %>% lapply(is.null) %>% unlist %>% any == TRUE) {
  
  print_help(input_arg_info)
  
  stop("Make sure you entered the arguments correctly", call. = FALSE)
  
}

# DEBUG #######
# input_file_path <- "Z:/PGNEXUS_kassem_MSC/Kassem_OB/analysis_scRNAseq/source/cluster_1_pseudorep3_SJ.out.tab"
# output_dir <- "Z:/PGNEXUS_kassem_MSC/Kassem_OB/analysis_scRNAseq/source/"
# output_name <- "test"
###############

input_file_path <- input_args$file
output_dir <- input_args$output_dir
output_name <- input_args$name

cat("input_file_path:", input_file_path, "\n")
cat("output_dir:", output_dir, "\n")
cat("output_name:", output_name, "\n")

if(!dir.exists(output_dir) ) {
  dir.create(output_dir, recursive = TRUE)}

tictoc::tic(msg = paste("Time elapsed for: ", input_args$name, sep = ""))

tibble_regtools_junctions <- read.delim(file = input_file_path, sep = "\t", header = FALSE, row.names = NULL, stringsAsFactors = FALSE) %>% 
  setNames(c("chr", "start", "end", "junc_id", "read_count", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")) %>% 
  as_tibble

print(head(tibble_regtools_junctions))

cat("\nconvert strand into STAR format - anything that's a question mark is assumed to be on the reverse strand.")
tibble_STAR.like_output_temp <- tibble_regtools_junctions

tibble_STAR.like_output_temp[, "strand"] <- tibble_STAR.like_output_temp$strand %>% 
  purrr::map(.f = function(a1) {
    if (a1 == "+") {
      return("1")
    } else if (a1 == "?") {
      return("2")
    }
  } ) %>% unlist

cat("\nchange the start and end to actually describe a junction.")
tibble_STAR.like_output_temp[, "start"] <- tibble_STAR.like_output_temp$start + (tibble_STAR.like_output_temp$blockSizes %>% strsplit(split = ",") %>% purrr::map(~.x[1]) %>% unlist %>% as.numeric)
tibble_STAR.like_output_temp[, "end"] <- tibble_STAR.like_output_temp$end - (tibble_STAR.like_output_temp$blockSizes %>% strsplit(split = ",") %>% purrr::map(~.x[2]) %>% unlist %>% as.numeric)

cat("\nadd the maximum spliced alignment overhang.")
tibble_STAR.like_output_temp <- add_column(tibble_STAR.like_output_temp, "max_overhang" = tibble_STAR.like_output_temp$blockSizes %>% strsplit(split = ",") %>% purrr::map(~max(.x %>% as.numeric)) %>% unlist)

print(head(tibble_STAR.like_output_temp))

print(colnames(tibble_STAR.like_output_temp))

cat("\ncreate final table")
tibble_STAR.like_output_final <- tibble_STAR.like_output_temp[, c("chr", "start", "end", "strand", "read_count", "max_overhang")] %>%
  add_column("intron_motif" = 0,
             "SJ_annotation" = 0, .after = "strand") %>%
  add_column("multi_mapping_read_count" = 0, .after = "read_count")

cat("\nwrite the final output table")
write.table(x = tibble_STAR.like_output_final,
            file = paste(output_dir, output_name, "_SJ.out.tab", sep = ""),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

tictoc::toc()

q()


