## FIRST, AN R SCRIPT TO GENERATE THE COMBINATION NAMES OF THE PAIRWISE COMPARISONS BETWEEN TIMEPOINTS.

# arg1: $sampleandaccessionlists the directory where the sample and accession lists are arg2: $results_dir the directory to which the end file is to be outputted arg3: $timepoint_name_list the list of timepoints

print(commandArgs())

args = commandArgs(trailingOnly=TRUE)

print(args)

print(paste("number of arguments: ", length(args), sep = ""))

# test if there are 3 arguments: if not, return an error
if (length(args)!=3) {
  stop("check the amount of arguments you have. should have specified 3", call.=FALSE)
}

sample_name_list <- read.delim(paste(args[3], sep = ""), sep = "\t", stringsAsFactors = FALSE, header = FALSE)

print(sample_name_list)

raw_comparison_table <- t(combn(sample_name_list[, 1], 2))

print(raw_comparison_table)

raw_comparison_table <- cbind(raw_comparison_table, paste(raw_comparison_table[, 1], "_vs_", raw_comparison_table[, 2], sep = ""))

# we now extract the list of comparisons between ALL samples IN COLUMN 3

list_of_timepoint_comparisons_final <- raw_comparison_table[, 3]

print(list_of_timepoint_comparisons_final)

print(paste(length(list_of_timepoint_comparisons_final), " comparisons will be considered.", sep = ""))

write.table(x = list_of_timepoint_comparisons_final, file = paste(args[2], "list_of_timepoint_comparisons_final.txt", sep = ""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

q()