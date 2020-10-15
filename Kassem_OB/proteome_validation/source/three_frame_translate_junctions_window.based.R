script_description <- "# THREE-FRAME TRANSLATION ######
Attempts to translate the nucleotides of all transcripts flanking all the splice junctions given.
Uses PARALLEL PURRR (FURRR) ^___^ is much faster than the last version
Recommended CPU/RAM connsumption: 6 cores/32GB for ~30,000 junctions. 4 cores/24GB for ~5000 junctions. 2 cores/16GB for ~1000 junctions. 
I personally would not use more than 12 cores because you get diminishing returns due to the longer time it takes to delegate the tasks to each worker.
RAM usage  scales by the number of cores you use. I do not recommend going lower than 16GB."

# arg1: path to a table containing the junction, reconstructed GTF file and desired output name information. should be in the format of a tab separated tibble to minimise the risk of error. each row describes A PATH the file to the junction table and reconstr. GTF to be translated. This program will automatically loop thru them. The column names must be: c("reconstructed_GTF_path", "junction_table_path", "output_database_name")
# junction_table_path: path to table containing junctions of interest. for example, UNION_junc_coor type files. MUST contain columns start, end, chr, strand, and fasta_header. one row per junction. should be in the form of a tibble. the fasta header you want included in the final FASTA file needs to be in the column called "fasta_header".
# reconstructed_GTF_path: path to the actual reconstructed GTF file (e.g. from Cufflinks, Strawberry). NOT THE CONTAINING DIRECTORY. tip: for better junction matching, combine the reconstructed GTF with reference GTF beforehand e.g. using StringTie
# output_database_name: a character string of what the final FASTA database file name and output table will be called

# arg2: path to the directory containing the genome FASTA files. Ideally from Ensembl... you need separate files by chromosomes, NOT the primary assembly. FORMATTING IMPORTANT!!!! MUST BE IN THE FORMAT: <_anything_><chr>.fa e.g. "Homo_sapiens.GRCh38.dna.chromosome.MT.fa" OR "chr16.fa" OR "Y.fa". What will not work: anything which does not have .fa extension e.g. "chr16.fasta", anything between the chromosome number and the .fa extension e.g. "chromosome1.ensembl.fa"
# arg3: how many nucleotides to translate upstream W.R.T. the transcript. for junction-based mode, the total nt. length to be translated will be arg3 + arg4. default for both is 50.
# arg4: how many nt to translate downstream W.R.T. the transcript. for exon-based mode, this will be 0 by default. the total arg3 + arg4 will indicate how many nt to translate beyond the exon.
# arg5: output directory
# arg6: number of cores to use. possible inputs: numbers 1 to any integer
# arg7: chromosomes to do (chrmode): select what chromosomes you want translated. possible inputs: numbers 0-2. 0 (default): nuclear chromosomes only i,e, 1:22, X & Y. 1: nuclear + mitochondrial i.e. 1:22, X & Y, M. 2: everything including haplotype/fusion chromosomes etc... this is possible provided the chromosome names 
# arg8: nonchromosomal file name. if you are doing haplotypes, please specify what the reference genome FASTA file for it is called or the script won't know. This single FASTA file must contain all the haplotype information. The script won't try to search for a second file.  In ensembl, this file is called "Homo_sapiens.GRCh38.dna.nonchromosomal.fa" or more generally, "*nonchromosomal.fa". So for this option, you would specify "--nonchrname nonchromosomal".

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
              help = "Compulsory. path to a table containing the junction, reconstructed GTF file and desired output name information. should be in the format of a tab separated tibble to minimise the risk of error. each row describes A PATH the file to the junction table and reconstr. GTF to be translated. This program will automatically loop thru them. The column names must be: c(\"reconstructed_GTF_path\", \"junction_table_path\", \"output_database_name\").
              - junction_table_path: path to table containing junctions of interest. for example, UNION_junc_coor type files. MUST contain columns start, end, chr, strand, and fasta_header. one row per junction. should be in the form of a tibble. the fasta header you want included in the final FASTA file needs to be in the column called \"fasta_header\".
              - reconstructed_GTF_path: path to the actual reconstructed GTF file (e.g. from Cufflinks, Strawberry). NOT THE CONTAINING DIRECTORY. tip: for better junction matching, combine the reconstructed GTF with reference GTF beforehand e.g. using StringTie
              - output_database_name: a character string of what the final FASTA database file name and output table will be called", metavar = "character"),
  "2" = make_option(c("-R", "--reference_genome_fasta_dir"), type = "character", default = NULL, 
              help = "Compulsory. path to the directory containing the genome FASTA files. Ideally from Ensembl... you need separate files by chromosomes, NOT the primary assembly. 
              FORMATTING IMPORTANT!!!! MAKE SURE THE REF. GENOME FASTA FILES ARE IN THE FORMAT: <_anything_><chr>.fa e.g. \"Homo_sapiens.GRCh38.dna.chromosome.MT.fa\" OR \"chr16.fa\" OR \"Y.fa\". What will not work: anything which does not have .fa extension e.g. \"chr16.fasta\", anything between the chromosome number and the .fa extension e.g. \"chromosome1.ensembl.fa\"", metavar = "character"),
  "3" = make_option(c("-U", "--upstream_window_size"), type = "double", default = 50, 
              help = "Optional. how many nucleotides to translate upstream W.R.T. the transcript. for junction-based mode, the total nt. length to be translated will be arg3 + arg4. default for both is 50.", metavar = "double"),
  "4" = make_option(c("-D", "--downstream_window_size"), type = "double", default = 50, 
              help = "Optional. how many nt to translate downstream W.R.T. the transcript. for exon-based mode, this will be 0 by default. the total arg3 + arg4 will indicate how many nt to translate beyond the exon.", metavar = "double"),
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
reference_genome_fasta_dir <- input_args$reference_genome_fasta_dir
upstream_window_size <- input_args$upstream_window_size
downstream_window_size <- input_args$downstream_window_size
output_dir <- input_args$output_dir

# manage parrallellisation rrlllRll

if (input_args$ncores != 0) {
  number_of_workers <- input_args$ncores
} else {
  number_of_workers <- future::availableCores()
} 

# LINUX STARTS COUNTING THE NUMBER OF CORES STARTING FROM 0. FIX THIS.
if(Sys.info()["sysname"] == "Windows") {
  
  cat("number of workers:", number_of_workers, "\n")
  
  future::plan(multiprocess)
  options(future.globals.maxSize = 30000000000, mc.cores = input_args$ncores)
  
} else {
  
  cat("number of workers:", number_of_workers, "\n")
  
  future::plan(multisession, workers = number_of_workers)
  # library(RhpcBLASctl)
  # RhpcBLASctl::omp_set_num_threads(number_of_workers)
  # RhpcBLASctl::blas_set_num_threads(number_of_workers)
  options(future.globals.maxSize = 30000000000, mc.cores = number_of_workers)
  
}

cat(future::availableCores(), "cores will be used\n")

# testing ########

# file_information_path <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/input_files/file_info_test.txt"
# reconstructed_gtf_path <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/GRAND_OBseries_ref_reconstructed_stringtiemerged.gtf"
# reference_genome_fasta_path <- "/srv/scratch/z3463471/hg38_ensembl_reference/raw_genome_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# output_dir <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/results_database_generation/"
# output_database_name <- "three_frame_translation_text"

# DEFINE FUNCTIONS ##########################

# FUNCTION TO RETURN A VECTOR OF ALL GENOMIC POSITIONS OCCURRING IN THE SPECIFIED WINDOW, SCANNING ALONG ONE TRANSCRIPT
# TO BE USED WITH PURRR - THIS PROCESSES THE TRANSLATION START/END COORDINATES FOR ONE TRANSCRIPT ONLY
# transcript.upstream_window: amount of nucleotides you want to consider in the 5' direction WITH RESPECT TO THE TRANSCRIPT (NOT THE GENOME!!)
# transcript.downstream_window: exactly the same but for 3' direction. the 2 nucleotide positions directly flanking the junction are included in this up/downstream window
# in contrast, up/downstream_scan_range: scans from the position of the theoretical median nucleotide between the two junction-flanking nucleotides. has half-integer value.

generate_all.genomic.positions_for_translation <- function(list_junction.adjacent_exon_start.end_coords_one_transcript, vector_all.genomic.positions_one_transcript, transcript.upstream_window = 50, transcript.downstream_window = 50) {
  
  
  # DEBUG ###################
  
  # transcript.upstream_window <- 50
  # transcript.downstream_window <- 50
  # 
  # list_junction.adjacent_exon_start.end_coords_one_transcript <- list_of_junction.adjacent_exon_start.end_coords_per_transcript[[2]]
  # vector_all.genomic.positions_one_transcript <- list_all.genomic.positions_per_transcript[[2]]
  
  ###########################
  
  nt_length_of_transcript <- length(vector_all.genomic.positions_one_transcript)
  
  upstream_scan_range <- transcript.upstream_window - 0.5
  downstream_scan_range <- transcript.downstream_window - 0.5
  
  transcript.upstream_flanking.nt.pos <- which(vector_all.genomic.positions_one_transcript == list_junction.adjacent_exon_start.end_coords_one_transcript$end_transcript.5prime_exon_coord)
  
  transcript.downstream_flanking.nt.pos <- which(vector_all.genomic.positions_one_transcript == list_junction.adjacent_exon_start.end_coords_one_transcript$start_transcript.3prime_exon_coord)
  
  transcript.median.nt.pos <- ((transcript.upstream_flanking.nt.pos %>% as.double) + (transcript.downstream_flanking.nt.pos %>% as.double)) / 2
  
  # CHECKPOINT - IF THE FLANKING POSITIONS ARE NOT DIRECTLY ADJACENT IN THE SPLICED TRANSCRIPT, THROW AN ERROR AND STOP.
  
  if (abs(transcript.upstream_flanking.nt.pos - transcript.downstream_flanking.nt.pos) != 1) {
    
    stop("ERROR IN FUNCTION: \"generate_all.genomic.positions_for_translation\": junction-flanking nucleotide positions are not adjacent in the mature spliced transcript")
    
  }
  
  # return vector of genomic positions to translate
  
  # translate as much as possible, up to the bounds of the transcript.
  
  if (list_junction.adjacent_exon_start.end_coords_one_transcript$sign == 1) {
    
    vector_nt.positions_to_translate <- vector_all.genomic.positions_one_transcript[max((transcript.median.nt.pos - upstream_scan_range), 0):min((transcript.median.nt.pos + downstream_scan_range), nt_length_of_transcript)]
    
  } else if (list_junction.adjacent_exon_start.end_coords_one_transcript$sign == -1) {
    
    vector_nt.positions_to_translate <- vector_all.genomic.positions_one_transcript[min((transcript.median.nt.pos + upstream_scan_range), nt_length_of_transcript):max((transcript.median.nt.pos - downstream_scan_range), 0)]
    
  }
  
  return(list("vector_genome.coords_to_translate" = vector_nt.positions_to_translate, "chr" = list_junction.adjacent_exon_start.end_coords_one_transcript$chr %>% trimws, "strand" = list_junction.adjacent_exon_start.end_coords_one_transcript$strand %>% trimws))
  
}

# END generate_all.genomic.positions_for_translation

# FUNCTION TO 3 FRAME TRANSLATE ONE LIST CONTAINING NUCLEOTIDE SEQUENCE AND STRAND

nt.sequence_strand_threeframetranslate <- function(list_nt.fwd.sequence_strand) {
  
  if (list_nt.fwd.sequence_strand$strand == "+") {
    
    translation_result <- list("translation_frame_0" = seqinr::translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 0, sens = "F"),
                               "translation_frame_1" = seqinr::translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 1, sens = "F"),
                               "translation_frame_2" = seqinr::translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 2, sens = "F"))
    
  } else if (list_nt.fwd.sequence_strand$strand == "-") {
    
    translation_result <- list("translation_frame_0" = seqinr::translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 0, sens = "R"),
                               "translation_frame_1" = seqinr::translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 1, sens = "R"),
                               "translation_frame_2" = seqinr::translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 2, sens = "R"))
    
  }
  
  return(translation_result)
  
}

# END nt.sequence_strand_threeframetranslate

# MAIN FUNCTION TO DO THREE FRAME TRANSLATION
# this function will find all flanking exons from ONE input splice junction range by matching with the reconstructed GTF file, then translate in 3 frames.
# INPUT: a) spliceregion_list, a list (e.g. a single element from a spliceregion_table %>% array tree) containing $chr, $start, $end, $strand, $junction_ID, b) reconstructed GTF table from rtracklayer, c) path to the reference genome fasta itself
# OUTPUT: a list of every transcript. inside, amino acid sequence from translation of 3 frames
# mode: can be one of the following: junction OR exon.
# junction mode e.g. JUM: spliceregion_table must contain columns: chr, start, end and strand of the splice junction region. one row per junction.
# exon mode e.g. Whippet, PSI-Sigma: spliceregion_table must contain columns: chr, start, end and strand of the differentially included exons. one row per exon.
# loop_number: used to print out the progress via purrr::imap (prints which junction it's currently looping through)

three_frame_translate_splicejunctions <- function(spliceregion_list, loop_number, reconstructed_gtf_table, reference_genome_fasta = reference_genome_fasta, mode = NULL, transcript.upstream_window = 50, transcript.downstream_window = 50) {
  
  # DEBUG ###################
  
  # transcript.upstream_window <- 50
  # transcript.downstream_window <- 50
  # 
  # reconstructed_gtf_table <- reconstructed_gtf
  # spliceregion_list <- list_junction_table_array.tree_temp[[1]]
  
  ###########################
  
  # filter the reconstructed GTF table for all exon entries that directly flank the splice junction
  
  cat("Now processing junction number:", loop_number, "\n")
  
  cat("finding which transcripts contain the specified junction...\n")
  
  if (spliceregion_list$strand == ".") {
    
    tibble_reconstructed_gtf_subset_flanking_exons <- reconstructed_gtf_table[reconstructed_gtf_table$seqnames == spliceregion_list$chr %>% trimws, ] %>% .[.$start <= ((spliceregion_list$end %>% as.numeric) + 2) & .$end >= ((spliceregion_list$start %>% as.numeric) - 2), ] %>% .[!(.$start <= ((spliceregion_list$end %>% as.numeric) - 2) & .$end >= ((spliceregion_list$start %>% as.numeric) + 2)), ] %>% .[.$type == "exon", ]
    
  } else if (spliceregion_list$strand == "+" | spliceregion_list$strand == "-") {
    
    tibble_reconstructed_gtf_subset_flanking_exons <- reconstructed_gtf_table[reconstructed_gtf_table$seqnames == spliceregion_list$chr %>% trimws, ] %>% .[.$strand == spliceregion_list$strand %>% trimws, ] %>% .[.$start <= ((spliceregion_list$end %>% as.numeric) + 2) & .$end >= ((spliceregion_list$start %>% as.numeric) - 2), ] %>% .[!(.$start <= ((spliceregion_list$end %>% as.numeric) - 2) & .$end >= ((spliceregion_list$start %>% as.numeric) + 2)), ] %>% .[.$type == "exon", ]
    
  } else {
    
    stop("Could not match the strand information in the transposed differential-only UNION_junc_coor_table. Make sure that the \"strand\" column in the UNION_junc_coor_table contains only +, - or .")
    
  }
  
  list_of_junction_associated_transcripts <- tibble_reconstructed_gtf_subset_flanking_exons$transcript_id %>% unique %>% array_tree %>% flatten
  
  # make a list for each transcript that directly flanks a junction.
  # then filter so that there are only a) exon PAIRS which b) are directly connected in the mature (spliced) transcript
  
  list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_reconstructed_gtf_subset_flanking_exons[tibble_reconstructed_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts) %>% keep(.x = ., .p = ~nrow(.x) == 2) %>% keep(.x = ., .p = ~abs((.x[2, "exon_number"] %>% paste %>% as.numeric) - (.x[1, "exon_number"] %>% paste %>% as.numeric)) == 1)
  
  #### these are the two main tables to be used for translation #########
  
  # to account for + or - strand transcription, the positive "sign" indicates the second exon has higher genome coordinate than the first. negative "sign" indicates the opposite.
  
  cat("defining junction-flanking nucleotide positions...\n")
  
  list_of_junction.adjacent_exon_start.end_coords_per_transcript <- purrr::map(.x = list_of_tibbles_flanking_exon_gtf.entries_per_transcript, .f = ~list("end_transcript.5prime_exon_coord" = .x[which(.x$exon_number %>% as.numeric == .x$exon_number %>% as.numeric %>% min), "end"] %>% paste %>% as.numeric, "start_transcript.3prime_exon_coord" = .x[which(.x$exon_number %>% as.numeric == .x$exon_number %>% as.numeric %>% max), "start"] %>% paste %>% as.numeric, "chr" = .x[1, "seqnames"] %>% paste, "strand" = .x[1, "strand"] %>% paste) %>% splice("sign" = (.$start_transcript.3prime_exon_coord - .$end_transcript.5prime_exon_coord)/abs(.$start_transcript.3prime_exon_coord - .$end_transcript.5prime_exon_coord)))
  
  cat("define genome-relative co-ordinates for each nucleotide in associated transcripts...\n")
  
  list_all.genomic.positions_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~reconstructed_gtf_table[reconstructed_gtf_table$transcript_id == .x & reconstructed_gtf_table$type == "exon", c("start", "end")] %>% array_tree %>% purrr::map(.x  = ., .f = ~.x[[1]]:.x[[2]]) %>% unlist %>% unique %>% sort) %>% set_names(list_of_junction_associated_transcripts) %>% .[names(list_of_tibbles_flanking_exon_gtf.entries_per_transcript)]
  
  #######################################################################
  
  # generate genome-relative coordinates flanking each junction for translation
  
  cat("define genomic positions of nucleotides to be translated..\n")
  
  list_genome.coords_for_translation <- purrr::map2(.x = list_of_junction.adjacent_exon_start.end_coords_per_transcript, .y = list_all.genomic.positions_per_transcript, .f = ~generate_all.genomic.positions_for_translation(.x, .y, transcript.upstream_window = transcript.upstream_window, transcript.downstream_window = transcript.downstream_window))
  
  cat("looking up reference FASTA to fetch nucleotides within specified transcript window...\n")
  
  list_forward_nucleotides_from_coords <- purrr::map(.x = list_genome.coords_for_translation, .f = ~list("forward_nucleotides" = reference_genome_fasta[[.x$chr %>% trimws]] %>% .[.x$vector_genome.coords_to_translate], "strand" = .x$strand))
  
  cat("initiating three-frame translation...\n")
  
  # three-frame translation
  
  list_threeframetranslate <- purrr::map(.x = list_forward_nucleotides_from_coords, .f = ~nt.sequence_strand_threeframetranslate(.x))
  
  list_threeframetranslate_fwd.nucleotides.genome.coords <- purrr::pmap(list(list_threeframetranslate, 
                                                                             list_forward_nucleotides_from_coords, 
                                                                             list_genome.coords_for_translation), 
                                                                        .f = ~list("threeframetranslate" = ..1, 
                                                                                   "forward_nucleotides_from_coords" = ..2["forward_nucleotides"] %>% unlist %>% paste(collapse = ","), 
                                                                                   "genome.coords_for_translation" = ..3["vector_genome.coords_to_translate"] %>% unlist %>% paste(collapse = ","),
                                                                                   "fasta_header" = spliceregion_list$fasta_header %>% paste))
  
  cat("...Success!!\n")
  
  return(list_threeframetranslate_fwd.nucleotides.genome.coords)
  
}

# END three_frame_translate_splicejunctions

# FUNCTION TO TRIM VIRTUAL PEPTIDES

# SOME RULES:
#   
# - STOP CODON CANNOT BE ON BOTH SIDES OF THE SPLICE JUNCTION
# - IF A STOP CODON IS DOWNSTREAM OF SPLICE JUNCTION, THE (AT LEAST) 16 AA BEFORE IT MUST HAVE BEEN ABLE TO BE TRANSLATED. THEREFORE, FILTER OUT LENGTH < 16 AA
# - IF ONE OR MORE STOP CODONS ARE UPSTREAM OF SJ, THE LAST STOP CODON MUST BE FOLLOWED BY A METHIONINE BEFORE THE SPLICE JUNCTION. it is entirely possible that a splice junction can be upstream of the start codon. BUT THAT MEANS THE SPLICED PEPTIDE ISNT SPANNING, IS IT
# - THEREFORE START CODON CAN BE, AT THE VERY LATEST, THE 16TH AA DUE TO THE DIVISION OF WINDOW LENGTH BY 3 FOR EACH FRAME.
# - for filtering intents and purposes, the splice junction should be considered the boundary of the 16th and 17th AA.
#
# BAHAHAHAHAHA I DON'T HAVE TO TAKE THE STRAND INTO ACCOUNT!! WOOP WOOP
# 
# THE SHORTEST HUMAN PROTEIN KNOWN IS 44 AA. THERE IS NO EXCUSE
#
# input: a vector of freshly in-silico translated peptides e.g. from seqinR, with one AA per element.
# output: either the successfully trimmed peptide in a string or a string if unsuccessful: "none".

trim_virtual_peptides <- function(vector_raw_virtual_peptide) {
  
  # DEBUG ######################################
  # vector_raw_virtual_peptide <- a[[2]][[2]][["threeframetranslate"]][[2]]
  # upstream_window_size <- 50
  # downstream_window_size <- 50
  ##############################################
  
  # divide the peptide into upstream and downstream
  upstream_segment <- vector_raw_virtual_peptide[1:(upstream_window_size %/% 3)]
  downstream_segment <- vector_raw_virtual_peptide[((upstream_window_size %/% 3) + 1):length(vector_raw_virtual_peptide)]
  
  # test conditions
  # stop codon on both sides
  if (("*" %in% upstream_segment == TRUE) & ("*" %in% downstream_segment == TRUE)) {
    
    string_trimmed_virtual_peptide <- "none"
    
    # stop codon in the downstream half only
  } else if (("*" %in% upstream_segment == FALSE) & ("*" %in% downstream_segment == TRUE)) {
    
    string_trimmed_virtual_peptide <- vector_raw_virtual_peptide %>% paste(collapse = "") %>% strsplit(split = "\\*") %>% unlist %>% .[1]
    
    # test for length > half the virtual peptide length
    if (nchar(string_trimmed_virtual_peptide) < upstream_window_size %/% 3) {
      string_trimmed_virtual_peptide <- "none"
    }
    
    # stop codon in the upstream half only
  } else if (("*" %in% upstream_segment == TRUE) & ("*" %in% downstream_segment == FALSE)) {
    
    # test whether there is an "M" before the end of the first segment after the last stop codon
    if ("M" %in% upstream_segment[last(which(upstream_segment == "*")) : (upstream_window_size %/% 3)]) {
      
      # if so, then the final peptide is going to be everything after the last stop codon.
      string_trimmed_virtual_peptide <- vector_raw_virtual_peptide[(last(which(upstream_segment == "*")) + 1) : length(vector_raw_virtual_peptide)] %>% paste(collapse = "")
      
    } else {
      
      # if not, then bad luck. it's thrown out.
      string_trimmed_virtual_peptide <- "none"
      
    }
    
    # no stop codon at all
  } else if ("*" %in% vector_raw_virtual_peptide == FALSE) {
    
    string_trimmed_virtual_peptide <- vector_raw_virtual_peptide %>% paste(collapse = "")
    
  }
  
  return(string_trimmed_virtual_peptide)
  
}

# END trim_virtual_peptides ###

# END FUNCTIONS #####################################################################################################
#####################################################################################################################
#####################################################################################################################

# BEGIN EXECUTION #################################

# DEBUG ################

# reconstructed_gtf_path <- "Z:/PGNEXUS_kassem_MSC/Kassem_OB_fastqc/analysis_strawberry/results_assemblyonly/merged/GRAND_OBseries_ref_denovo_reconstructed_stringtiemerged.gtf"
# reconstructed_gtf <- rtracklayer::import(reconstructed_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
# 
# reference_genome_fasta_dir <- "Z:/hg38_ensembl_reference/raw_genome_fasta/genome_fasta_extract2/"
# 
# junction_table_path <- "Z:/PGNEXUS_kassem_MSC/Kassem_OB_fastqc/proteome_validation/results_database_generation/angel_3FT_junctions/junction_table_OB_diff_any_qvalue0.01_any_deltaPSI_greaterthan_0.2.txt"
# 
# tibble_junction_table <- read.delim(junction_table_path, sep = "\t", stringsAsFactors = FALSE) %>% as_tibble
# upstream_window_size <- 50
# downstream_window_size <- 50
########################

vector_ref_genome_paths_by_chr <- paste(reference_genome_fasta_dir, list.files(reference_genome_fasta_dir)[list.files(reference_genome_fasta_dir) %>% grep(., pattern = ".*.fa$")], sep = "")

file_information_table <- read_delim(file_information_path, delim = "\t")

# reference_genome_fasta <- seqinr::read.fasta(file = reference_genome_fasta_path, forceDNAtolower = FALSE)

total_peptide_window_size <- upstream_window_size + downstream_window_size

cat("total translation window size:", total_peptide_window_size, "\n")

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
  
  list_junction_table_array.tree <- tibble_junction_table %>% array_tree
  
  # main execution step of three-frame translation
  # list_three_frame_translate_result <- future_imap(.x = list_junction_table_array.tree, .f = ~three_frame_translate_splicejunctions(.x, .y, reconstructed_gtf_table = reconstructed_gtf, reference_genome_fasta = reference_genome_fasta, mode = NULL, transcript.upstream_window = upstream_window_size, transcript.downstream_window = downstream_window_size), .progress = TRUE, .options = future_options(globals = c("reference_genome_fasta", "reconstructed_gtf")))
  
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
  
  # i couldn't do this any other way. This way i loop thru each chromosome and economically subset the junction table into smaller, more manageable tables.
  for (chr in chr_to_run) {
    
    # start counting
    tictoc::tic("Last chromosome")
    
    cat("now running chromosome ...", chr, "\n")
    
    # get positions of the vector where the path of the ref. genome fasta file is
    vector_ref_genome_paths_by_chr_position <- grep(x = vector_ref_genome_paths_by_chr, pattern = paste("(\\D|^)", chr, ".fa$", sep = ""))
    
    # temporary allocation to ref genome fasta list
    reference_genome_fasta_chr_temp <- seqinr::read.fasta(file = paste(vector_ref_genome_paths_by_chr[vector_ref_genome_paths_by_chr_position]), forceDNAtolower = FALSE)
    
    # temporary subsetting of reconstructed GTF to reduce overhead when exporting
    reconstructed_gtf_temp <- reconstructed_gtf[reconstructed_gtf$seqnames == chr, ]
    
    # subset the junction table for the current chromosome and array_tree it for furrr
    list_junction_table_array.tree_temp <- tibble_junction_table[tibble_junction_table$chr == chr, ] %>% array_tree
    
    assign(x = paste("list_three_frame_translate_result_temp_chr_", chr, sep = ""), 
           value = future_imap(.x = list_junction_table_array.tree_temp, .f = ~three_frame_translate_splicejunctions(.x, .y, reconstructed_gtf_table = reconstructed_gtf_temp, reference_genome_fasta = reference_genome_fasta_chr_temp, mode = NULL, transcript.upstream_window = upstream_window_size, transcript.downstream_window = downstream_window_size), .progress = TRUE, .options = future_options(globals = c("reconstructed_gtf_temp", "list_junction_table_array.tree_temp", "reference_genome_fasta_chr_temp", "three_frame_translate_splicejunctions", "generate_all.genomic.positions_for_translation", "nt.sequence_strand_threeframetranslate", "trim_virtual_peptides", "upstream_window_size", "downstream_window_size"))))
    
    tictoc::toc()
    
  }
  
  tictoc::toc()
  
  # make a list combining the 3FT results of each chromosome
  # trim off all those elements (per chromosome) that did not have any 3FT result.
  list_three_frame_translate_result <- ls(pattern = "list_three_frame_translate_result_temp_chr_") %>% array_tree %>% purrr::map(.x = ., .f = ~get(.x)) %>% compact
  # flatten into the correct format - each level 1 element pertains to ONE junction.
  list_three_frame_translate_result <- list_three_frame_translate_result %>% flatten
  
  # DEBUG ###########################################
  
  # test <- three_frame_translate_splicejunctions(list_junction_table_array.tree[[11]], reconstructed_gtf_table = reconstructed_gtf, mode = NULL, transcript.upstream_window = 50, transcript.downstream_window = 50)
  # 
  # list_three_frame_translate_result_chr12 <- purrr::map(.x = list_junction_table_array.tree, .f = ~three_frame_translate_splicejunctions(.x, reconstructed_gtf_table = reconstructed_gtf, mode = NULL, transcript.upstream_window = 50, transcript.downstream_window = 50))
  # 
  # test <- list_three_frame_translate_result_chr12 %>% compact %>% rbindlist
  
  ###################################################
  
  # use FASTA headers as the top level list names (for manual navigation and debug purposes)
  # complex filtering of virtual peptides using trim_virtual_peptides()
  cat("complex filtering of virtual peptides using trim_virtual_peptides()\n")
  # b <- purrr::map_depth(.x = list_three_frame_translate_result, .depth = 2, .ragged = FALSE, .f = ~purrr::modify_at(.x = .x, .at = "threeframetranslate", .f = ~purrr::map(.x = .x, .f = ~.x %>% trim_virtual_peptides)))
  b <- future_map(.x = list_three_frame_translate_result, 
                  .f = ~purrr::map(.x = .x, 
                                   .f = ~purrr::modify_at(.x = .x, .at = "threeframetranslate", 
                                                          .f = ~purrr::map(.x = .x, .f = ~.x %>% trim_virtual_peptides)
                                                          )
                                   ), .progress = TRUE, .options = future_options(globals = c("list_three_frame_translate_result", "trim_virtual_peptides", "upstream_window_size", "downstream_window_size", "last", "purrr::map", "purrr::modify_at"))
                  )
  # trim off resulting virtual AAs which are "none"
  cat("trim off resulting virtual AAs which are \"none\"\n")
  # c <- purrr::map_depth(.x = b, .depth = 2, .ragged = FALSE, .f = ~purrr::modify_at(.x = .x, .at = "threeframetranslate", .f = ~purrr::discard(.x, .p = ~.x %>% paste == "none")))
  c <- future_map(.x = b, 
                  .f = ~purrr::map(.x = .x, 
                                   .f = ~purrr::modify_at(.x = .x, .at = "threeframetranslate", 
                                                          .f = ~purrr::discard(.x, .p = ~.x %>% paste == "none")
                                                          )
                                   ), .progress = TRUE, .options = future_options(globals = c("b", "purrr::map", "purrr::modify_at", "purrr::discard"))
                  )
  # remove transcript-level elements that do not have any more extant three frame translate results 
  cat("remove transcript-level elements that do not have any more extant three frame translate results\n")
  d <- purrr::map(.x = c, .f = ~purrr::keep(.x = .x, .p = ~length(.x[["threeframetranslate"]]) >= 1))
  # remove junction entries that do not have any extant transcript-level results after filtering
  cat("remove junction entries that do not have any extant transcript-level results after filtering\n")
  e <- purrr::keep(.x = d, .p = ~length(.x) >= 1)
  
  # SANITY CHECK
  cat("CHECK: If the junction table used in translation was prefiltered by the same GTFs, then there must be no junctions that could not be translated:\n")
  cat("Number of junctions to be translated:", tibble_junction_table %>% nrow, "\n")
  
  list_translated_unsuccessfully <- purrr::map_depth(.x = list_three_frame_translate_result, .depth = 2, .ragged = FALSE, .f = ~any(length(.x[["threeframetranslate"]]) == 0) == TRUE) %>% purrr::map(~any(.x) == TRUE) %>% unlist
  
  cat("Number of junctions which could not be translated:", which(list_translated_unsuccessfully == TRUE) %>% length, "\n")
  
  # message of shame
  if (which(list_translated_unsuccessfully == TRUE) %>% length != 0) {
    
    cat("Uh-oh. There were some junctions which weren't translated :(
          
Fasta headers names of those which could not be translated:\n")
    
    indices_unsuccessful_junctions <- which(list_translated_unsuccessfully == TRUE)
    
    unsuccessful_junction_names <- purrr::map(.x = list_three_frame_translate_result[indices_unsuccessful_junctions], .f = ~.x[[1]][["fasta_header"]]) %>% unlist
    
    cat(unsuccessful_junction_names, sep = "\n")
    
    cat("- If you used the same GTF files for matching junctions (if you use something like stringtiemerge, it should still give you the same GTF as if you considered them separately), then consider whether all the chromosome names match. Inside the GTF file, must be in the format \"1\" not \"chr1\". 
- If you selected the option to translate all chromosomes incl. haplotype, fusion chromosomes, consider if the names in the junction table, GTF and reference fa. file all match.
- If you didn't filter using the same GTF, then this message may not mean there a problem. It just means that these junctions simply did not find any flanking exons.\n")
    
  }
  
  cat("Additionally, the number of junctions which didn't have any valid sequences due to filtering criteria:", ((list_three_frame_translate_result %>% length) - (e %>% length)), "\n")
  
  # COLLAPSING TIME
  # first collapse the 3FT virtual peptides into its own tibble
  cat("first collapse the 3FT virtual peptides into its own tibble\n")
  # f <- purrr::map_depth(.x = e, .depth = 2, .ragged = FALSE, .f = ~purrr::modify_at(.x = .x, .at = "threeframetranslate", .f = ~.x %>% as_tibble %>% t %>% as_tibble(rownames = "translation_frame") %>% setNames(c("translation_frame", "virtual_peptide_sequence"))) %>% flatten)
  
  f <- future_map(.x = e, .f = ~purrr::map(.x = .x, .f = ~purrr::modify_at(.x = .x, .at = "threeframetranslate", .f = ~.x %>% as_tibble %>% t %>% as_tibble(rownames = "translation_frame") %>% setNames(c("translation_frame", "virtual_peptide_sequence"))) %>% flatten), .progress = TRUE, .options = future_options(globals = c("e", "as_tibble", "purrr::map", "purrr::modify_at", "flatten")))

#   collapse3FTfunction <- function(junction_of_transcripts) {
# 
#     if (length(junction_of_transcripts) < 500) {
# 
#       purrr::map(.x = junction_of_transcripts,
#                  .f = ~purrr::modify_at(.x = .x, .at = "threeframetranslate",
#                                         .f = ~.x %>% as_tibble %>% t %>% as_tibble(rownames = "translation_frame") %>% setNames(c("translation_frame", "virtual_peptide_sequence"))
#                                         ) %>% flatten
#                  )
# 
#     } else {
# 
#       future_map(.x = junction_of_transcripts,
#                  .f = ~purrr::modify_at(.x = .x, .at = "threeframetranslate",
#                                         .f = ~.x %>% as_tibble %>% t %>% as_tibble(rownames = "translation_frame") %>% setNames(c("translation_frame", "virtual_peptide_sequence"))
# 
#                                         ) %>% flatten, .progress = TRUE
#                  )
# 
#     }
# 
#   }
# 
#   f <- purrr::map(.x = e,
#                   .f = ~collapse3FTfunction(.x)
#                   )
#   
  # collapse inside every junction
  cat("collapse inside every junction\n")
  g <- future_map(.x = f, .f = ~purrr::map(.x = .x, .f = ~.x %>% as_tibble) %>% rbindlist, .progress = TRUE, .options = future_options(globals = c("f", "as_tibble", "purrr::map", "rbindlist")))
  
  # ALL COLLAPSE, remove duplicated peptide sequences + header
  cat("ALL COLLAPSE, remove duplicated peptide sequences + header\n")
  h <- g %>% rbindlist %>% as_tibble %>% dplyr::distinct(., virtual_peptide_sequence, fasta_header, .keep_all = TRUE)
  # tidy up the translation frame column
  cat("tidy up the translation frame column\n")
  h[, "translation_frame"] <- gsub(x = h$translation_frame, pattern = "translation_frame_", replacement = "")
  
  # FINALLY! WE WRITE THE FASTA!
  # since our window was 100 nt, that means max. AA length of peptide = 33
  write.fasta(sequences = h$virtual_peptide_sequence %>% array_tree %>% flatten, names = h$fasta_header, file.out = paste(output_dir, output_database_name, ".fasta", sep = ""), open = "w", nbchar = 40, as.string = TRUE)
  
  # ADDITIONALLY! WE WRITE A TIBBLE CONTAINING THE GENOME COORD-PEPTIDE MAPPING
  j <- dplyr::left_join(h, tibble_junction_table, by = "fasta_header")
  # write a table
  write.table(x = j, file = paste(output_dir, output_database_name, "_3FT.summary.info.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
}
  
# finish counting
tictoc::toc()
  
q()
