# THREE-FRAME TRANSLATION ######
# Will now attempt to translate the nucleotides of all transcripts flanking all the splice junctions given.
# requires at least 60 GB of RAM for ~2500 junctions. 80 GB to be safe. Takes around 8 hours for 2500 junctions.
# 4 cores minimum. 6 recommended.

# arg1: path to a table containing the junction and reconstructed GTF file information. should be in the format of a tab separated tibble to minimise the risk of error. each row describes A PATH the file to the junction table and reconstr. GTF to be translated. This program will automatically loop thru them. Set col names as: c("reconstructed_GTF_path", "junction_table_path")
# junction table: path to table containing junctions of interest. for example, UNION_junc_coor type files. MUST contain columns start, end, chr, strand. one row per junction. should be in the form of a tibble.
# reconstructed GTF table: path to the actual reconstructed GTF file (e.g. from Cufflinks, Strawberry). NOT THE CONTAINING DIRECTORY. tip: for better junction matching, combine the reconstructed GTF with reference GTF beforehand e.g. using StringTie

# arg2: path to the actual genome FASTA file. Ideally from Ensembl... you need the whole primary assembly. Not separate files by chromosome.
# arg3: how many nucleotides to translate upstream W.R.T. the transcript. for junction-based mode, the total nt. length to be translated will be arg4 + arg5
# arg4: how many nt to translate downstream W.R.T. the transcript. for exon-based mode, this will be 0 by default. the total arg4 + arg5 will indicate how many nt to translate beyond the exon.
# arg5: output directory
# arg6: output_database_name: a string of what the final FASTA database file name and output table will be called

print(commandArgs())

args = commandArgs(trailingOnly=TRUE)

print(args)

print(paste("number of arguments: ", length(args), sep = ""))

# test if there are 6 arguments: if not, return an error
if (length(args) != 6) {
  stop("check the amount of arguments you have. should have specified 6", call.=FALSE)
}

# SET ENVIRONMENT ##########

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("seqinr", "tidyverse", "purrr", "dplyr", "rtracklayer"))

library(seqinr)
library(tidyverse)
library(purrr)
library(dplyr)
library(rtracklayer)

file_information_path <- args[1]

reference_genome_fasta_path <- args[2]
output_dir <- args[5]

# testing ########

# UNION_junc_coor_differential_path <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/input_files/"
# reconstructed_gtf_path <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/GRAND_OBseries_ref_reconstructed_stringtiemerged.gtf"
# reference_genome_fasta_path <- "/srv/scratch/z3463471/hg38_ensembl_reference/raw_genome_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# output_dir <- "/srv/scratch/z3463471/PGNEXUS_kassem_MSC/Kassem_OB/proteome_validation/results_database_generation/"
# output_database_name <- "three_frame_translation_angel_table.txt"

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
  
  return(list("vector_genome.coords_to_translate" = vector_nt.positions_to_translate, "chr" = list_junction.adjacent_exon_start.end_coords_one_transcript$chr, "strand" = list_junction.adjacent_exon_start.end_coords_one_transcript$strand))
  
}

# END generate_all.genomic.positions_for_translation

# FUNCTION TO 3 FRAME TRANSLATE ONE LIST CONTAINING NUCLEOTIDE SEQUENCE AND STRAND

nt.sequence_strand_threeframetranslate <- function(list_nt.fwd.sequence_strand) {
  
  if (list_nt.fwd.sequence_strand$strand == "+") {
    
    translation_result <- list("translation_frame_0" = translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 0, sens = "F"),
                               "translation_frame_1" = translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 1, sens = "F"),
                               "translation_frame_2" = translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 2, sens = "F"))
    
  } else if (list_nt.fwd.sequence_strand$strand == "-") {
    
    translation_result <- list("translation_frame_0" = translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 0, sens = "R"),
                               "translation_frame_1" = translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 1, sens = "R"),
                               "translation_frame_2" = translate(list_nt.fwd.sequence_strand$forward_nucleotides, frame = 2, sens = "R"))
    
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

three_frame_translate_splicejunctions <- function(spliceregion_list, loop_number, reconstructed_gtf_table, mode = NULL, transcript.upstream_window = 50, transcript.downstream_window = 50) {
  
  # DEBUG ###################
  
  # transcript.upstream_window <- 50
  # transcript.downstream_window <- 50
  # 
  # reconstructed_gtf_table <- reconstructed_gtf
  # spliceregion_list <- UNION_junc_coor_differential_shiftedback_arraytree[[2629]]
  
  ###########################
  
  # filter the reconstructed GTF table for all exon entries that directly flank the splice junction
  
  print(paste("Now processing junction number:", loop_number))
  
  print("finding which transcripts contain the specified junction...")
  
  if (spliceregion_list$strand == ".") {
    
    tibble_reconstructed_gtf_subset_flanking_exons <- reconstructed_gtf_table[reconstructed_gtf_table$seqnames == spliceregion_list$chr, ] %>% .[.$start <= ((spliceregion_list$end %>% as.numeric) + 1) & .$end >= ((spliceregion_list$start %>% as.numeric) - 1), ] %>% .[!(.$start <= (spliceregion_list$end %>% as.numeric) & .$end >= (spliceregion_list$start %>% as.numeric)), ] %>% .[.$type == "exon", ]
    
  } else if (spliceregion_list$strand == "+" | spliceregion_list$strand == "-") {
    
    tibble_reconstructed_gtf_subset_flanking_exons <- reconstructed_gtf_table[reconstructed_gtf_table$seqnames == spliceregion_list$chr, ] %>% .[.$strand == spliceregion_list$strand, ] %>% .[.$start <= ((spliceregion_list$end %>% as.numeric) + 1) & .$end >= ((spliceregion_list$start %>% as.numeric) - 1), ] %>% .[!(.$start <= (spliceregion_list$end %>% as.numeric) & .$end >= (spliceregion_list$start %>% as.numeric)), ] %>% .[.$type == "exon", ]
    
  } else {
    
    stop("Could not match the strand information in the transposed differential-only UNION_junc_coor_table. Make sure that the \"strand\" column in the UNION_junc_coor_table contains only +, - or .")
    
  }
  
  list_of_junction_associated_transcripts <- tibble_reconstructed_gtf_subset_flanking_exons$transcript_id %>% unique %>% array_tree %>% flatten
  
  # make a list for each transcript that directly flanks a junction.
  # then filter so that there are only a) exon PAIRS which b) are directly connected in the mature (spliced) transcript
  
  list_of_tibbles_flanking_exon_gtf.entries_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~tibble_reconstructed_gtf_subset_flanking_exons[tibble_reconstructed_gtf_subset_flanking_exons$transcript_id == .x, ] %>% dplyr::arrange(exon_number %>% as.numeric)) %>% set_names(list_of_junction_associated_transcripts) %>% keep(.x = ., .p = ~nrow(.x) == 2) %>% keep(.x = ., .p = ~abs((.x[2, "exon_number"] %>% paste %>% as.numeric) - (.x[1, "exon_number"] %>% paste %>% as.numeric)) == 1)
  
  #### these are the two main tables to be used for translation #########
  
  # to account for + or - strand transcription, the positive "sign" indicates the second exon has higher genome coordinate than the first. negative "sign" indicates the opposite.
  
  print("defining junction-flanking nucleotide positions...")
  
  list_of_junction.adjacent_exon_start.end_coords_per_transcript <- purrr::map(.x = list_of_tibbles_flanking_exon_gtf.entries_per_transcript, .f = ~list("end_transcript.5prime_exon_coord" = .x[which(.x$exon_number %>% as.numeric == .x$exon_number %>% as.numeric %>% min), "end"] %>% paste %>% as.numeric, "start_transcript.3prime_exon_coord" = .x[which(.x$exon_number %>% as.numeric == .x$exon_number %>% as.numeric %>% max), "start"] %>% paste %>% as.numeric, "chr" = .x[1, "seqnames"] %>% paste, "strand" = .x[1, "strand"] %>% paste) %>% splice("sign" = (.$start_transcript.3prime_exon_coord - .$end_transcript.5prime_exon_coord)/abs(.$start_transcript.3prime_exon_coord - .$end_transcript.5prime_exon_coord)))
  
  print("define genome-relative co-ordinates for each nucleotide in associated transcripts...")
  
  list_all.genomic.positions_per_transcript <- purrr::map(.x = list_of_junction_associated_transcripts, .f = ~reconstructed_gtf_table[reconstructed_gtf_table$transcript_id == .x & reconstructed_gtf$type == "exon", c("start", "end")] %>% array_tree %>% purrr::map(.x  = ., .f = ~.x[[1]]:.x[[2]]) %>% unlist %>% unique %>% sort) %>% set_names(list_of_junction_associated_transcripts) %>% .[names(list_of_tibbles_flanking_exon_gtf.entries_per_transcript)]
  
  #######################################################################
  
  # generate genome-relative coordinates flanking each junction for translation
  
  print("define genomic positions of nucleotides to be translated..")
  
  list_genome.coords_for_translation <- purrr::map2(.x = list_of_junction.adjacent_exon_start.end_coords_per_transcript, .y = list_all.genomic.positions_per_transcript, .f = ~generate_all.genomic.positions_for_translation(.x, .y, transcript.upstream_window = transcript.upstream_window, transcript.downstream_window = transcript.downstream_window))
  
  print("looking up reference FASTA to fetch nucleotides within specified transcript window...")
  
  list_forward_nucleotides_from_coords <- purrr::map(.x = list_genome.coords_for_translation, .f = ~list("forward_nucleotides" = reference_genome_fasta[[.x$chr]] %>% .[.x$vector_genome.coords_to_translate], "strand" = .x$strand))
  
  print("initiating three-frame translation...")
  
  # three-frame translation
  
  list_threeframetranslate <- purrr::map(.x = list_forward_nucleotides_from_coords, .f = ~nt.sequence_strand_threeframetranslate(.x))
  
  print("...Success!!")
  
  return(list_threeframetranslate)
  
}

# END three_frame_translate_splicejunctions

# BEGIN EXECUTION #################################

reference_genome_fasta <- seqinr::read.fasta(file = reference_genome_fasta_path, forceDNAtolower = FALSE)

file_information_table <- read_delim(reference_genome_fasta_path, delim = "\t")

## loop thru reconstructed GTF path slowly, loop thru junction files quickly.

for (reconstructed_gtf_path in file_information_table$reconstructed_GTF_path) {
  
  reconstructed_gtf <- rtracklayer::import(reconstructed_gtf_path) %>% as_tibble %>% dplyr::mutate_if(is.factor, as.character)
  
}
     
     
     
     c("reconstructed_GTF_path", "junction_table_path")




UNION_junc_coor_differential <- read_delim(UNION_junc_coor_differential_path, delim = "\t")

UNION_junc_coor_differential_path <- args[1]
reconstructed_gtf_path <- args[2]

###################################################

UNION_junc_coor_differential_arraytree <- UNION_junc_coor_differential %>% array_tree

# test <- three_frame_translate_splicejunctions(UNION_junc_coor_differential_shiftedback_arraytree[[11]], reconstructed_gtf_table = reconstructed_gtf, mode = NULL, transcript.upstream_window = 50, transcript.downstream_window = 50)
# 
# list_three_frame_translate_result_chr12 <- purrr::map(.x = UNION_junc_coor_differential_shiftedback_arraytree, .f = ~three_frame_translate_splicejunctions(.x, reconstructed_gtf_table = reconstructed_gtf, mode = NULL, transcript.upstream_window = 50, transcript.downstream_window = 50))
# 
# test <- list_three_frame_translate_result_chr12 %>% compact %>% rbindlist

list_three_frame_translate_result <- purrr::imap(.x = UNION_junc_coor_differential_arraytree, .f = ~three_frame_translate_splicejunctions(.x, .y, reconstructed_gtf_table = reconstructed_gtf, mode = NULL, transcript.upstream_window = 50, transcript.downstream_window = 50))

list_three_frame_translate_result %>% length

list_three_frame_translate_result %>% compact %>% length

# use FASTA headers as the top level list names

a <- list_three_frame_translate_result %>% set_names(UNION_junc_coor_differential$fasta_header)

b <- purrr::map_depth(.x = a, .depth = 3, .ragged = FALSE, .f = ~.x %>% unlist %>% paste(collapse = ""))

c <- purrr::map_depth(.x = b, .depth = 3, .ragged = FALSE, .f = ~.x %>% paste %>% strsplit(., split = "\\*") %>% .[[1]] %>% .[1])

d <- purrr::map_depth(.x = c, .depth = 2, .ragged = FALSE, .f = ~purrr::discard(.x, .p = ~.x %>% paste %>% nchar < 6))

e <- purrr::map_depth(.x = d, .depth = 2, .ragged = FALSE, .f = ~.x %>% unlist %>% setNames(., NULL))

f <- purrr::map(.x = e, .f = ~.x %>% unlist %>% setNames(., NULL) %>% unique) %>% compact

g <- stack(f) %>% as_tibble %>% setNames(c("peptide_sequence", "fasta_header"))

# FINALLY! WE SAVE THE TIBBLE!

write.table(x = g, file = paste(output_dir, output_database_name, ".txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# FINALLY! WE WRITE THE FASTA!

# since our window was 100 nt, that means max. AA length of peptide = 33

write.fasta(sequences = g$peptide_sequence %>% array_tree %>% flatten, names = g$fasta_header, file.out = paste(output_dir, output_database_name, ".fasta", sep = ""), open = "w", nbchar = 40, as.string = TRUE)

q()