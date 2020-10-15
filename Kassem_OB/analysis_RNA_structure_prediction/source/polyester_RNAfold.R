
library(tidyverse)
library(polyester)
library(rtracklayer)

# gtf_path <- "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/alltimepoints_denovo_reconstructed_stringtiemerged.gtf"
gtf_path <- "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/alltimepoints_denovo_reconstructed_stringtiemerged_chromosomal_and_MT_only.gtf"

raw_genome_fasta_dir <- "/mnt/Tertiary/sharedfolder/hg38_ensembl_reference/raw_genome_fasta/dna_by_chr_short_names/"

tibble_gtf <- rtracklayer::import(con = gtf_path, format = "gtf") %>% as_tibble

rtracklayer::export(tibble_gtf %>% dplyr::filter(seqnames %in% c(1:22, "X", "Y", "MT")), format = "gtf", con = "/mnt/Tertiary/sharedfolder/PGNEXUS_kassem_MSC/Kassem_OB/analysis_strawberry/results_assemblyonly/merged/alltimepoints_denovo_reconstructed_stringtiemerged_chromosomal_and_MT_only.gtf")

DNAStringSet_gtf_nt_sequences <- polyester::seq_gtf(gtf = gtf_path, 
                                                    seqs = raw_genome_fasta_dir, 
                                                    feature = "transcript", 
                                                    exononly = TRUE,
                                                    idfield = "transcript_id", 
                                                    attrsep = "; ")
