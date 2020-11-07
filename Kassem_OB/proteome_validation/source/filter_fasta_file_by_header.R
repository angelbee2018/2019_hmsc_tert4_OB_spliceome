

fasta_raw <- seqinr::read.fasta(file = "/mnt/Helium_8TB_1/uniprot_trembl.fasta", seqtype = "DNA", forceDNAtolower = FALSE, whole.header = TRUE)
fasta_human <- fasta_raw[grep(x = fasta_raw %>% names, pattern = "homo sapiens", ignore.case = TRUE)]

seqinr::write.fasta(sequences = fasta_human, names = names(fasta_human), file.out = "/mnt/Tertiary/sharedfolder/hg38_ensembl_reference/proteome_fasta/uniprot_sprot_isoforms_human.fasta")
