library(tidyverse)

tibble_sraruntable <- read.delim("/mnt/Helium_8TB_1/aarhus_smallseq_timeseries_mirna_BM_MSC_to_OB_GSE107279/data/SraRunTable.txt", sep = ",", stringsAsFactors = FALSE) %>% as_tibble

tibble_series_matrix <- tibble("Sample.Name" = "GSM2863648 GSM2863649 GSM2863650 GSM2863651 GSM2863652 GSM2863653 GSM2863654 GSM2863655 GSM2863656 GSM2863657 GSM2863658 GSM2863659 GSM2863660 GSM2863661 GSM2863662 GSM2863663 GSM2863664 GSM2863665 GSM2863666 GSM2863667 GSM2863668 GSM2863669 GSM2863670 GSM2863671 " %>% strsplit(split = " ") %>% unlist, 
                               "timepoint_replicate" = "OB_G1-12h OB_G1-24h OB_G1-6h OB_G1-d0 OB_G1-d10 OB_G1-d13 OB_G1-d3 OB_G1-d7 OB_G2-12h OB_G2-24h OB_G2-6h OB_G2-d0 OB_G2-d10 OB_G2-d13 OB_G2-d3 OB_G2-d7 OB_G3-12h OB_G3-24h OB_G3-6h OB_G3-d0 OB_G3-d10 OB_G3-d13 OB_G3-d3 OB_G3-d7" %>% strsplit(split = " ") %>% unlist)

# table join
tibble_sraruntable_with_timeseries_info <- dplyr::left_join(tibble_sraruntable, tibble_series_matrix, by = "Sample.Name")

# extract SRR to timepoint/replicate mapping
tibble_SRR_to_timeseries_info <- tibble_sraruntable_with_timeseries_info[, c("Run", "timepoint_replicate")] %>%
  # edit the timeseries info to our liking
  dplyr::mutate("final_fastqfilenames" = gsub(x = `timepoint_replicate`, pattern = "OB", replacement = "BM_MSC_to_OB") %>% 
                  gsub(pattern = "(.*)G([1-3])\\-(.*)", replacement = "\\1\\3\\_r\\2") %>%
                  gsub(pattern = "(.*_)d([^_]+)(.*)", replacement = "\\1\\2d\\3") %>%
                  gsub(pattern = "_0d_", replacement = "_MSC_") %>%
                  gsub(pattern = "24h", replacement = "1d")) 

# write the table
write.table(x = tibble_SRR_to_timeseries_info, file = "/mnt/Helium_8TB_1/aarhus_smallseq_timeseries_mirna_BM_MSC_to_OB_GSE107279/data/SRR_to_timeseries_info.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
