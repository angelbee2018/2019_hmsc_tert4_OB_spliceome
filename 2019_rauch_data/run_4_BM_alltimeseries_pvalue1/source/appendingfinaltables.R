library(dplyr)


for (eventname in c("A3SS_events", "A5SS_events", "cassette_exon_events", "composite_events", "intron_retention", "MXE_events"))
{
temp_oneevent_table <- read.delim(file = paste("/srv/scratch/z3463471/2019_rauch_data/bam_bai_sj_sam_files/run_1/JUM_diff/FINAL_JUM_OUTPUT_pvalue_0.05/AS_differential_JUM_output_", eventname, "_pvalue_0.05_final_detailed.txt", sep = ""), header = TRUE, sep = "\t")
temp_oneevent_table[, "event_name"] <- paste(eventname)
write.table(x = temp_oneevent_table, file = paste("/srv/scratch/z3463471/2019_rauch_data/bam_bai_sj_sam_files/run_1/JUM_diff/FINAL_JUM_OUTPUT_pvalue_0.05/", eventname, "_appended.txt", sep = ""), sep = "\t")
}





#temp_alleventstable <- data.frame(matrix(ncol = 30, nrow = 0))
#temp_alleventstable <- dplyr::bind_rows(temp_alleventstable, temp_oneevent_table)
#for (colnumber in c(5:15, 17:22, 29))
#{
#temp_oneevent_table[, colnumber] <- as.numeric(temp_oneevent_table[, colnumber])
#}
