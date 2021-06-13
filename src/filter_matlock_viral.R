#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)

###########
# GLOBALS #
###########

matlock_output <- snakemake@input[["matlock_output"]]
viral_scaffold_list <- snakemake@input[["viral_scaffold_list"]]

########
# MAIN #
########

matlock_table <- fread(matlock_output)
viral_scaffolds <-  fread(viral_scaffold_list, header=FALSE)

matlock_viral_scaffolds <- matlock_table[matlock_table$V2 %in% viral_scaffolds$V1 | matlock_table$V6 %in% viral_scaffolds$V1]
fwrite(matlock_viral_scaffolds, snakemake@output[["viral_matlock"]])

# write log
sessionInfo()
