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

matlock_bam <- snakemake@input[["bam"]]

########
# MAIN #
########

matlock_interactions <- fread(matlock_bam)
##make matrix of interactions
matlock_interaction_counts <- matlock_interactions[,.N,by=.(V2, V6)]
##write interaction matrix to analyse at later date
fwrite(matlock_interaction_counts, snakemake@output[["interaction_matrix"]])

interaction_locations <- matlock_interactions[,c(2,3,6,7)]
fwrite(interaction_locations, snakemake@output[["interaction_locations"]])
##may need to filter out viral interactions

# write log
sessionInfo()
