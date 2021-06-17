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
library(dplyr)

###########
# GLOBALS #
###########

matlock_bam <- snakemake@input[["bam"]]
viral_contigs_file <- snakemake@input[["viral_contigs_file"]]

########
# MAIN #
########

viral_contigs_table <- fread(viral_contigs_file)
viral_contigs <- subset(viral_contigs_table, plot_label == "Viral contig")

matlock_interactions <- fread(matlock_bam)
##make matrix of interactions
matlock_interaction_counts <- matlock_interactions[,.N,by=.(V2, V6)]
##write interaction matrix to analyse at later date
fwrite(matlock_interaction_counts, snakemake@output[["interaction_matrix"]])

##self interactions
self_interactions <- subset(matlock_interaction_counts, matlock_interaction_counts$V2 == matlock_interaction_counts$V6)
self_interactions$scaffold_no <- tstrsplit(self_interactions$V2, "_", keep=c(2))
self_interactions$scaffold_no <- as.character(self_interactions$scaffold_no)
self_interactions$scaffold_no <- as.numeric(self_interactions$scaffold_no)
setorder(self_interactions, scaffold_no)
##only hi-c self interactions
scaffolds <- self_interactions[c(1,2,3,4,5,6,7,8,9,10,11,12),]
sum(scaffolds$N)

##filter interaction matrix for viral contigs
V2_viral <- subset(matlock_interaction_counts, V2 %in% viral_contigs$`#Name`)
V6_viral <- subset(matlock_interaction_counts, V6 %in% viral_contigs$`#Name`)
viral_interactions <- full_join(V2_viral, V6_viral)
fwrite(viral_interactions, snakemake@output[['viral_interactions']])

##filter for viral contig interactions with hi-c scaffolds
hic_scaffolds <- subset(viral_contigs_table, plot_label == "Hi-C scaffold")
hic_viral_scaffolds <- subset(viral_contigs_table, plot_label == "Hi-C scaffold and viral")
all_hic<-full_join(hic_scaffolds, hic_viral_scaffolds)
viral_hic_V2 <- subset(viral_interactions, V2 %in% all_hic$`#Name`)
viral_hic_V6 <- subset(viral_interactions, V6 %in% all_hic$`#Name`)
viral_hic_interactions <- full_join(viral_hic_V2, viral_hic_V6)
sum(viral_hic_interactions$N)

interaction_locations <- matlock_interactions[,c(2,3,6,7)]
fwrite(interaction_locations, snakemake@output[["interaction_locations"]])
##may need to filter out viral interactions

# write log
sessionInfo()
