library(data.table)
library(dplyr)

##read in filtered matlock output - only lines that included the viral scaffolds
viral_matlock <- fread("output/mh_hic_matlock/matlock_viral_scaffolds.csv")
##list of viral scaffolds
viral_scaffold_ids <- fread("data/mh_genome/viral_scaffold_ids.txt", header=FALSE)

##get full interation table (each scaffold and counts of every scaffold it interacts with)
#group rows by V2 and V6 and then counts the number of rows in each group
viral_matlock_interaction_counts <- viral_matlock[, .N, by = .(V2,V6)]
##interaction matrix?
x <- dcast(viral_matlock_interaction_counts)

##filter for lines where scaffolds do not interact with themselves
viral_matlock_nonself <- viral_matlock[V2 != V6]
##self interactions - have multiple hits for each scaffold due to different positions interacting
viral_matlock_self <- viral_matlock[V2 == V6]


##how to see whether what they interact with is on list of viral scaffolds
viral_matlock_nonself$V2_viral_status <- ifelse(viral_matlock_nonself$V2 %in% viral_scaffold_ids$V1, 'viral', 'not viral')
viral_matlock_nonself$V6_viral_status <- ifelse(viral_matlock_nonself$V6 %in% viral_scaffold_ids$V1, 'viral', 'not viral')

##number where V2 or V6 is non-viral
##6775
sum(viral_matlock_nonself$V2_viral_status == 'not viral')
##7313
sum(viral_matlock_nonself$V6_viral_status == 'not viral')
##14,088 of 17,680 non-self interactions include one non-viral scaffold

##3592 interactions between 2 viral scaffolds
sum(viral_matlock_nonself$V2_viral_status == 'viral' & viral_matlock_nonself$V6_viral_status == 'viral')

##how many unique non-viral scaffolds have interactions with them
non_viral_interactions <- viral_matlock_nonself[viral_matlock_nonself$V2_viral_status == "not viral" | viral_matlock_nonself$V6_viral_status == "not viral"]

##alter table to be viral scaffold vs non-viral scaffold interactions??

##list of hi-c genome clusters by scaffold id
hic_clusters <- fread("data/mh_hic_genome/clusters.by_name_simple.txt", skip=16, fill=TRUE)
scaffold_to_cluster <- hic_clusters[,c(2,4)]
setnames(scaffold_to_cluster, old=c("V2", "V4"), new=c("scaffold_id", "cluster_id"))
scaffold_to_cluster$cluster_id <- tstrsplit(scaffold_to_cluster$cluster_id, "\\(", keep=c(2))
scaffold_to_cluster$cluster_id <- tstrsplit(scaffold_to_cluster$cluster_id, "\\)", keep=c(1))

##merge table with scaffold to cluster table to see if non-viral
##interacting scaffolds have been clustered into hic
v2_cluster_id <- merge(non_viral_interactions, scaffold_to_cluster, by.x="V2", by.y="scaffold_id", all.x=TRUE)
setnames(v2_cluster_id, old=c("cluster_id"), new=c("V2_cluster_id"))
nv_interactions_hic_status <- merge(v2_cluster_id, scaffold_to_cluster, by.x="V6", by.y="scaffold_id", all.x=TRUE)
setnames(nv_interactions_hic_status, old=c("cluster_id"), new=c("V6_cluster_id"))
nv_interactions_hic_status <- nv_interactions_hic_status[,c(2,4,16,18,1,7,17,19)]


##get table to viral scaffold vs interacting scaffold
v2_nv <- subset(nv_interactions_hic_status, nv_interactions_hic_status$V2_viral_status == "not viral")
##reorder table so that viral scaffold info first
v2_nv <- v2_nv[,c(5,6,8,1,2,4)]
setnames(v2_nv, old=c("V6", "V7", "V6_cluster_id", "V2", "V3", "V2_cluster_id"), new=c("viral_scaffold", "viral_scaffold_position", "viral_scaffold_cluster_id", "nv_interacting_scaffold", "nv_interacting_scaffold_position", "nv_interacting_scaffold_cluster_id"))
##get other half of table to merge with
v2_v <- subset(nv_interactions_hic_status, nv_interactions_hic_status$V2_viral_status == "viral")
v2_v <- v2_v[,c(1,2,4,5,6,8)]
setnames(v2_v, old=c("V2", "V3", "V2_cluster_id", "V6", "V7", "V6_cluster_id"), new=c("viral_scaffold", "viral_scaffold_position", "viral_scaffold_cluster_id", "nv_interacting_scaffold", "nv_interacting_scaffold_position", "nv_interacting_scaffold_cluster_id"))
##join two tables now that viral scaffolds all in first column, with non-viral interacting all in column 4
full_nv_interaction_table <- full_join(v2_nv, v2_v)
fwrite(full_nv_interaction_table, "output/mh_hic_matlock/viral_nonviral_scaffold_interactions.csv")

##how many times does each viral scaffold interact with something?
##all 21 viral scaffolds interact
viral_interaction_counts <- count(full_nv_interaction_table, full_nv_interaction_table$viral_scaffold)
##scaffold3939 was another I suspected might not be viral - it has 4809 interactions
##scaff28315 next highest - also has a non-viral gene on it
##both had hits to Oryctes rhinoceros nudivirus
only_viral_interaction_counts <- viral_interaction_counts[viral_interaction_counts$`full_nv_interaction_table$viral_scaffold`!="Scaffold28315",]
only_viral_interaction_counts <- only_viral_interaction_counts[only_viral_interaction_counts$`full_nv_interaction_table$viral_scaffold`!="Scaffold3939",]
##mean no interactions for scaffolds with only viral genes
mean(only_viral_interaction_counts$n)


##sum of nv interacting partners where they ARE clustered into hi-c - 7654 (54.3%)
sum(!is.na(full_nv_interaction_table$nv_interacting_scaffold_cluster_id))
##sum of nv partners where they are NOT clustered into hi-c - 6434 (45.7%)
sum(is.na(full_nv_interaction_table$nv_interacting_scaffold_cluster_id))


##no. interactions per viral scaffold
##compare to mean no. interactions for non viral scaffolds?
##what hic clusters is each scaffold interacting with and how many times for each?
##- looks like they all interact with all clusters so far?

