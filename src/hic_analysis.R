##filter for interactions near start or end of contigs
##clusters by name has length of contigs for those that clustered
##what about contigs that did not cluster?

##should still see accurate interactions across whole scaffolds not just at ends - so probably shouldn't filter

##dont have any kind of bias correction of hi-c data
##is low no interactions because we may have tissue specificity and so only detect interactions in small no. tissues?

##, interaction frequency between loci in cis decreases, on average, as their genomic distance increase

##then how to filter for interactions near ends of scaffolds?

#samtools idxstats in.sam|in.bam|in.cram

#Retrieve and print stats in the index file corresponding to the input file. Before calling idxstats, the input BAM file should be indexed by samtools index.
#If run on a SAM or CRAM file or an unindexed BAM file, this command will still produce the same summary statistics, but does so by reading through the entire file. This is far slower than using the BAM indices.
#The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads. It is written to stdout.