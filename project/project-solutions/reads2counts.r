#!/usr/bin/env Rscript
library("dada2")

#############################################################################################################
## This script implements an Amplicon Sequence Variant (ASV) approach to find unique (organism-specific)   ## 
##                        sequence variants of the 16S gene in each given sample.                          ##
##  The script takes forward & reverse reads (in no specific order) as input arguments on the command line ##
##                       and should be run only once over all samples in a dataset.                        ##
##           Depending on the available computational power, this script can take a while to run.          ##
#############################################################################################################

reads = commandArgs(trailingOnly=TRUE)

# find and extract the forward and reverse reads
forward_reads <- reads[grep("_R1_", reads)]
reverse_reads <- reads[grep("_R2_", reads)]

# learn the error rate of the base calling
err_forward_reads <- learnErrors(forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(reverse_reads, multithread=TRUE)

# find unique sequences
dada_forward <- dada(forward_reads, err=err_forward_reads, pool="pseudo", multithread=TRUE)
dada_reverse <- dada(reverse_reads, err=err_reverse_reads, pool="pseudo", multithread=TRUE)

# merge forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, forward_reads, dada_reverse,
                    reverse_reads, trimOverhang=TRUE, minOverlap=20)

# make a sequence abundance table
seqtab <- makeSequenceTable(merged_amplicons)
names <- sprintf("ASV_%s", seq(1:ncol(seqtab)))
colnames(seqtab) <- names

# write out this table to a csv file
write.table(seqtab, "counts_matrix.csv", sep=",", quote=F, col.names=NA)

# hierarchically cluster the samples on euclidean distance
euc_dist <- dist(seqtab)
euc_clust <- hclust(euc_dist, method="ward.D2")
euc_dend <- as.dendrogram(euc_clust, hang=0.1)

# plot and save
png("dendrogram.png")
plot(euc_dend, ylab="Euc. dist.")
dev.off() 