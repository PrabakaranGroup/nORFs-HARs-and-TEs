# Code to read in data from featureCounts and combine them into a single counts file

getwd()
# workingDir should contain feature counts files
workingDir = "~/"
setwd(workingDir)
library(dplyr)
library(tidyr)
library(textshape)
# Load counts files
all_counts<-read.table("outfeatureCounts.txt",check.names = FALSE,header = TRUE)
# Data frame with TE coordinates
TEcoord<-all_counts[,c(1:6)]
all_counts<-all_counts[rowSums(all_counts[,-c(1:6)]) > ncol(all_counts[,-c(1:6)]),]
# Modify sample names to remove filetype
new_colnames<-sub('.Aligned.sortedByCoord.out.bam','',colnames(all_counts))
colnames(all_counts)<-new_colnames
# Prepare featureCounts raw count (element-wise) for correlation
TEraw<-all_counts[,-c(2:6)]
TEraw<-column_to_rownames(TEraw)
TEraw<-TEraw[,colnames(tes)]
save(TEraw,TEcoord,file = "featureCounts.RData")
save(TEcoord,file = "TEcoord.RData")
