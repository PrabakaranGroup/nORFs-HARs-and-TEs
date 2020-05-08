getwd()
workingDir = "~/"
setwd(workingDir)
library(WGCNA)
library(stringr)
library(edgeR)
library(dplyr)
options(stringsAsFactors = FALSE)
tbn0 <- read.table("2.txt",check.names = FALSE,header = TRUE,row.names = NULL)
nORFtranscriptids<- read.table("norf_transcripts_hg19_matches.txt",header = TRUE,check.names = FALSE)
tstochange<-tbn0$row.names%in%nORFtranscriptids$tsc_id
tstochange1<-match(tbn0$row.names[tstochange],nORFtranscriptids$tsc_id)
tbn0$row.names[tstochange]<-nORFtranscriptids$norf_id[tstochange1]
hhh<-duplicated(pull(tbn0,row.names))|duplicated(pull(tbn0,row.names),fromLast=TRUE)
hhhh<-tbn0[hhh,]
tbn<- tbn0 %>% group_by(row.names) %>% summarise_all(list(sum))
# Select BrainGVEX samples only
tbn1 = as.data.frame(tbn[,2:428])
row.names(tbn1)<- tbn$row.names
sam = names(tbn1)
traitData = read.csv("ClinicalData0.csv")
relevantSamples1 = match(sam,traitData$individualID)
datTraits1 = traitData[relevantSamples1,-1]
rownames(datTraits1) = sam
# normalise counts
tra<- factor(datTraits1$diagnosis)
tra<-relevel(tra,ref = "Bipolar Disorder")
p<-DGEList(counts = tbn1,group = tra)
keep<-filterByExpr(p)
p<-p[keep, , keep.lib.sizes=FALSE]
p<-calcNormFactors(p)
q<-as.data.frame(log2(cpm(p)+1))
# Isolate nORF counts
nORFdata<-q[!with(q,grepl("ENS",row.names(q))),]
datExpr1 = as.data.frame(t(nORFdata))
gsg1 = goodSamplesGenes(datExpr1, verbose = 3)
gsg1$allOK
if (!gsg1$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg1$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr1)[!gsg1$goodGenes], collapse = ", ")))
  if (sum(!gsg1$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg1$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg1$goodSamples, gsg1$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 37, col = "blue")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 37, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr1[keepSamples,]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
traitData$Suicide<-0
traitData$Suicide[grep("sui|SUIC|Suic",traitData$causeDeath)]<-1
relevantSamples = match(rownames(datExpr),traitData$individualID)
dim(traitData)
names(traitData)
allTraits = traitData[, -c(5,6,9,11,15,16)]
allTraits = allTraits[,c(4:9,11:13,16)]
datTraits = allTraits[c(unique(relevantSamples)),-1]
rownames(datTraits) = allTraits[unique(relevantSamples),1]

collectGarbage()
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# Choose a set of soft-thresholding powers
powers = c(c(1:10),seq(from = 12,to=20,by=2))
# Call netwrok topology analysis function
sft = pickSoftThreshold(datExpr,powerVector = powers, verbose = 5)
#Plot results
sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free toppology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers,cex = cex1,col = "red")
# R^2 cutoff
abline(h=0.90, col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels = powers,cex = cex1,col="red")
# Similarity and adjacency
softPower = 7
adjacency = adjacency(datExpr,power = softPower)
# Topological Overlap Matrix
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
# Clustering using TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12,9)
plot(geneTree,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity",
     labels=FALSE, hang = 0.04)
minModuleSize = 30
# Module identification
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colours")
# Merge similar modules
# Eigengenes
MEList = moduleEigengenes(datExpr,colors = dynamicColors)
MEs = MEList$eigengenes
# Eigengene dissimilarity
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot
sizeGrWindow(7,6)
plot(METree, main = "Clustering of module eigengenes",
     xlab="",sub="")
MEDissThres=0.25
abline(h=MEDissThres)
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# new colours
mergedColors = merge$colors
# new eigengenes
mergedMEs = merge$newMEs
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Rename moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

save(datExpr, datTraits, file = "01-dataInput.RData")
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Identify DE nORF transcripts and modules they belong to 
DE = read.table('de_norfs.txt',header = TRUE, row.names = 1,check.names = FALSE)
DEmatch = match(row.names(DE),row.names(nORFdata))
DEmatch
fullrecords <-  DEmatch[!complete.cases(DEmatch)]
droprecords <-  DEmatch[complete.cases(DEmatch)]
man = row.names(nORFdata)
man = as.data.frame(man[droprecords])
man[2] = moduleLabels[droprecords]
man[3] = moduleColors[droprecords]

# Identify nORFs near HARs
# HARs
allHARs <- read.table("Year 4/mergedHARs.txt",header = TRUE)
# nORFs
allnORFs <- read.table("norf_coordinates.bed")
allnORFs <- allnORFs %>% distinct()
Ihar2 <- GRanges(seqnames = Rle(allHARs$Chr,1),ranges = IRanges(allHARs$Start-100000,end = allHARs$End+100000,names = allHARs$ID),strand = NULL)
InORFs <- GRanges(seqnames = Rle(allnORFs$V1,1),ranges = IRanges(allnORFs$V2,end = allnORFs$V3,names = allnORFs$V4),strand = allnORFs$V6)
HARs_nORFs <- findOverlaps(Ihar2,InORFs)
nORFs_near<-allnORFs[unique(HARs_nORFs@to),]
nearmatch = match(nORFs_near$V4,row.names(nORFdata))
fullrecords1 <-  nearmatch[!complete.cases(nearmatch)]
droprecords1 <-  nearmatch[complete.cases(nearmatch)]
man1 = row.names(nORFdata)
man1 = as.data.frame(man1[droprecords1])
man1[2] = moduleLabels[droprecords1]
man1[3] = moduleColors[droprecords1]

# nORF not nears
man2 = row.names(nORFdata)
man2 = as.data.frame(man2[-c(droprecords1)])
man2[2] = moduleLabels[-c(droprecords1)]
man2[3] = moduleColors[-c(droprecords1)]

man4 = row.names(nORFdata)
man5 = data.frame(nORF = man4,module = moduleColors, nearHAR = 1, DE = 1, category = 1)
man5$nearHAR[droprecords1]="near HAR"
man5$nearHAR[-c(droprecords1)]="not near HAR"
man5$DE[droprecords]="DE"
man5$DE[-c(droprecords)]="non DE"
man5$category=paste(man5$DE,man5$nearHAR,sep = " ")

library(ggplot2)
man5_base <- ggplot(data = man5, aes(x=module,fill=category))
man5_base + geom_bar() + theme_classic(base_size = 13) + scale_y_continuous(expand = c(0,0))