library('WGCNA')
library('stringr')
library('GenomicRanges')
library('dplyr')
library(edgeR)
options(stringsAsFactors = FALSE)
getwd()
# Set the working directory to a folder containing files in this code to be loaded
workingDir = "~/"
setwd(workingDir)
# Import TE GTF
repeats<- read.table('GRCh37_GENCODE_rmsk_TE.gtf')
# Load transcript counts and switch transcript ids for nORF ids
allTranscripts <- read.table("2.txt",check.names = FALSE,header = TRUE)
nORF_transcript_matches<- read.table("norf_transcripts_hg19_matches.txt",header = TRUE,check.names = FALSE)
row.names(allTranscripts)[row.names(allTranscripts)%in%nORF_transcript_matches$tsc_id]<-make.unique(paste(nORF_transcript_matches$norf_id))
# Check that transcript IDs have been switched with nORF IDs
nORFtranscriptmatches<-match(nORF_transcript_matches$norf_id,row.names(allTranscripts))
# Load in TE raw counts and coordinates, keep relevant samples of transcripts
load("featureCount.RData")
allTranscripts<-allTranscripts[,colnames(TEraw)]
allTranscripts<-allTranscripts[rowSums(allTranscripts) > 427,]
# Normalise the counts by TMM normalisation
transcriptome<-rbind(allTranscripts,TEraw)
traitData = read.csv("ClinicalData0.csv")
sampleNames<-names(TEraw)
relevantTraits<-subset(traitData,traitData$individualID%in%sampleNames)
TraitFactor<- factor(relevantTraits$diagnosis)
TraitFactor<-relevel(TraitFactor,ref = "Bipolar Disorder")
DGE<-DGEList(transcriptome, group = TraitFactor)
keep<-filterByExpr(DGE)
DGE<-DGE[keep, , keep.lib.sizes=FALSE]
DGE<-calcNormFactors(DGE)
norm_counts<-as.data.frame(t(cpm(DGE)))
# Split up TE and nORF expression
TEexpr0<-norm_counts[,with(norm_counts, grep("\\:",colnames(norm_counts)))]
nORFexpr0<-norm_counts[,with(norm_counts, grep("\\:|ENST",colnames(norm_counts),invert=TRUE))]
save(norm_counts,TEexpr0,nORFexpr0,file = "normalised-counts.RData")
# load("normalised-counts.RData")
TEcoord$Family<-repeats[match(str_c(str_split_fixed(TEcoord$Geneid,":",4)[,c(4)]),repeats$V10),"V16"]
TEcoord$Class<-repeats[match(str_c(str_split_fixed(TEcoord$Geneid,":",4)[,c(4)]),repeats$V10),"V19"]
TEGRange <- GRanges(seqnames = Rle(TEcoord$Chr,1),ranges = IRanges(start = TEcoord$Start ,end = TEcoord$End,names = TEcoord$Geneid),strand = TEcoord$Strand)
# Load in nORF coordinates and remove duplicates
allnORFs <- read.table("norf_coordinates.bed")
allnORFs <- allnORFs %>% distinct()
# GRanges object of the 2000 bp upstream of nORFs
nORFwindow <- GRanges(seqnames = Rle(allnORFs$V1,1),ranges = IRanges(allnORFs$V2-2000,end = allnORFs$V2,names = allnORFs$V4),strand = allnORFs$V6)
nORFwindow2  <- GRanges(seqnames = Rle(allnORFs$V1,1),ranges = IRanges(allnORFs$V3,end = allnORFs$V3+2000,names = allnORFs$V4),strand = allnORFs$V6)
TE_nORFoverlap <- findOverlaps(TEGRange,nORFwindow)
#Identify TE-associated nORFs and nORF-associated TEs and subset expression of these
nORFswithTEs <- as.character(allnORFs$V4[unique(TE_nORFoverlap@to)])
TEwithnORFs <- as.character(TEcoord$Geneid[unique(TE_nORFoverlap@from)])
nORFexpr<-nORFexpr0[,colnames(nORFexpr0)%in%nORFswithTEs]
TEexpr<-TEexpr0[,colnames(TEexpr0)%in%TEwithnORFs]
log_allexpr<-log2(cbind(nORFexpr,TEexpr)+1)
normalitytest<-normality(log_allexpr)

# Correlations and identification of strong correlations
spearmanCorr<-cor(log2(TEexpr+1),log2(nORFexpr+1),method = "spearman")
pearsonCorr<-cor(log2(TEexpr+1),log2(nORFexpr+1))
spearmanPvalue<-corPvalueStudent(spearmanCorr,length(sampleNames))
spearmanPvalue<-matrix(p.adjust(as.vector(as.matrix(spearmanPvalue)),method = "fdr"),ncol=ncol(spearmanPvalue))
pearsonPvalue<-corPvalueStudent(pearsonCorr,length(sampleNames))
pearsonPvalue<-matrix(p.adjust(as.vector(as.matrix(pearsonPvalue)),method = "fdr"),ncol=ncol(pearsonPvalue))
# Strong correlations of TE-nORFs and nORF-TEs only
TEtonORF<-data.frame(TE=as.character(TEcoord$Geneid[TE_nORFoverlap@from]),
                     Class=TEcoord$Class[TE_nORFoverlap@from],
                     Family=TEcoord$Family[TE_nORFoverlap@from],
                     nORF=allnORFs$V4[TE_nORFoverlap@to],stringsAsFactors = FALSE)
TEtonORF<-TEtonORF[!duplicated(TEtonORF),]
TEtonORF<-TEtonORF[!is.na(match(TEtonORF$TE,row.names(spearmanCorr))),]
TEtonORF<-TEtonORF[!is.na(match(TEtonORF$nORF,colnames(spearmanCorr))),]
CorrIndex1<-data.frame(row=match(TEtonORF$TE,row.names(spearmanCorr)),col=match(paste(TEtonORF$nORF),paste(colnames(spearmanCorr))))
TEtonORF$r_s<-spearmanCorr[as.matrix(CorrIndex1)]
TEtonORF$p_s<-spearmanPvalue[as.matrix(CorrIndex1)]
TEtonORF$r_p<-pearsonCorr[as.matrix(CorrIndex1)]
TEtonORF$p_p<-pearsonPvalue[as.matrix(CorrIndex1)]
TEtonORF<-subset(TEtonORF,abs(r_s)>0.5 | abs(r_p)>0.5)
# All strong correlations - to avoid crashing R, pick a suitably high number for "chunks"
chunks=1000
quotient<-length(names(TEexpr0))%/%chunks
remainder<-length(names(TEexpr0))%%chunks
allStrong<-data.frame(TE=character(),
                      nORF=character(),
                      Class=character(),
                      Family=character(),
                      r_s=numeric(),
                      p_s=numeric(),
                      r_p=numeric(),
                      p_p=numeric(),
                      stringsAsFactors = FALSE)
for (val in c(1:chunks)) {
  startcol = (val-1)*quotient+1
  if (val < chunks) {
    endcol = val*quotient
  } else {
    endcol = val*quotient+remainder
  }
  colRange = names(TEexpr0)[c(startcol:endcol)]
  TEexprSubset = TEexpr0[,colRange]
  spearmanCorrSub = cor(log2(TEexprSubset+1),log2(nORFexpr0+1),method = "spearman")
  spearmanPvalueSub = corPvalueStudent(spearmanCorrSub,length(sampleNames))
  spearmanPvalueSub = matrix(p.adjust(as.vector(as.matrix(spearmanPvalueSub)),method = "fdr"),ncol=ncol(spearmanPvalueSub))
  pearsonCorrSub = cor(log2(TEexprSubset+1),log2(nORFexpr0+1))
  pearsonPvalueSub = corPvalueStudent(pearsonCorrSub,length(sampleNames))
  pearsonPvalue = matrix(p.adjust(as.vector(as.matrix(pearsonPvalueSub)),method = "fdr"),ncol=ncol(pearsonPvalueSub))
  spearmanIndices = data.frame(which(abs(spearmanCorrSub)>0.8,arr.ind = TRUE))
  spearmanFrame = data.frame(TE=row.names(spearmanCorrSub)[spearmanIndices$row],
                             Class=TEcoord$Class[match(row.names(spearmanCorrSub)[spearmanIndices$row],TEcoord$Geneid)],
                             Family=TEcoord$Family[match(row.names(spearmanCorrSub)[spearmanIndices$row],TEcoord$Geneid)],
                             nORF=colnames(spearmanCorrSub)[spearmanIndices$col],
                             r_s=spearmanCorrSub[as.matrix(spearmanIndices)],
                             p_s=spearmanPvalueSub[as.matrix(spearmanIndices)],
                             r_p=pearsonCorrSub[as.matrix(spearmanIndices)],
                             p_p=pearsonPvalueSub[as.matrix(spearmanIndices)])
  pearsonIndices = data.frame(which(abs(pearsonCorrSub)>0.8,arr.ind = TRUE))
  pearsonFrame = data.frame(TE=row.names(pearsonCorrSub)[pearsonIndices$row],
                            Class=TEcoord$Class[match(row.names(pearsonCorrSub)[pearsonIndices$row],TEcoord$Geneid)],
                            Family=TEcoord$Family[match(row.names(pearsonCorrSub)[pearsonIndices$row],TEcoord$Geneid)],
                             nORF=colnames(pearsonCorrSub)[pearsonIndices$col],
                             r_s=spearmanCorrSub[as.matrix(pearsonIndices)],
                             p_s=spearmanPvalueSub[as.matrix(pearsonIndices)],
                             r_p=pearsonCorrSub[as.matrix(pearsonIndices)],
                             p_p=pearsonPvalueSub[as.matrix(pearsonIndices)])
  combinedFrame = dplyr::union(spearmanFrame,pearsonFrame)
  allStrong = dplyr::union(allStrong,combinedFrame)
}
write.csv(allStrong,"allStrongCorrelations.csv")
