library(GenomicRanges)
library(tidyr)
library(dplyr)
# Prep of files for INRICH analysis
# SNP files were renamed for ease of use and identification
bip18snps<-read.table('bip18snps',header = TRUE)
bip18snpsI<-bip18snps[,c(1,3)]
write.table(bip18snpsI,"bip18snps.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
scz18snps<-read.table('scz18snps.txt',header = TRUE)
scz18snpsI<-scz18snps[,c(3,4)]
write.table(scz18snpsI,"scz_rsf.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
refnORFs <- read.table("norf_coordinates.bed")
refnORFs <- refnORFs %>% distinct(V1,V2,V3, .keep_all=TRUE)
refnORFs <- separate(refnORFs,V1,into = c("V1","chr"),sep="chr")
refnORFs$chr <- sub("X","23",refnORFs$chr)
refnORFsa <- refnORFs[,c(2,3,4,5,5)]
write.table(refnORFsa,"INRICH/norfs_rgf.gene.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
HARlist <- read.table("mergedHARs.txt",header = TRUE)
HARlist <- HARlist %>% separate(Chr,into = c("N","chr"),sep="chr")
HARlist <- HARlist[,-1]
HARlist$chr <- sub("X","23",HARlist$chr)
HARR <- GRanges(seqnames = Rle(HARlist$chr,1),ranges = IRanges(HARlist$Start-100000,HARlist$End+100000,names = HARlist$ID))
nORFR <- GRanges(seqnames = Rle(refnORFsa$chr,1),ranges = IRanges(refnORFsa$V2,refnORFsa$V3,names = refnORFsa$V4))
HAN <- findOverlaps(HARR,nORFR)
HaAsNo<-data.frame(nORF=refnORFsa$V4[HAN@to],HARs=HARlist$Overlaps[HAN@from])
HaAsNo<- HaAsNo %>% group_by(nORF) %>% summarise(HARs=paste(HARs,collapse = ","))
HaAsNo$vHARs<-grepl("Pol|PRA|ANC",HaAsNo$HARs)
HaAsNo$mHARs<-grepl("BUSH|2xH",HaAsNo$HARs)
HaAsNo$pHARs<-grepl("2xP|haD",HaAsNo$HARs)
NvHARs<-data.frame(gene_id=HaAsNo$nORF[HaAsNo$vHARs],gene_set_id="vHARs",gene_set_description="nORFs near vHARs")
NmHARs<-data.frame(gene_id=HaAsNo$nORF[HaAsNo$mHARs],gene_set_id="mHARs",gene_set_description="nORFs near mHARs")
NpHARs<-data.frame(gene_id=HaAsNo$nORF[HaAsNo$pHARs],gene_set_id="pHARs",gene_set_description="nORFs near pHARs")
NHARs<-data.frame(gene_id=HaAsNo$nORF,gene_set_id="HARs",gene_set_description="nORF near HARs")
write.table(NvHARs,"INRICH/NvHARs.set",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(NmHARs,"INRICH/NmHARs.set",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(NpHARs,"INRICH/NpHARs.set",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(NHARs,"INRICH/NHARs.set",quote = FALSE,row.names = FALSE,col.names = FALSE)

# Prep of files for PLINK
james<- bip18snps[-(bip18snps$INFO<0.9),]
jamesc<-james[,c("SNP","P")]
write.table(jamesc,"PLINK/bip18pvalclean",quote = FALSE,row.names = FALSE)
iago<-separate(scz18snps,SNP,c("SNP",2,3,4),sep="\\:")
iago<-iago[,-c(2:4)]
iagoc<-iago[,c("SNP","P")]
write.table(iagoc,"PLINK/scz18pvalclean",quote = FALSE,row.names = FALSE)
# Clump in PLINK then return
bipclump<-read.table("PLINK/bip_myclump1.clumped",header = TRUE)
sczclump0<-read.table("PLINK/scz_myclump1.clumped",header = TRUE)
# Remove all MHC regions snps but 1
MHCscz<-subset(sczclump0,CHR==6 & BP>28477797 & BP<33448354)
MHCscz<-MHCscz[-which.min(MHCscz$P),]
removeSNPs<-as.numeric(row.names(MHCscz))
sczclump<-sczclump0[-c(removeSNPs),]
# SNPlists clumped
# 1e-2
write.table(bipclump[bipclump$P<1e-2,c("SNP")],"PLINK/bipsnplist2",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-2,c("SNP")],"PLINK/sczsnplist2",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(bipclump[bipclump$P<1e-2,c("CHR","BP")],"PLINK/bip2.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-2,c("CHR","BP")],"PLINK/scz2.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
# 1e-3
write.table(bipclump[bipclump$P<1e-3,c("SNP")],"PLINK/bipsnplist3",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-3,c("SNP")],"PLINK/sczsnplist3",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(bipclump[bipclump$P<1e-3,c("CHR","BP")],"PLINK/bip3.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-3,c("CHR","BP")],"PLINK/scz3.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
# 1e-4
write.table(bipclump[bipclump$P<1e-4,c("SNP")],"PLINK/bipsnplist4",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-4,c("SNP")],"PLINK/sczsnplist4",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(bipclump[bipclump$P<1e-4,c("CHR","BP")],"PLINK/bip4.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-4,c("CHR","BP")],"PLINK/scz4.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
# 1e-5
write.table(bipclump[bipclump$P<1e-5,c("SNP")],"PLINK/bipsnplist5",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-5,c("SNP")],"PLINK/sczsnplist5",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(bipclump[bipclump$P<1e-5,c("CHR","BP")],"PLINK/bip5.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-5,c("CHR","BP")],"PLINK/scz5.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
# 1e-6
write.table(bipclump[bipclump$P<1e-6,c("SNP")],"PLINK/bipsnplist6",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-6,c("SNP")],"PLINK/sczsnplist6",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(bipclump[bipclump$P<1e-6,c("CHR","BP")],"PLINK/bip6.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-6,c("CHR","BP")],"PLINK/scz6.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
# 1e-7
write.table(bipclump[bipclump$P<1e-7,c("SNP")],"PLINK/bipsnplist7",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-7,c("SNP")],"PLINK/sczsnplist7",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(bipclump[bipclump$P<1e-7,c("CHR","BP")],"PLINK/bip7.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(sczclump[sczclump$P<1e-7,c("CHR","BP")],"PLINK/scz7.snp.map",quote = FALSE,row.names = FALSE,col.names = FALSE)