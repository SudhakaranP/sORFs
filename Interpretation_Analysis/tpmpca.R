sorfs <- read.table("/media/david/genomicsdrive/NovelPeptideCharacterisation/data/4667_sorf_samplewise_tpmvalues.txt",header=T)
known <- read.table("/media/david/genomicsdrive/NovelPeptideCharacterisation/data/knowntranscripts_tpmvalues.txt", header=T)
bcelltranscripts <- read.csv("/media/david/genomicsdrive/NovelPeptideCharacterisation/Bcells_transcripts.csv",header=F)
tcelltranscripts <- read.csv("/media/david/genomicsdrive/NovelPeptideCharacterisation/Tcells_transcripts.csv",header=F)
known <- known[is.na(match(known[,1],bcelltranscripts[,1])&match(known[,2],bcelltranscripts[,2])&match(known[,3],bcelltranscripts[,3])),]
known <- known[is.na(match(known[,1],tcelltranscripts[,1])&match(known[,2],tcelltranscripts[,2])&match(known[,3],tcelltranscripts[,3])),]

sorfs <- t(sorfs)
sorfs <- as.numeric(sorfs[2:13,])
sorfs <- matrix(sorfs, nrow=12)
sorfs <- sorfs + 1
sorfs <- log2(sorfs)

known <- t(known)
known <- as.numeric(known[5:16,])
known <- matrix(known, nrow=12)
known <- apply(known,2,function(x) sapply(x,function(y) if(is.na(y)){return(0);}else{return(y)}))
known <- known + 1
known <- log2(known)

sorfmedians <- apply(sorfs,2,median)
knownmedians <- apply(known,2,median)

sorfs <- apply(sorfs,1,FUN=function(x) x-sorfmedians)
known <- apply(known,1,FUN=function(x) x-knownmedians)

sorfs <- t(sorfs)
known <- t(known)

knownpca <- prcomp(known, scale=F, center=F)
sorfspca <- prcomp(sorfs, scale=F, center=F)

BT <- c(rep("red",6),rep("black",6))
mf <- c(rep("red",3),rep("black",3),rep("red",3),rep("black",3))

png("/media/david/genomicsdrive/NovelPeptideCharacterisation/analysis/BTknowntranscriptpca.png")
plot(knownpca$x[,1],knownpca$x[,2],col=BT, xlab="PC1", ylab="PC2",main="PCA of known transcript expression levels")
legend(-125,75,c("B cell","T cell"),col=c("red","black"),pch=1)
dev.off()

png("/media/david/genomicsdrive/NovelPeptideCharacterisation/analysis/MFknowntranscriptpca.png")
plot(knownpca$x[,1],knownpca$x[,2],col=mf, xlab="PC1", ylab="PC2",main="PCA of known transcript expression levels")
legend(-125,75,c("male","female"),col=c("red","black"),pch=1)
dev.off()

png("/media/david/genomicsdrive/NovelPeptideCharacterisation/analysis/BTsorftranscriptpca.png")
plot(sorfspca$x[,1],sorfspca$x[,2],col=BT, xlab="PC1", ylab="PC2",main="PCA of sORF transcript expression levels")
legend(-15,-15,c("B cell","T cell"),col=c("red","black"),pch=1)
dev.off()

png("/media/david/genomicsdrive/NovelPeptideCharacterisation/analysis/MFsorftranscriptpca.png")
plot(sorfspca$x[,1],sorfspca$x[,2],col=mf, xlab="PC1", ylab="PC2",main="PCA of sORF transcript expression levels")
legend(-15,-15,c("male","female"),col=c("red","black"),pch=1)
dev.off()

sorfTstats <- apply(sorfspca$x,1,function(x) sum((x/sorfspca$sdev)^2))