library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)
library(plyr)
getloc <- function(x) {return(c(nrow(x[x$Loc=="_",]),nrow(x[x$Loc=="M",]),nrow(x[x$Loc=="S",])))}
getmostprobable <- function(x) {match(max(x[seq(2,10)]),x)} 
eucdist <- function(x,y) {return(sqrt(sum((x-y)^2)))}

setwd("/media/david/genomicsdrive/NovelPeptideCharacterisation/analysis")
#data <- read.csv(header=T,"../HS_genbank_GRCh38_v1-loc.csv")
#multiloc <- read.csv(header=T,"../HS_genbank_GRCh38_v1-multiscore.csv")
#altorfdist <- count(as.factor(read.csv(header=F,"../data/HS_genbank_GRCh38_v1.2_altorf_gene_dist.csv")$V1))

data <- read.csv(header=T,"../PrabakaransORFs8175-loc.csv")
multiloc <- read.csv(header=T,"../PrabakaransORFs8175-multiscore.csv")

#pdf("altorf_subcellular_location_analysis.pdf")
pdf("PrabakaransORFs8175_subcellular_location_analysis.pdf")
RC1 <- data[data$RC=="1",]
RC2 <- data[data$RC=="2",]
RC3 <- data[data$RC=="3",]
RC4 <- data[data$RC=="4",]
RC5 <- data[data$RC=="5",]
RCheaders <- c("RC1","RC2","RC3","RC4","RC5")
#Summary of RC dist
RCdist <- data.frame(cbind(matrix(RCheaders,ncol=1),rbind(nrow(RC1),nrow(RC2),nrow(RC3),nrow(RC4),nrow(RC5))))
colnames(RCdist) <- c("RC","Freq")
RCdistbp <- ggplot(RCdist, aes(x=RC,y=Freq,label=Freq)) + geom_col() + geom_text(nudge_y=0.5) + labs(title="altORF RC Score Distribution (TargetP)",y="Frequency",x="RC Score")
#Summary of RC vs altORF length
RClens <- data.frame(quantile(RC1$Len),quantile(RC2$Len),quantile(RC3$Len),quantile(RC4$Len),quantile(RC5$Len))
colnames(RClens) <- c("RC1 Seq Len","RC2 Seq Len","RC3 Seq Len","RC4 Seq Len","RC5 Seq Len")
rownames(RClens) <- c("Min","Lower","Median","Upper","Max")
RClenstbl <- tableGrob(RClens)
grid.arrange(RCdistbp, RClenstbl, ncol=1)
#Subcellular locations
altORFlocs <- data.frame(matrix(c(c("Other","Mito","Secrete"),getloc(RC1),getloc(RC2),getloc(RC3),getloc(RC4),getloc(RC5)),nrow=3,ncol=6))
colnames(altORFlocs) <- c("Loc", RCheaders)
altORFlocs <- melt(altORFlocs,id.vars=1)
altORFlocs$value <- as.numeric(altORFlocs$value)
altORFlocsbp <- ggplot(altORFlocs, aes(x=Loc, y=value)) + geom_col(aes(fill=variable),position="dodge") + labs(title="altORF Subcellular Locations by RC score (TargetP)", x="Subcellular Location", y="Freq")
plot(altORFlocsbp)
#Multiloc and distribution of altorfs on genes
#altorfdistdf <- data.frame(quantile(altorfdist$freq))
#colnames(altorfdistdf) <- c("altORF count per Gene")
#rownames(altorfdistdf) <- c("Min","Lower","Median","Upper","Max")
#altorfdisttbl <- tableGrob(altorfdistdf)
#altorfdisthist <- ggplot(altorfdist, aes(x=freq)) + geom_histogram(breaks=seq(0,50)) + labs(title="Histogram of altORF per Gene (up to 50 altORFs)",x="altORF count per Gene",y="Frequency")
#grid.arrange(altorfdisthist,altorfdisttbl,ncol=2)

colnames(multiloc) <- c("ID","Nuc","Cyt","Per","Mit","Gol","Pl.m","ER","E.cel","Lys")
mostprobable <- colnames(multiloc)[apply(multiloc,1,getmostprobable)]
mostprobabledist <- count(mostprobable)
mostprobablebp <- ggplot(mostprobabledist, aes(x=x,y=freq,label=freq)) + geom_col() + geom_text(nudge_y=0.5) + labs(title="Most probable location of sORF peptides(Multiloc2)",y="Frequency",x="Subcellular Location")
grid.arrange(mostprobablebp,ncol=1)
dev.off()

