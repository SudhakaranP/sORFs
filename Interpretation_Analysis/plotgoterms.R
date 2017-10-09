term <- "exonicsorfsGOTermsList"

all <- read.csv(paste0(term,".csv"))
cc <- all[all[,1]=="GOTERM_CC_DIRECT",]
mf <- all[all[,1]=="GOTERM_MF_DIRECT",]
bp <- all[all[,1]=="GOTERM_BP_DIRECT",]

cc <- cc[order(cc[,5]),]
mf <- mf[order(mf[,5]),]
bp <- bp[order(bp[,5]),]

topcc <- cc[1:10,5]
topmf <- mf[1:10,5]
topbp <- bp[1:10,5]

names(topcc) <- cc[1:10,2]
names(topmf) <- mf[1:10,2]
names(topbp) <- bp[1:10,2]

png(paste0(term,"cc.png"))
par(mar=c(5.1,10,4.1,2.1))
barplot(-log(topcc), ylab="-log(benjamini)", main="Top 10 Cellular Components of genes\nwith exonic sorfs expressed overall", las=2, cex.names=0.6, horiz = T)
dev.off()
png(paste0(term,"mf.png"))
par(mar=c(5.1,10,4.1,2.1))
barplot(-log(topmf), ylab="-log(benjamini)", main="Top 10 Molecular Functions of genes\nwith exonic sorfs expressed overall", las=2, cex.names=0.6, horiz = T)
dev.off()
png(paste0(term,"bp.png"))
par(mar=c(5.1,10,4.1,2.1))
barplot(-log(topbp), ylab="-log(benjamini)", main="Top 10 Biological Profile of genes\nwith exonic sorfs expressed overall", las=2, cex.names=0.6, horiz = T)
dev.off()