library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)

motifcount <- function(x) {
  motifcountvec <- sapply(as.vector(x),function(y){if(y>0) return(1) else return(0)})
  return(sum(motifcountvec))
}
setwd("/media/david/genomicsdrive/NovelPeptideCharacterisation/analysis")
gopspair <- read.csv(header=T,"../ps_scan/gopspair.csv")

#prosite <- read.csv(header=T,"../HS_genbank_GRCh38_v1-prosite.csv")
#seqlens <- read.csv(header=T,"../HS_genbank_GRCh38_v1-seqlen.csv")

prosite <- read.csv(header=T,"../PrabakaransORFs2-prosite.csv")
seqlens <- read.csv(header=T,"../PrabakaransORFs2-seqlen.csv")

prosite <- na.omit(prosite)

motifabundance <- apply(prosite[,seq(2,2500)],2,sum)
motifabundance <- motifabundance[motifabundance > 0]
motifabundance <- sort(motifabundance,decreasing=T)

#Filter for sequences with likely domains
prosite[,names(motifabundance)[motifabundance > 1000]] <- 0
prosite <- prosite[order(prosite[,1]),]
seqlens <- seqlens[order(seqlens$ID),]
seqlenfield <- seqlens[seqlens$ID %in% prosite[,1],2]
seqidfield <- seqlens[seqlens$ID %in% prosite[,1],1]
motifcountperseq <- apply(prosite[,seq(2,2500)],1,motifcount)
motiflenratio <- motifcountperseq/seqlenfield

summary(motiflenratio>0.05)
summary(motifcountperseq > 2)

topmlr <- cbind(seqidfield[motiflenratio>0.05],seqlenfield[motiflenratio>0.05],motiflenratio[motiflenratio>0.05])
topmlr <- topmlr[order(topmlr[,3],decreasing=T),]

topmcps <- cbind(cbind(seqidfield[motifcountperseq>2],seqlenfield[motifcountperseq>2],motifcountperseq[motifcountperseq>2]))
topmcps <- topmcps[order(topmcps[,3],decreasing=T),]

#write.csv(topmlr,"altorfmotiflenratio.csv")
#write.csv(topmcps,"altorfmotifcount.csv")

write.csv(topmlr,"PrabakaransORFsmotiflenratio.csv")
write.csv(topmcps,"PrabakaransORFsmotifcount.csv")
