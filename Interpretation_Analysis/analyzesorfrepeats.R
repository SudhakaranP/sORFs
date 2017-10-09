setwd("/media/david/genomicsdrive/NovelPeptideCharacterisation/")
data <- read.csv("data/8175_allinfo.csv",header=T)
malesorfs <- data[(!is.na(data$tpm_abund_bm) | !is.na(data$tpm_abund_tm | !is.na(data$prot_abund_bm) | !is.na(data$prot_abund_tm))) & (is.na(data$tpm_abund_tf) & is.na(data$tpm_abund_bf) & is.na(data$prot_abund_tf) & is.na(data$prot_abund_bf)),c(1,12,14,20,22)]
femalesorfs <- data[(!is.na(data$tpm_abund_bf) | !is.na(data$tpm_abund_tf | !is.na(data$prot_abund_bf) | !is.na(data$prot_abund_tf))) & (is.na(data$tpm_abund_tm) & is.na(data$tpm_abund_bm) & is.na(data$prot_abund_tm) & is.na(data$prot_abund_bm)),c(1,13,15,21,23)]
tcellsorfs <- data[(!is.na(data$tpm_abund_tf) | !is.na(data$tpm_abund_tm | !is.na(data$prot_abund_tf) | !is.na(data$prot_abund_tm)))  & (is.na(data$tpm_abund_bm) & is.na(data$tpm_abund_bf) & is.na(data$prot_abund_bm) & is.na(data$prot_abund_bf)),c(1,14,15,22,23)]
bcellsorfs <- data[(!is.na(data$tpm_abund_bf) | !is.na(data$tpm_abund_bm | !is.na(data$prot_abund_bf) | !is.na(data$prot_abund_bm)))  & (is.na(data$tpm_abund_tm) & is.na(data$tpm_abund_tf) & is.na(data$prot_abund_tm) & is.na(data$prot_abund_tf)),c(1,12,13,20,21)]
transcriptsorfs <- data[!is.na(data$tpm_abund_bf) | !is.na(data$tpm_abund_bm) | !is.na(data$tpm_abund_tf) | !is.na(data$tpm_abund_tm),c(1,12:15)]
translatesorfs <- data[!is.na(data$prot_abund_bf) | !is.na(data$prot_abund_bm) | !is.na(data$prot_abund_tf) | !is.na(data$prot_abund_tm),c(1,20:23)]
data <- read.csv("sORFsTranslationAndTranscription-mm10repeatscoordsandclasses-intersect.csv", header=F)
repeatstartpos <- apply(data,1,FUN=function(x) (as.numeric(x[7])-as.numeric(x[2])) / (as.numeric(x[3])-as.numeric(x[2])))
repeatendpos <- apply(data,1,FUN=function(x) (as.numeric(x[8])-as.numeric(x[3])) / (as.numeric(x[3])-as.numeric(x[2])))
data <- cbind(data,repeatstartpos)
data <- cbind(data,repeatendpos)
plotsorfoverlap <- function(sorfid) {
  temp <- data[data[,5]==sorfid,]
  repeats <- apply(temp,MARGIN=1,FUN=function(x) {c(as.numeric(x[7]),as.numeric(x[8]))})
  repeats <- cbind(repeats,c(as.numeric(temp[1,2]),as.numeric(temp[1,3])))
  colnames(repeats) <- c(as.character(temp[,11]),as.character(sorfid))
  boxplot(repeats, horizontal = T, cex.names=0, main=sorfid,xlab=paste(temp[1,1],"coords"), las=1, cex.axis=0.5)
}

isunique <- function(repeattype) {
   <- 
  
}

categories <- c("DNA","LINE","Low_complexity","LTR","RC","RNA","SINE","Simple_repeat","Unknown","Other")
uniquesorfs <- unique(data[,5])
uniqueltrsorfs <- unique(data[grep("LTR",data[,11]),5])

print(paste("Number of repeat classes:",length(summary(data[,10]))))
print(paste("Number of repeat classes with 1 instance: ",sum(summary(data[,10])==1)))
print(paste("Number of sorfs with an LTR:",length(uniqueltrsorfs)))
print(paste("Number of sorfs with an LTR from Bcell", sum(!is.na(match(bcellsorfs[,1],uniqueltrsorfs)))))
print(paste("Number of sorfs with an LTR from Tcell", sum(!is.na(match(tcellsorfs[,1],uniqueltrsorfs)))))
print(paste("Number of sorfs with an LTR from male cells", sum(!is.na(match(malesorfs[,1],uniqueltrsorfs)))))
print(paste("Number of sorfs with an LTR from female cells", sum(!is.na(match(femalesorfs[,1],uniqueltrsorfs)))))

sorfrepeatprofile <- sapply(uniqueltrsorfs,FUN=function(x) as.character(data[data[,5]==x,10]))
uniquerepeatsinsorfs <- sapply(levels(data[,10]),FUN=function(y) sapply(sorfrepeatprofile,FUN=function(x) if(y %in% x){return(1);} else{return(0);}))
summary(apply(uniquerepeatsinsorfs,1,sum) == 1)
uniqueltrsorfs[apply(uniquerepeatsinsorfs,1,sum) == 1]
specialsorfs <- transcriptsorfs[!is.na(match(transcriptsorfs[,1],translatesorfs[,1])),]

#pdf("/media/david/genomicsdrive/NovelPeptideCharacterisation/analysis/sORFrepeatsgraph.pdf")
#barplot(sapply(categories,FUN=function(x){length(grep(x,data[),11]))}), las=3, cex.names=0.7, ylab="Freq",main="Count of repeat types")
#barplot(summary(data[,12]), cex.names=0.7, main="Type of repeat overlap with reference to sORF", ylab="Freq")
#sapply(uniquesorfs,FUN=plotsorfoverlap)
#dev.off()