gomat <- read.csv("/media/david/genomicsdrive/NovelPeptideCharacterisation/ps_scan/gopspair.csv",header=T)
gomat <- gomat[apply(gomat,1,sum)!=0,]
gomat <- gomat[,apply(gomat,2,sum)!=0]

data <- read.csv("/media/david/genomicsdrive/NovelPeptideCharacterisation/HS_genbank_GRCh38_v1-prosite.csv",header=T)
data <- data[,c(T,apply(data[,2:ncol(data)],2,sum)!=0)]
data <- cbind(data[,1],data[,!is.na(match(colnames(data),rownames(gomat)))])
data <- data[apply(data[,2:ncol(data)],1,sum)!=0,]

goindices <- match(colnames(data[,2:ncol(data)]),rownames(gomat))
temp <- apply(data[,2:ncol(data)],1,FUN = function(x) goindices[x>0])
goterms <- unique(unlist(apply(gomat[unique(unlist(temp)),],1,function(x) colnames(gomat)[x>0])))