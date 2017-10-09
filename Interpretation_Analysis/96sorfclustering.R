refseq <- read.csv("/media/david/genomicsdrive/NovelPeptideCharacterisation/PrabakaransORFs2refseqTopPDBEVfoldrmsdmat.csv",header=T)
uniref <- read.csv("/media/david/genomicsdrive/NovelPeptideCharacterisation/PrabakaransORFs2TopPDBEVfoldrmsdmat.csv",header=T)

refseqcor <- 1 - cor(refseq)
unirefcor <- 1 - cor(uniref)

refseqcorclust <- hclust(as.dist(refseqcor))
unirefcorclust <- hclust(as.dist(unirefcor))