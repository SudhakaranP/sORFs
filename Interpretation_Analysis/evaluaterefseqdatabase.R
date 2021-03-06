data <- read.csv("/media/david/genomicsdrive/NovelPeptideCharacterisation/PrabakaransORFs2-PrabakaransORFs2refseq-compare.csv",header=T)

identicalcouplings <- data[data$skippedcouplings == 0,]
identicalcouplingsandstructpredicted <- data[data$skippedcouplings == 0 & !is.na(data$rmsd),]
bothpredicted <- data[!is.na(data$struct1) & !is.na(data$struct2),]

print(paste("Number of structures predicted using uniprot100:", summary(!is.na(data$struct1))["TRUE"]))
print(paste("Number of structures predicted using refseq transcript:",summary(!is.na(data$struct2))["TRUE"]))
print(paste("Number of structures predicted by both:",summary(!is.na(data$struct1) & !is.na(data$struct2))["TRUE"]))
print(paste("Number of structures with no prediction:",summary(is.na(data$struct1) & is.na(data$struct2))["TRUE"]))
print(paste("Average alphabetascore for uniref100 structures:",mean(data$alphabetascore1,na.rm=T)))
print(paste("Average alphabetascore for refseq structures:",mean(data$alphabetascore2,na.rm=T)))
print(paste("Average RMSD:",mean(data$rmsd,na.rm=T)))
print(paste("Average Normalised Coupling Rank Change:",mean(identicalcouplings$normalisedcouplingrank,na.rm=T)))
print(paste("Average Coupling Rank Difference:",mean(identicalcouplings$couplingsrank,na.rm=T)))
print(paste("Average Couplings Probability Difference:",mean(identicalcouplings$couplingsprob,na.rm=T)))

png("/media/david/genomicsdrive/NovelPeptideCharacterisation/documents/rmsdalphabeta.png")
plot(abs(bothpredicted$alphabetascore1-bothpredicted$alphabetascore2),bothpredicted$rmsd, xlab="AlphaBeta Score Difference", ylab="RMSD", main="Relationship between RMSD and AlphaBeta Score\n Difference of Structures Predicted with \nUniref100 and Refseq Transcript")
fit <- lm(bothpredicted$rmsd~abs(bothpredicted$alphabetascore1-bothpredicted$alphabetascore2))
fitcor <- cor(abs(bothpredicted$alphabetascore1-bothpredicted$alphabetascore2),bothpredicted$rmsd)
abline(fit)
text(0.3,27.5,paste("Cor",fitcor),adj=c(0,0.5))
text(0.21,25,paste0("y=",fit$coefficients[2],"x+",fit$coefficients[1]),adj=c(0,0.5))
dev.off()

png("/media/david/genomicsdrive/NovelPeptideCharacterisation/documents/rmsdcouplings.png")
plot(identicalcouplings$couplingsprob,identicalcouplings$rmsd, xlab="Average Difference in Couplings Probability", ylab="RMSD", main="Relationship between RMSD and Average Difference\n in Coupling Probabilities of Structures Predicted with \nUniref100 and Refseq Transcript and identical couplings")
fit <- lm(identicalcouplings$rmsd~identicalcouplings$couplingsprob)
fitcor <- cor(identicalcouplings$rmsd,identicalcouplings$couplingsprob)
abline(fit)
text(0.075,20,paste("Cor",fitcor),adj=c(0,0.5))
text(0.06,17.5,paste0("y=",fit$coefficients[2],"x+",fit$coefficients[1]),adj=c(0,0.5))
dev.off()

png("/media/david/genomicsdrive/NovelPeptideCharacterisation/documents/rmsdcouplingrank.png")
plot(identicalcouplingsandstructpredicted$couplingsrank,identicalcouplingsandstructpredicted$rmsd, xlab="Average Change in Couplings Rank ", ylab="RMSD", main="Relationship between RMSD and Average Change\nin Coupling Rank of Structures Predicted with \nUniref100 and Refseq Transcript and identical couplings")
fit <- lm(identicalcouplingsandstructpredicted$rmsd~identicalcouplingsandstructpredicted$couplingsrank)
fitcor <- cor(identicalcouplingsandstructpredicted$rmsd,identicalcouplingsandstructpredicted$couplingsrank)
abline(fit)
text(100,25,paste("Cor",fitcor),adj=c(0,0.5))
text(100,22.5,paste0("y=",fit$coefficients[2],"x+",fit$coefficients[1]),adj=c(0,0.5))
dev.off()

png("/media/david/genomicsdrive/NovelPeptideCharacterisation/documents/rmsdcouplingranknormalised.png")
plot(identicalcouplingsandstructpredicted$normalisedcouplingrank,identicalcouplingsandstructpredicted$rmsd, xlab="Average Normalised Change in Couplings Rank ", ylab="RMSD", main="Relationship between RMSD and Average Normalised\nChangein Coupling Rank of Structures Predicted with \nUniref100 and Refseq Transcript and identical couplings")
fit <- lm(identicalcouplingsandstructpredicted$rmsd~identicalcouplingsandstructpredicted$normalisedcouplingrank)
fitcor <- cor(identicalcouplingsandstructpredicted$rmsd,identicalcouplingsandstructpredicted$normalisedcouplingrank)
abline(fit)
text(0.125,25,paste("Cor",fitcor),adj=c(0,0.5))
text(0.125,22.5,paste0("y=",fit$coefficients[2],"x+",fit$coefficients[1]),adj=c(0,0.5))
dev.off()

png("/media/david/genomicsdrive/NovelPeptideCharacterisation/documents/normalisedcouplingrankhist.png")
hist(data$normalisedcouplingrank, xlab="Average Change in Normalised Coupling Rank",main="Distribution of average change in normalised coupling\nrank for couplings calculated using Uniref100 and refseq")
dev.off()

png("/media/david/genomicsdrive/NovelPeptideCharacterisation/documents/couplingrankhist.png")
hist(data$couplingsrank, xlab="Average Change in Coupling Rank",main="Distribution of average change in coupling rank\nfor couplings calculated using Uniref100 and refseq")
dev.off()