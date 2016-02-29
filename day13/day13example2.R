# Pew et al. (2015) related: an R package for analysing pairwise relatedness from codominant molecular markers. Molecular Ecology Resources.
# https://frasierlab.wordpress.com/software/

library(related)
## single locus -  Ritland (1996) 
mydat <- data.frame(ID=as.character(c("1", "2", "3")), SNP1a= c(1, 2, 1), SNP1b=c(1, 1, 1))
mydat[,1] <- as.character(mydat[,1])
mydat
output <- coancestry(mydat , ritland=1)
output$relatedness

## two loci - Ritland (1996) 
mydat <- data.frame(ID=as.character(c("1", "2", "3")), SNP1a= c(1, 2, 1), SNP1b=c(1, 1, 1), SNP2a=c(1,1,1), SNP2b=c(1,2,2))
mydat[,1] <- as.character(mydat[,1])
mydat
output <- coancestry(mydat , ritland=1)
output$relatedness


