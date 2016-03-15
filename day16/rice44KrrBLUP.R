# rice data: http://ricediversity.org/data/sets/44kgwas/

##################################
# data cleaning 
##################################
library(BGLR)
out<- read_ped("RiceDiversity_44K_Genotypes_PLINK/sativas413.ped")
p=out$p
n=out$n
out=out$x
#Recode snp to 0,1,2 format using allele 1
# 0 --> 0
# 1 --> 1
# 2 --> NA
# 3 --> 2
out[out==2]=NA
out[out==3]=2
X=matrix(out,nrow=p,ncol=n,byrow=TRUE)
X=t(X)

# genotype imputation
for (j in 1:ncol(X)){
  X[,j] <- ifelse(is.na(X[,j]), mean(X[,j], na.rm=TRUE), X[,j])
}

# accession ID
fam <-read.table("RiceDiversity_44K_Genotypes_PLINK/sativas413.fam", header = FALSE, stringsAsFactors = FALSE)  
head(fam)
rownames(X) <- fam$V2 # 413 x 36901

# phenotypes
rice.pheno <- read.table("http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt", header=TRUE, stringsAsFactors = FALSE, sep = "\t")
table(rownames(X) == rice.pheno$NSFTVID)
y <- rice.pheno$Flowering.time.at.Arkansas # # use the first trait 
index <- !is.na(y)
y <- y[index] # 374
X <- X[index,] # 374x36901


##################################
# mixed.solve()
# ridge regression BLUP 
##################################
library(rrBLUP)
fit <- mixed.solve(y = y, Z=X)
# marker additive genetic variance
fit$Vu
# residual variance
fit$Ve
# intercept 
fit$beta
# marker effects 
fit$u


#########################################################
# GWAS()
# single marker regression + polygenic effects (G matrix) 
##########################################################
# create data frame pheno 
my.pheno <- data.frame(NSFTV_ID=rice.pheno[,2], y=rice.pheno[,3]) # use the first trait 
my.pheno <- my.pheno[index,]
table(rownames(X) == my.pheno[,1])
# read map file
map <-read.table("RiceDiversity_44K_Genotypes_PLINK/sativas413.map", header = FALSE, stringsAsFactors = FALSE)  
my.geno <- data.frame(marker=map[,2], chrom=map[,1], pos=map[,4], t(X-1), check.names=FALSE) # X = \in{-1, 0, 1}
scores <- GWAS(my.pheno, my.geno, fixed=min.MAF=0.05, P3D=TRUE, plot=FALSE)
scores[,4] # returns -log10p




