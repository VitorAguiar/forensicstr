# This script calculates confidence intervals for the inbreeding coefficient

if(!"compiler"%in%rownames(installed.packages())) {install.packages("compiler", repos="http://cran.rstudio.com/")}
library(compiler)
enableJIT(3)

# ObsMat and ObsGen are functions to create a list of matrices of genotype frequencies at each locus 
ObsMat <- function(grp, dat) {
	cols <- grep(grp, names(dat))
	Names <- levels(factor(c(dat[ ,cols[1]], dat[ ,cols[2]])))
	table(factor(dat[ ,cols[1]], levels = Names[Names > 0]), factor(dat[ ,cols[2]], levels = Names[Names > 0]))
}

ObsGen <- function(dat, Names) {
	Obsdf <- lapply(Names, ObsMat, data.frame(dat))
	ObsdfD <- lapply(Obsdf, function(tab) array(as.double(tab), dim = dim(tab), dimnames = dimnames(tab)))
	names(ObsdfD) <- Names
	return(ObsdfD)
}

# Function to apply a threshold of minimum allele frequency
# Alleles with frequency below count of 5 or a certain percentage are turned into "-5"
## which is a code for missing data that is turned into NA afterwards
popMAF <- function(dat) {
	datx <- as.matrix(dat)
	
	for (i in seq(1, ncol(datx), 2)) {
		allcount <- table(datx[ ,i:(i + 1)])
		ac <- allcount[as.numeric(names(allcount)) > 0]
		bMAF <- as.numeric(names(ac[ac <= 5]))
		datx[datx[ ,(i)] %in% bMAF, i] <- -5
		datx[datx[ ,(i + 1)] %in% bMAF, i + 1] <- -5
	}
	
	return(datx)
}

# function calculates the moment estimator of the single inbreeding coefficient for a locus
finbred <- function(matobs) {
	n <- sum(matobs)
	numer1 <- numer2 <- denom1 <- numeric()
	for(j in rownames(matobs)) {
		
		Puu <- matobs[j,j]/n
		pu <- sum(sum(matobs[j,]), sum(matobs[,j]))/(2*n)
		n1 <- Puu - pu^2
		numer1 <- c(numer1, n1)
		n2 <- pu - Puu
		numer2 <- c(numer2, n2)
		d1 <- pu*(1-pu)
		denom1 <- c(denom1, d1)
	}
	denom2 <- numer1
	
	f <- ( sum(numer1) + sum(numer2)/(2*n) )/( sum(denom1) - sum(denom2)/(2*n) )
	return(round(f,5))
}

# bootstrap to get the confidence intervals for the f coefficient
boot <- function(dat, nreps=10000) {
	
	N <- dim(dat)[1]
	vecFboot <- numeric(nreps)
	
	for (i in 1:nreps) {
		vec <- 1:N
		use <- sample(vec, replace=TRUE)
		
		bootsample <- dat[use, ]
		row.names(bootsample) <- NULL 
		
		bootObsMat <- ObsGen(bootsample, colnames(dat)[1])[[1]]
		
		vecFboot[i] <- finbred(bootObsMat)
	}
	
	# ci
	boot.ci <- round(quantile(vecFboot, c(0.025,0.975)), digits=5)
	return(boot.ci)
}


# Read the data
pop <- read.csv("path_to_input_data.csv", stringsAsFactors=FALSE, row.names=1)

# Take the loci names
lnames <- unique(gsub("(.*)\\.\\d{1}$", "\\1", names(pop)))

nloci <- length(lnames)

# Apply popMAF function to the data
pop <- popMAF(pop)
pop[pop<0] <- NA

# Make a list where each entry is a matrix of observed genotypes for a locus
matobs <- ObsGen(pop, lnames)

# allocate memmory for vector filled in the loop
f <- ci.low <- ci.upp <- numeric(nloci)

L <- 1
for(l in seq(1,ncol(pop)-1,2)) {
	
	dat <- pop[ ,c(l,l+1)]
	dat <- dat[complete.cases(dat),]
	
	f[L] <- finbred(matobs[[L]])
	
	boot.ci <- boot(dat)
	ci.low[L] <- boot.ci[1]
	ci.upp[L] <- boot.ci[2]
	L <- L+1
}

#save the results
cidf <- data.frame(lnames, f, ci.low, ci.upp)
write.table(cidf, "cidf.csv", sep=",", row.names=FALSE)
