if(!"DNAtools"%in%rownames(installed.packages())) {install.packages("DNAtools", repos="http://cran.rstudio.com/")}

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
		# Threshold = 0.01%:
		bMAF <- as.numeric(names(ac[ac < sum(ac) * (0.01/100)]))
		# Threshold = count of 5:
		#bMAF <- as.numeric(names(ac[ac <= 5]))
		datx[datx[ ,(i)] %in% bMAF, i] <- -5
		datx[datx[ ,(i + 1)] %in% bMAF, i + 1] <- -5
	}
	
	return(datx)
}

# Read the data
pop <- read.csv("sample_input.csv", stringsAsFactors=FALSE, row.names=1)

# Turn negative values into missing data
pop[pop<0] <- NA

# Take the loci names
lnames <- unique(gsub("(.*)\\.\\d{1}$", "\\1", names(pop)))

# Number of loci
nloci <- length(lnames)

# Apply popMAF funcion to the data
popMAFdf <- popMAF(pop)

# Create the matrices of genotypes
Obs.list <- ObsGen(popMAFdf, lnames)

# Allocate memory for the vectors and lists filled in the loop below
# nind: N
# vech: observed heterozygosity
# vech.exp: expected heterozygosity
# vecpM: Probability of match
# vecpD: Power of Discrimination
# vecpE: Power of Exclusion
# vecTPI: Typical Paternity Index
# vecPIC: Polymorphic Information Content
# vecf: The Moment estimate of Single Inbreeding coefficient for a locus with more than 2 alleles (Weir's book Data Analysis II page 79)
nind <- vech <- vech.exp <- vecpM <- vecpD <- vecpE <- vecTPI <- vecPIC <- vecf <- numeric(nloci)
#dfall: list of data frames with alleles frequencies at each locus
dfall <- vector("list", nloci)

#initialize counter
L <- 1
# Loop to perform all analyses
while(L <= nloci) {
	
	# Matrix of Observed Genotypes at the current locus
	obs <- as.matrix(Obs.list[[L]])
	
	# Genotypes at the current locus
	currpop <- popMAFdf[ ,(L*2 - 1):(L*2)]
	
	# data for the current locus without missing data
	complete <- currpop
	complete <- complete[complete.cases(complete), ]
	
	# Number of individuals analyzed in the current locus 
	n <- dim(complete)[1]
	nind[L] <- n
	# Allele Frequencies
	alleles <- c(complete[,1], complete[,2])
	alFreq <- table(alleles)/length(alleles)
	
	aldf <- data.frame(as.numeric(names(alFreq)), round(array(alFreq), digits = 6))
	dfall[[L]] <- aldf
	
	# Observed Homozygosity
	H <- sum(diag(obs))/sum(obs)
	
	# Observed and Expected Heterozygosity
	h <- 1 - H
	vech[L] <- h
	
	h.exp <- 1 - sum(alFreq^2)
	vech.exp[L] <- h.exp
	
	# Matching Probability
	pM <- sum((obs/sum(obs))^2)
	vecpM[L] <- pM
	
	# Power of Discrimination
	pD <- 1 - pM
	vecpD[L] <- pD
	
	# Power of Exclusion
	pE <-  (h^2) * (1 - (2*h*H^2))
	vecpE[L] <- pE
	
	# Typical Paternity Index
	TPI <- (H+h)/(2*H)
	vecTPI[L] <- TPI
	
	# Polymorphic information content
	PIC <- 1 - sum((obs/sum(obs))^2) - ( (sum((obs/sum(obs))^2)^2) + (sum((obs/sum(obs))^4)) )
	vecPIC[L] <- PIC 
	
	# Inbreeding Coefficient
	numer1 <- numer2 <- denom1 <- numeric()
	for(i in rownames(obs)) {
		
		Puu <- obs[i,i]/n
		pu <- alFreq[i]
		n1 <- Puu - pu^2
		numer1 <- c(numer1, n1)
		n2 <- pu - Puu
		numer2 <- c(numer2, n2)
		d1 <- pu*(1-pu)
		denom1 <- c(denom1, d1)
	}
	denom2 <- numer1
	
	vecf[L] <- ( sum(numer1) + sum(numer2)/(2*n) )/( sum(denom1) - sum(denom2)/(2*n) )
	
	L <- L+1
}

# Naming vectors and list with the loci names
names(nind) <- names(vech) <- names(vech.exp) <- names(vecpM) <- names(vecpD) <- names(vecpE) <- names(vecTPI) <- names(vecPIC) <- names(vecf) <- lnames

# Creating the allele frequency table
# Table 1
for (i in seq_along(dfall)) names(dfall[[i]]) <- c('alleles', paste0('freq.', i))
alfreqdf <- Reduce(function(x,y) merge(x,y, by = 'alleles', all = TRUE), dfall)
names(alfreqdf) <- c("alleles", lnames)
alfreqdf[is.na(alfreqdf)] <- ""

# Creating a data frame with the statistic of forensic relevance
# This is the bottom part of the Table 1
pars <- rbind(nind, round(vech,6), round(vech.exp,6), round(vecf,6), round(vecpM, 6), round(vecpD, 6), round(vecpE,6), round(vecTPI,6), round(vecPIC, 6))
rownames(pars) <- NULL
parsnames <- c("N", "het.obs", "het.exp", "f", "p.M", "p.dis", "p.exc", "TPI", "PIC")
pars <- data.frame(alleles=parsnames, pars)

tab1 <- rbind(alfreqdf, pars)

# Saving the table to the working directory
write.table(tab1, "table1_1.csv", sep=",", quote=FALSE, row.names=FALSE)

# Formating the data frame to DNAtools format
pop <- data.frame(id=rownames(pop), pop)
rownames(pop) <- NULL

# Newer version of DNAtools doesn't accept missing values.
# Code the missing alleles by e.g. 999 or another numeric value not in the set of alleles
pop[is.na(pop)] <- 999

# DNAtools analysis of match
# threads=0 can be used on linux and mac only.
require(DNAtools)
res <- dbCompare(pop, hit=12, threads=0, trace=FALSE)

# Save the results
mat <- as.data.frame(res$m)
mat <- data.frame("partial/match"=rownames(mat), mat)
colnames(mat) <- c("match/partial", 0:nloci)
write.table(mat, "table2.csv", sep=",", row.names=FALSE)

pop[pop==999] <- NA
# missing data plot
missing <- numeric()
for(i in seq(2,ncol(pop),2)) {
	
	m <- dim(pop[is.na(pop[,i]) | is.na(pop[,i+1]), ])[1] / dim(pop[,c(i,i+1)])[1]
	missing <- c(missing, m)
}
names(missing) <- lnames

# Saving the plot as .png with 1200 dpi of resolution to the working directory
png("figure1.png", width=3.25, height=3.25, units="in", res=1200, pointsize=6, family="Arial")
par(mar = c(7,4,2,1))
barplot(missing*100, space=.1, ylim=c(0,100), axisnames=TRUE, col="grey20", cex.names=1.2, las=2, names.arg=names(missing), border=NA, axis.lty="solid", ylab="missing data percentage (%)")
dev.off()
