forensicstr
===========

### R scripts to analyze forensic STR datasets. 

##### Rforensic.R consumes a csv file with the individuals' genotypes and produces 2 tables and a histogram.

Table 1 contains:

* Allele frequencies for each locus
* Number of individuals
* Observed Heterozygosity
* Expected Heterozygosity
* The moment estimate of the Inbreeding Coefficient
* Matching Probability
* Power of Exclusion
* Power of Discrimination
* Typical Paternity Index
* Polymorphic Information Content

Table 2:

table of complete and partial adventitious matches as computed by the DNAtools R package. If you run this part of the code, the program can take several minutes to finish, depending on the amount of individuals in the dataset.

Histogram:
Percentages of missing data across loci.

##### fconfint.R calculates confidence intervals for the inbreeding coefficient estimates.

##### see paper: http://dx.doi.org/10.1016/j.fsigen.2014.04.015
 
