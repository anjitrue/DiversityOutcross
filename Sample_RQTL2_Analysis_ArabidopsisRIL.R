#### Sample RQTL2 Analysis: Arabidopsis recombinant inbred lines (RIL) ####
################################################################################
# qtl2 mapping example.
# Moore et al. (2013) Genetics 195:1077-1086
# Anji Trujillo
# etrujillo2@wisc.edu
# July 12, 2017
################################################################################

##############################
# Load and install packages. #
##############################

install.packages("qtl2", repos="http://rqtl.org/qtl2cran/bin/windows/contrib/3.4/") #install R/qtl2 via mini-CRAN at rqtl.org
install_github("rqtl/qtl2")
install.packages("qtl2", repos="http://rqtl.org/qtl2cran")

options(stringsAsFactors = F)

library(devtools)
library(qtl2) # Loads qtl2geno, qtl2scan & qtl2plot.
library(qtl2convert)
library(RSQLite)
library(dplyr)
library(qtl)

#####################
# Load in the data. #
#####################
# Data are in qtl2geno/extdata/grav2.zip

grav2 <- read_cross2( system.file("extdata", "grav2.zip", package="qtl2geno") )
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DOex/DOex.zip")
DOex <- read_cross2(file)

####################################
# Calculate genotype probabilities #
####################################
# First task in QTL analysis is to calculate conditional genotype probabilities
# given observed marker data, at each putative QTL position.
# Use calc_genoprob() in glt2geno package. 
# Result is returned as a list of 3-D arrays (one per chromosome)

iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2geno") )
str(iron) #with chromosome and marker

map <- insert_pseudomarkers(iron$gmap, step=1) # insert psuedomarkers between markers
pr <- calc_genoprob(iron, map, err=0.002) #calculate QTL genotype probabilites at each marker and psuedomarker


pr <- calc_genoprob(DOex, error_prob=0.002) # calculate genotype probabilities for DO
apr <- genoprob_to_alleleprob(pr) # convert to allele probabilities

############################
# Calculate kinship matrix #
############################
# By default genotype probabilites are converted to allel probabilities
# kinship matrix calculates the portion of shared allels
# To eliminate the effect of varying marker density accross the genome only use probabilites 
# along the grid of psudedomarker (defined by the step argument in insert_psuedomarkers())

kinship <- calc_kinship(pr, use_allele_probs = FALSE, omit_x = TRUE) 

grid <- calc_grid(iron$gmap, step=1) # determine the grid of pseudomarkers
pr_grid <- probs_to_grid(pr, grid) # determine probabilities for positions that are not on the grid
kinship_grid <- calc_kinship(pr_grid)

kinship_loco <- calc_kinship(pr, "loco") # for linearl mixed model genome scan
kinship_loco[[1]]

k <- calc_kinship(apr, "loco") # calculate kinship for for DO

##################################
# Covariates for the X chromosome#
##################################

Xcovar <- get_x_covar(iron)

sex <- (DOex$covar$Sex == "male")*1
names(sex) <- rownames(DOex$covar) # include individual IDs as names

#########
# Scan1 #
#########

out <- scan1(pr, iron$pheno, Xcovar=Xcovar)

out <- scan1(apr, DOex$pheno, k, sex)

#################
# Plot the data #
#################

par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(out, DOex$gmap)

DOex$gmap


plot(sug) # 
sug <- calc.genoprob(sug, step = 1) # insert the QTL genotype probabilites along with the step (density of mapping units)

out.em <-  scanone(sug) # performs a single-QTL genome scan
summary(out.em, threshold = 3) # return chromosomes with LOD scores greater than 3
plot(out.em) # plots LOD curves

out.hk <- scanone(sug, method = "hk") #genome scan via Haley Knott regression
plot(out.em, out.hk, col = c("blue", "red")) # plot out.hk
plot(out.hk - out.em, ylim = c(-0.3, 0.3), ylab = "LOD(HK) - LOD(EM)") # plot difference between two genome scans (hk - single)

sug <- sim.geno(sug, step = 1, n.draws = 64) # perform a genome scan by multiple imputations using sim.geno function, ex. 64 imputations
out.imp <- scanone(sug, method = "imp") #

plot(out.em, out.hk, out.imp, col = c("blue", "red", "green")) # plot the three curves
plot(out.em, out.hk, out.imp, col = c("blue", "red", "green"), chr = c(7,15)) # plot the three curves for chromosomes 7 and 15
plot(out.imp - out.em , out.hk - out.em, col = c("blue", "red", "green"), ylim = c(-1,1)) # plot difference between genome scans

operm <- scanone(sug, method = "hk", n.perm = 1000) #
plot(operm) #1000 genome wide scans with a maximum LOD Scores
summary(operm, perms = operm, alpha = 0.2)
summary(operm)


