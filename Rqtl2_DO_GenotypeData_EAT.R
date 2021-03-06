########## RQTL2 Analysis: DO Mice ##########
################################################################################
# R/QTL2 Genotype Data.
# Data provided by Dan Gatti.
# Anji Trujillo
# etrujillo2@wisc.edu
# July 10, 2017
################################################################################

##############################
# Load and install packages. #
##############################

install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen"))
install_github("rqtl/qtl2")
library(devtools)
library(qtl2) # Loads qtl2geno, qtl2scan & qtl2plot.
library(qtl2convert)
library(RSQLite)
library(dplyr)
library(qtl)
library(qtl2scan)


#####################
# Load in the data. #
#####################

# Files provided by Dan Gatti
setwd("C:/Users/etrujillo/Desktop/DOProjectFolder")
load('Attie_DO_genoprobs.Rdata') # Rdata file provided by Dan Gatti
load('GM_snps.Rdata') # Rdata file provided by Dan Gatti

# Output files from QTL analysis; this is specific for Liver Metabolites
load('qtl_LiverMetabolites_20170824.rData') # Qtl output for Liver Data
#load('qtl_LiverMetabolites_20170824_findpeaks.rData') #find peaks for Liver Data

#Required for qtl parameters
load('kinshipMatrix_DO_DanGattie_20170824.rData') # kinship Matrix generated
load('genoprobs_subsetmarkers_DO_DanGattie_QTl2convert_20170824.rData') # qtl2-style data structure with correct number of markers (120789)
load('genoProbs_withcorrectMouseandMarkers_DO_20170824.rData') # overwrite genoprobs with the correct number of mice and markers

# Phenotype file and covariate file
Liver_Metabolites <- read.csv(file = "RModified/21Aug2017DOLiverMetabolites_ComBatNormalized.csv") # Phenotype data
covar <- read.csv(file = "DOWave1through4_Covariates.csv")

#####################
# Explore the data. #
#####################

str(probs)
str(GM_snps)

rownames(probs) # print mice IDs
row.names(GM_snps) #142,259 Markers

row.names(Liver_Metabolites) <- Liver_Metabolites$Mouse.ID #set row names for phenotype file
Liver_Metabolites <- Liver_Metabolites[,-c(2:6)] # remove covariate information from phenotype file, from Normalization step

# Match the mice in phenotype file and covariate file
LiverMetabolites_Mice <- row.names(Liver_Metabolites) # save mouse id's from phenotype data
subSet.covar <- covar[covar$Mouse.ID %in% LiverMetabolites_Mice,] #120789 obs. of 16 variables

#Psuedo wave of mice with AA prefix
length(grep("AA", row.names(probs))) # 99 mice with AA prefix
length(grep("f", row.names(probs))) # 49 are females
length(grep("m", row.names(probs))) # 50 are females

#####################
# Data Preparation. #
#####################
# access attributes from probs (loaded earlier)
# create new object with dimnames
# extract and save markers

row.names(probs) <- correctName.probs <- c(substr(row.names(probs)[1:99],4,8),row.names(probs)[100:479]) # remove appended AA. and sex ID

attributes(probs) 
x <- dimnames(probs) 
markers <- x[[3]] # 120789 markers used
mice <- x[[1]] # 479 mice

###################
# Subset Markers. #
###################
# Subset markers from total marker count 143259 to the the correct number 120789
# Create a genetic map with marker, chromosome and position

subSet.GM_snps <- GM_snps[GM_snps$marker %in% markers,] #120789 obs. of 16 variables

DOgmap <- subSet.GM_snps[,1:3] # subset the marker, chromosome, and position (bp)

#################
# Convert data. #
#################
# Convert DOQTL-style allele probabilities to qtl2-style data structure.
# Calculate genetic similarity among individuals (kinship matrix) from conditional genotype probabilities
# pos_column = "pos" will provide QTL positions in Mb
# loco provides kinship matrix leaving out one chromosome at a time

genoprobs = qtl2convert::probs_doqtl_to_qtl2(probs = probs, map = subSet.GM_snps,  
              chr_column = "chr", pos_column = "pos", marker_column = "marker")
save(genoprobs, file =  "genoprobs_subsetmarkers_DO_DanGattie_QTl2convert_20170824.rData") # save genoprobs as rData object

K = calc_kinship(probs = genoprobs, type = "loco")
save(K, file =  "kinshipMatrix_DO_DanGattie_20170824.rData") # save K as rData object

################
# Subset Mice. #
################
# Subset mice to match the mice where phenotype and genotype available
# Use for loop to subset the mice only present in phenotype data
# array indexing [,,]


MouseNames <- row.names(genoprobs[[1]]) # save mouse id's from genoprob data
rowsKeep <- which(MouseNames %in% Liver_Metabolites$Mouse.ID) # identify which mice are present in both phenotype data and genoprob data

overWriteGenoProbs <- genoprobs # create a new genoprobs object

  for(i in names(genoprobs))
    { 
      overWriteGenoProbs[[i]] <- overWriteGenoProbs[[i]][rowsKeep,,]
    }
save(overWriteGenoProbs, file =  "genoProbs_withcorrectMouseandMarkers_DO_20170824.rData") # save overWriteGenoProbs as rData object

######################
# Set up covariates. #
######################

addcovar = model.matrix(~sex + Batch + Wave, data = subSet.covar)[,-1] #create a matrix by expanding factors to a set of dummy variables
row.names(addcovar) <- subSet.covar$Mouse.ID

sex <- (subSet.covar$sex == "M")*1
names(sex) <- subSet.covar$Mouse.ID

#########
# Scan1 #
#########
# Perform a genome scan by Haley-Knott regression using scan1()
# scan1() takes genotype probabilites, matrix of phenotypes (here Live rLipids)
# need to add

#qtl = scan1(genoprobs = overWriteGenoProbs, pheno = Liver_Metabolites[,-c(1)] , kinship = K, cores = 1)
#save(qtl, file =  "qtl_LiverLipid_20170824.rData")
#str(qtl)

qtl = scan1(genoprobs = overWriteGenoProbs, pheno = Liver_Metabolites[,-c(1)],
            kinship = K, addcovar = addcovar, cores = 1) # add covariate
save(qtl, file =  "qtl_LiverMetabolites_20170828.rData")
################
# Plotting QTL #
################
# Create map. Markers with position (bp)

DOgmapVector <- DOgmap[,3] # save position of markers into own vector
names(DOgmapVector) <- DOgmap[,1] # set names of object as the marker IDs 
map = split(DOgmapVector, DOgmap$chr) #
map = map[order(as.numeric(names(map)))]

#names(DOgmap) <- rownames(DOgmap$marker) # include individual IDs as names
# pdf("qtl_LiverLipids_plots.pdf", width = 14)
# plot_scan1(qtl, map, lodcolumn = 1, main = colnames(qtl)[1])
# dev.off()
# pdf("qtl_LiverLipids_plot2.pdf", width = 14)
# plot_scan1(qtl, map, lodcolumn = 2, main = colnames(qtl)[2])
# dev.off()

# Overlaying plots for all lipids in Liver_Lipid file
pdf("Figures/qtl_LiverMetabolites_plotloop_20170824.pdf", width = 14)
plot(qtl, map, lodcolumn = 1, main = "QTL for Liver Lipids")
  
for(i in 2:ncol(qtl))
    {
      plot(qtl, map, lodcolumn = i, col = i, add = TRUE)
    }
    dev.off()

# Overlaying plots  for lipids 1:i
pdf("Figures/qtl_LiverMetabolites_plot1through20_20170824.pdf", width = 14)
  plot(qtl, map, lodcolumn = 1)
    for(i in 2:5)
    {
      plot(qtl, map, lodcolumn = i, col = i, add = TRUE)
      #legend("topleft", legend = qtl , col=palette, names(palette), ncol=2, lwd=2, bg="gray95")
    }
dev.off()

pdf("qtl_LiverLipids_plot_TG.56.6.pdf", width = 14)
plot(qtl, map, lodcolumn = 366, col = "blue")
legend("topleft", col="blue", lwd=2, bg="gray95")
dev.off()

legend("topright", inset=.05, title="Number of Cylinders",
       legend =c("4","6","8"), fill=terrain.colors(3), horiz=TRUE)
legend("bottomleft", legend = levels(batch), pch = 16, col = 1:length(levels(batch)),
       x.intersp = 0.7, y.intersp = 0.7)
#####################
# Finding LOD Peaks #
#####################
# Function find_peaks() will identify a set of LOD peaks that exceed some threshold

findPeaks <- find_peaks(qtl, map, threshold=4, drop=1.5)
save(findpeaks, file =  "qtl_LiverMetabolites_20170824_findpeaks.rData")

subSetLiver_Lip <- Liver_Lip[, 2] #subset only one phenotype
names(subSetLiver_Lip)<- LiverLip_Mice # append names of mice to each measurement

chr = 15
qtl.coef = scan1coef(genoprobs = overWriteGenoProbs[,chr], subSetLiver_Lip, kinship = K[[chr]])

# Plot only coefficients.
pdf("FounderStrainEffects_chrom15.pdf", width = 14)
plot(x = qtl.coef, map = map[chr], columns = 1:8, col = CCcolors,
     main = colnames(qtl)[366])
legend("bottomleft", col=CCcolors, names(CCcolors), ncol=2, lwd=2, bg="gray95")
dev.off()

# Plot coefficients with LOD cruve.
pdf("CoefficieintWithLODCurve.pdf", width = 14)
plot(x = qtl.coef, map = map[chr], columns = 1:8, col = CCcolors,
     scan1_output = qtl, main = colnames(qtl)[1])
dev.off()

# Calculate BLUP coefficients on Chr
qtl.blup = scan1blup(genoprobs = overWriteGenoProbs[,chr], pheno = subSetLiver_Lip,
                     kinship = K[[chr]])
