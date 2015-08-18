# For genetic analysis of my WGS Plasmodium population(s)
# Started 10 Dec 2014
# Updated 01 May 2015
# Basics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf
# Genomics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf
# Extra Commands: http://www.inside-r.org/packages/cran/adegenet/docs/.rmspaces


###################################################
################# Load Libraries ##################
###################################################

library(adegenet)
library(stringr)
library(pegas)


#####################################################
################# Define Functions ##################
#####################################################

## Function to create genlight from VCF.
genlight.maker <- function(infile) {
  loci <- read.vcf(infile)
  genlight <- new("genlight", loci) # convert data frame into genlight object
  ploidy(genlight) <- 1 # add back population information
  return(genlight)
}

## Function to assign samples to pops
## based on a list of their names
pop.definer <- function(ind_names) {
  kp <- as.numeric(str_detect(ind_names, "BB"))*1 # assign KP pop number
  bb <- as.numeric(str_detect(ind_names, "KP"))*2 # assign BB pop number
  om <- as.numeric(str_detect(ind_names, "OM"))*3 # assign OM pop number
  sn <- as.numeric(str_detect(ind_names, "SN"))*3 # assign SN pop number
  tb <- as.numeric(str_detect(ind_names, "TB"))*3 # assign TB pop number
  srr <- as.numeric(str_detect(ind_names, "SRR"))*4 # assign SRR pop number
  err <- as.numeric(str_detect(ind_names, "ERR"))*4 # assign ERR pop number
  pops <- kp + bb + om + sn + tb + srr + err
  return(pops)
}

# Function to mark hi IC50 samples, given a vector of samples and of IC50s
ic50.marker <- function(ind_names, hi_ic50s) {
  ic <- as.numeric(str_detect(ind_names, paste(hi_ic50s, sep = "", collapse = "|")))
  return(ic)
}

## Function to plot eigenplot
eig.plotter <- function(pca) {
  barplot(pca$eig, xlab = "", ylab = "Variance")
}

## Function to plot PCA
pca.plotter <- function(pca, pops, x, y) {
  plot(jitter(pca$scores[,y], factor=700) ~ jitter(pca$scores[,x], factor=700), 
     col=pops, 
     pch=19, 
     axes=FALSE, 
     xlab=paste("PC", x, " - ", round(pca$eig[x]/sum(pca$eig)*100), "% of the Variance", sep = ""),
     ylab=paste("PC", y, " - ", round(pca$eig[y]/sum(pca$eig)*100), "% of the Variance", sep = ""),
  )
  axis(1)
  axis(2)
}


##############################################################
############## Plasmodium falciparum Cambodia ################
##############################################################

### PREPARE DATA FOR ANALYSIS ###
gl <- genlight.maker("/run/user/1001/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/pf/variants/our_goods_UG.pass.vcf") # make genlight
hi <- read.table("hi_ic50s.txt", header = FALSE)
pop(gl) <- pop.definer(indNames(gl)) # define pops OR
pop(gl) <- as.factor(str_extract(indNames(gl), "[A-Z]+")) # define pops as factor OR
pop(gl) <- ic50.marker(indNames(gl), hi$V1) # mark the high IC50s
pca <- glPca(gl) # calculate PCA




### PLOT EIGENVALUES ###
eig.plotter(pca)
title(substitute(paste("Cambodia ", italic('P. falciparum'), " Eigenvalues" )), line = 0.5, cex.main = 1.5)

### PLOT PCA PICTURE ###
pca.plotter(pca, pop(gl), 1, 2)

### CHOOSE A LEGEND AND TITLE ###
title(substitute(paste(italic('P. falciparum'), " PCA vs. PPQ IC50" )), line = -0.5, cex.main=1.5)
legend(-15, -5, legend = c("IC50 top 25%", "IC50 bottom 75%"), col = c("red", "black"), pch=19, bty="n", cex=1)
### OR ###
title(substitute(paste("Cambodian ", italic('P. falciparum'), " Isolates by Province" )), line = -0.5, cex.main=1.2)
legend(-12, -5, legend = c("Battambang", "Oddar Meanchey", "Kampot"), col = c("red", "blue", "green"), pch=19, bty="n", cex=1.2)
### OR ###
title(substitute(paste("PPQ and MQ Resistance vs. Genetic Background")), line = -0.5, cex.main=1.2)
legend(-12, -5, legend = c("PPQ Resistant", "MQ Resistant", "PPQ + MQ Resistant"), col = c("red", "blue", "green"), pch=19, bty="n", cex=1.2)

