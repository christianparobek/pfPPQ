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


#####################################################
################# Define Functions ##################
#####################################################

## Function to create a genlight object from
## STRUCTURE file. VCF->STR using PGDSpider2
genlight.maker <- function(infile) {
  table <- read.table(infile, skip=1, na.strings="-9") # read in data, missing is "-9" in str format
  sorted <- table[order(table[,1]),] # sort
  inds <- sorted$V1 # grab the indiv names
  pops <- sorted$V2 # grab the pop names
  sorted <- sorted[-c(1,2)] # remove ind and pop columns from data frame
  genlight <- new("genlight", sorted) # convert data frame into genlight object
  indNames(genlight) <- inds # add back individual information
  ploidy(genlight) <- 1 # add back population information
  return(genlight)
}

## Function to assign samples to pops
## based on a list of their names
pop.definer <- function(ind_names) {
  library(stringr)
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
  library(stringr)
  ic <- as.numeric(str_detect(ind_names, paste(hi_ic50s, sep = "", collapse = "|")))
  return(ic)
}

## Function to record eigenplots
eig.plotter <- function(pca) {
  barplot(pca$eig, xlab = "", ylab = "Variance")
}

## Function to record PCAs
pca.plotter <- function(pca, pops, x, y) {
  plot(jitter(pca$scores[,y], factor=300) ~ jitter(pca$scores[,x], factor=300), 
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
pf_cam_gl <- genlight.maker("/run/user/1001/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/adegenet/our_goods_pf.pass.str") # make genlight
hi <- read.table("hi_ic50s.txt", header = FALSE)
pf_cam_pops <- pop.definer(indNames(pf_cam_gl)) # define pops OR
pf_cam_pops <- ic50.marker(indNames(pf_cam_gl), hi$V1) # mark the high IC50s
pf_cam_pca <- glPca(pf_cam_gl) # calculate PCA

### PLOT EIGENVALUES ###
eig.plotter(pf_cam_pca)
title(substitute(paste("Cambodia ", italic('P. falciparum'), " Eigenvalues" )), line = 0.5, cex.main = 1.5)

### PLOT PCA PICTURE ###
pca.plotter(pf_cam_pca, pf_cam_pops + 1, 1, 2)
legend(-80, -30, legend = c("IC50 top 25% ", "IC50 bottom 75%"), col = c("red", "black"), pch=19, bty="n", cex=1.5)
title(substitute(paste(italic('P. falciparum'), " PCA vs. PPQ IC50" )), line = -0.5, cex.main=1.5)






library(pegas)
loci <- read.vcf("/run/user/1001/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/pf/variants/our_goods_UG.pass.vcf", nloci=100)

