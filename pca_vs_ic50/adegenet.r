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
#library(stringr)
library(pegas)
library(ggplot2)
library(gridExtra)

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


#####################################################
#################### INPUT DATA #####################
#####################################################

## Read in VCF
gl <- genlight.maker("../variants/slimmed.recode.vcf") # make genlight
indNames(gl) <- str_extract(indNames(gl), "[A-Z]+...") # fix the names

## Read in metadata
data <- read.table("PPQ_res_Aug20_v2.txt", sep="\t", header=TRUE)


#####################################################
################## CALCULATE PCA ####################
#####################################################

## Calculate PCA
pca <- glPca(gl)


#####################################################
################# ASSIGN GROUPINGS ##################
#####################################################

## Define based on PPQ and MQ IC50s
id <- as.character(data$WGS_ID[data$WGS_ID %in% indNames(gl)])
ppq <- data$ppq_ic50[data$WGS_ID %in% indNames(gl)]
mq <- data$mq_ic50[data$WGS_ID %in% indNames(gl)]
ord <- match(indNames(gl), id) # got to match the order of JonP's metadata to the VCF
pop(gl) <-
  as.numeric(ppq[ord] > quantile(ppq[ord], na.rm=TRUE)[4]) + 
  as.numeric(mq[ord] > quantile(mq[ord], na.rm=TRUE)[4])*2


#####################################################
################## COLOR GATING #####################
#####################################################

cp1 <- as.numeric((pca$scores[,1] > 20) & (pca$scores[,2] > -5))
cp2 <- as.numeric((pca$scores[,1] > -5) & (pca$scores[,1] < 10) & 
                    (pca$scores[,2] > -5) & (pca$scores[,2] < 10))*2
cp3 <- as.numeric((pca$scores[,1] < -5) & (pca$scores[,2] < 0))*3
cp4 <- as.numeric((pca$scores[,1] < -5) & (pca$scores[,2] > 0))*4

cols <- cp1 + cp2 + cp3 + cp4

#####################################################
############# ggplot TO MAKE SPECTRUM ###############
#####################################################

ppq200 <- ppq
ppq200[ppq200 > 200] <- 200

df <- as.data.frame(jitter(pca$scores, factor=700))

mqplot <- ggplot(df, aes(x=PC1, y=PC2)) + 
  geom_point(aes(size = mq[ord], color = factor(cols)), alpha=0.6) + 
  theme_bw() +
  ggtitle("Mefloquine IC50") +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    legend.position = "none") +
  scale_size_continuous(range = c(2,9)) +
  labs(
    x = paste("PC", 1, " - ", round(pca$eig[1]/sum(pca$eig)*100), "% of the Variance", sep = ""),
    y = paste("PC", 2, " - ", round(pca$eig[2]/sum(pca$eig)*100), "% of the Variance", sep = "")
  )

ppqplot <- ggplot(df, aes(x=PC1, y=PC2)) + 
  geom_point(aes(size = ppq200[ord], color = factor(cols)), alpha=0.6) + 
  theme_bw() +
  ggtitle("Piperaquine IC50") +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    legend.position = "none") +
  scale_size_continuous(range = c(2,9)) +
  labs(
    x = paste("PC", 1, " - ", round(pca$eig[1]/sum(pca$eig)*100), "% of the Variance", sep = ""),
    y = paste("PC", 2, " - ", round(pca$eig[2]/sum(pca$eig)*100), "% of the Variance", sep = "")
  )


grid.arrange(ppqplot, mqplot, ncol=2)

