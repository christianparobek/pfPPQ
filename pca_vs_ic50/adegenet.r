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

## Define on location
pop(gl) <- pop.definer(indNames(gl))

## Define based on PPQ and MQ IC50s
id <- as.character(data$WGS_ID[data$WGS_ID %in% indNames(gl)])
ppq <- data$ppq_ic50[data$WGS_ID %in% indNames(gl)]
mq <- data$mq_ic50[data$WGS_ID %in% indNames(gl)]
ord <- match(indNames(gl), id) # got to match the order of JonP's metadata to the VCF
pop(gl) <-
  as.numeric(ppq[ord] > quantile(ppq[ord], na.rm=TRUE)[4]) + 
  as.numeric(mq[ord] > quantile(mq[ord], na.rm=TRUE)[4])*2


#####################################################
##################### PLOT PCA ######################
#####################################################

scores <- pca$scores
scores <- jitter(pca$scores, factor=700)

plot(scores[,2] ~ scores[,1],
     type = 'n',
     axes = FALSE,
     xlab=paste("PC", 1, " - ", round(pca$eig[1]/sum(pca$eig)*100), "% of the Variance", sep = ""),
     ylab=paste("PC", 2, " - ", round(pca$eig[2]/sum(pca$eig)*100), "% of the Variance", sep = ""),
     )
points(scores[,2][is.na(pop(gl))] ~ scores[,1][is.na(pop(gl))]) # plot any with undetermined IC50
points(scores[,2][pop(gl) == 0] ~ scores[,1][pop(gl) == 0]) # plot any with undetermined IC50
points(scores[,2][pop(gl) == 1] ~ scores[,1][pop(gl) == 1], col = "red", pch=19) # plot high PPQ IC50
points(scores[,2][pop(gl) == 2] ~ scores[,1][pop(gl) == 2], col = "blue", pch=19) # plot high MQ IC50
points(scores[,2][pop(gl) == 3] ~ scores[,1][pop(gl) == 3], col = "green", pch=19) # plot high MQ IC50
axis(1)
axis(2)


#####################################################
############# CHOOSE TITLE AND LEGEND ###############
#####################################################

legend(-11, -4, legend = c("PPQ IC50 UQ", "MQ IC50 UQ", "PPQ & MQ IC50 UQ", "Missing Data or non-UQ"), col = c("red", "blue", "green", "black"), pch=c(19,19,19,1), cex=1)
title(substitute(paste("PPQ and MQ Resistance vs. Genetic Background")), line = -0.5, cex.main=1.2)
### OR ###
title(substitute(paste(italic('P. falciparum'), " PCA vs. PPQ IC50" )), line = -0.5, cex.main=1.5)
legend(-15, -5, legend = c("IC50 top 25%", "IC50 bottom 75%"), col = c("red", "black"), pch=19, bty="n", cex=1)
### OR ###
title(substitute(paste("Cambodian ", italic('P. falciparum'), " Isolates by Province" )), line = -0.5, cex.main=1.2)
legend(-12, -5, legend = c("Battambang", "Oddar Meanchey", "Kampot"), col = c("red", "blue", "green"), pch=19, bty="n", cex=1.2)
### OR ###
title(substitute(paste("PPQ and MQ Resistance vs. Genetic Background")), line = -0.5, cex.main=1.2)
legend(-12, -5, legend = c("PPQ Resistant", "MQ Resistant", "PPQ + MQ Resistant"), col = c("red", "blue", "green"), pch=19, bty="n", cex=1.2)

#####################################################
################## COLOR GATING #####################
#####################################################

cp1 <- as.numeric((pca$scores[,1] > 20) & (pca$scores[,2] > -5))
cp2 <- as.numeric((pca$scores[,1] > -5) & (pca$scores[,1] < 10) & 
                    (pca$scores[,2] > -5) & (pca$scores[,2] < 10))*2
cp3 <- as.numeric((pca$scores[,1] < -5) & (pca$scores[,2] < 0))*3
cp4 <- as.numeric((pca$scores[,1] < -5) & (pca$scores[,2] > 0))*4


#####################################################
############# ggplot TO MAKE SPECTRUM ###############
#####################################################

ppq200 <- ppq
ppq200[ppq200 > 200] <- 200

library(ggplot2)
library(gridExtra)

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

