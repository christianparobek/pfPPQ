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
library(abind)
library(ggplot2)
library(reshape)
library(grid)
library(gridExtra)
library(strataG)


#####################################################
################# Define Functions ##################
#####################################################

## Function to create genlight from VCF.
## The AAKM contigs can screw with this
## so `grep -v "=AAKM" infile.vcf > outfile.vcf`
genlight.maker <- function(infile) {
  loci <- read.vcf(infile, from = 1, to = 1000000) # reads first million sites
  genlight <- new("genlight", loci) # convert data frame into genlight object
  ploidy(genlight) <- as.integer(1) # add back population information
  return(genlight)
}

## Function to assign samples to pops
## based on their CP group number
pop.definer <- function(ind_names) {
  library(stringr)
  cp1 <- as.numeric(str_detect(ind_names, "OM352|SN003|SN019|SN032|SN043|SN060|SN066|SN072|SN082|SN083|SN093|SN099|KP054|KP065|BB084|BB059|BB080"))*1 # assign KP pop number
  cp2 <- as.numeric(str_detect(ind_names, "KP001|KP004|BB085|KP059|KP073|SN076|SN079|SN091|SN103|SN109|SN078|KP027|KP030|KP062|SN064|SN097|SN107|SN111"))*2# assign BB pop number
  cp3 <- as.numeric(str_detect(ind_names, "BB052|BB082|BB068|BB069|SN022|SN035|SN038|SN039|SN044|SN048|SN052|SN057|SN058|SN061|SN084|SN085|SN092|SN095|SN105|SN117"))*3 # assign OM pop number
  cp4 <- as.numeric(str_detect(ind_names, "SN015|SN016|SN030|SN031|SN042|SN046|SN047|SN062|SN063|SN071|SN074|SN075|SN077|SN081|SN086|SN096|SN098|SN101|SN102|SN106|SN108|SN114|SN116"))*4 # assign SN pop number
  srr <- as.numeric(str_detect(ind_names, "SRR|ERR|TB"))*5 # assign SRR and ERR pop number
  pops <- cp1 + cp2 + cp3 + cp4 + srr
  return(pops)
}

## Function to plot PCAs
pca.plotter <- function(pca, pops, x, y) {
  plot(pca[,y] ~ pca[,x], 
       col=pops, 
       pch=19, 
       axes=FALSE, 
       xlab=paste("PC", x, sep = ""),
       ylab=paste("PC", y, sep = ""),
  )
  axis(1)
  axis(2, las = 2)
}


###############################################
########## PCA - READ & PROCESS DATA ##########
###############################################

## read in PF data
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/cp-kh/")
pf_gl <- genlight.maker("joint.vcf") # make genlight
pf_pops <- pop.definer(indNames(pf_gl)) # define pops
pf_pca <- glPca(pf_gl, nf = 8) # calculate PCA
pf_pca_jit <- as.data.frame(jitter(pf_pca$scores, factor=100))
# and make a copy with some jitter

# remove the 65 that are from the Gambia, which seem to have a PC1 of >15
sansGambia <- pf_pca_jit[!pf_pca_jit$PC1 > 17,]
sansGambiaCol <- pf_pops[!pf_pca_jit$PC1 > 17]

# replace colors
sansGambiaCol[sansGambiaCol == 1] <- "deepskyblue"
sansGambiaCol[sansGambiaCol == 2] <- "brown1"
sansGambiaCol[sansGambiaCol == 3] <- "darkolivegreen3"
sansGambiaCol[sansGambiaCol == 4] <- "darkgoldenrod1"
sansGambiaCol[sansGambiaCol == 5] <- "gray"


######################################################
################# PLOT PCA & BOXPLOT #################
######################################################

svg("cp_kh.svg", width = 11, height = 5.5)
par(mfrow = c(2,2), mar = c(5,5.5,4,2))

# plot the boxplot
boxplot(boots[9,2,], boots[5,2,], boots[6,2,], boots[7,2,], boots[8,2,], axes = FALSE, ylab = "CP2 vs. KH Group\n(Euclidean Distance)", xlab = "KH Groups")
axis(1, at = 1:5, labels = c("KHA", "KH1", "KH2", "KH3", "KH4"))
axis(2, las = 2)
mtext("B", 2, las = 2, cex = 2, at = 5.6, line = 3)

# plot the pca
pca.plotter(sansGambia, sansGambiaCol, 1, 2)
text(c(-12, -18, -5, 2), c(21, -8, 2, -7), labels = c("CP4", "CP3", "CP2", "CP1"), col = c("darkgoldenrod1", "darkolivegreen3", "brown1", "deepskyblue"))
mtext("A", 2, las = 2, cex = 2, at = 21, line = 3)

# plot the pca
pca.plotter(sansGambia, sansGambiaCol, 1, 3)
text(c(-12, -18, -5, 2), c(21, -8, 2, -7), labels = c("CP4", "CP3", "CP2", "CP1"), col = c("darkgoldenrod1", "darkolivegreen3", "brown1", "deepskyblue"))
mtext("A", 2, las = 2, cex = 2, at = 21, line = 3)

# plot the pca
pca.plotter(sansGambia, sansGambiaCol, 2, 3)
text(c(-12, -18, -5, 2), c(21, -8, 2, -7), labels = c("CP4", "CP3", "CP2", "CP1"), col = c("darkgoldenrod1", "darkolivegreen3", "brown1", "deepskyblue"))
mtext("A", 2, las = 2, cex = 2, at = 21, line = 3)

dev.off()


######################################################
################ TRY IT ggplot STYLE #################
######################################################

## define theme
theme <- theme_bw() + 
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.position = "none",
        plot.title = element_text(size = rel(2)))

#####
## Define the palette order
#####
palette <- c("brown1", "darkgoldenrod1", "darkolivegreen3", "deepskyblue", "grey")

#####
## Define a useful plotting function
#####
plot.pca <- function(df, plotname, PCx, PCy) {
  pic <- ggplot(df, aes_string(x=colnames(df)[PCx], y=colnames(df)[PCy])) + 
    geom_point(aes(size = 2, color = sansGambiaCol), alpha=0.9) +
    theme +
    scale_color_manual(values=palette) +
     labs(x = paste("PC", PCx, " - ", round(pf_pca$eig[PCx]/sum(pf_pca$eig)*100), "% variance", sep = ""),
          y = paste("PC", PCy, " - ", round(pf_pca$eig[PCy]/sum(pf_pca$eig)*100), "% variance", sep = "")) +
    theme(plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"))
  return(pic)
}

pf_1_2 <- plot.pca(sansGambia, "", 1, 2)
pf_1_3 <- plot.pca(sansGambia, "", 1, 3)
pf_2_3 <- plot.pca(sansGambia, "", 2, 3)

## Prepare the data for the boxplot
## Remember, you'll need to run the code in the boxplotter.r script
com <- as.data.frame(t(boots[5:9,2,]))
com1 <- com[c(5,1,2,3,4)]
names(com1) <- c("KHA", "KH1", "KH2", "KH3", "KH4")

## Plot the boxplot
box <- ggplot(data = melt(com1), aes(variable, value)) + 
  geom_boxplot(aes(fill = "brown1")) +
  theme +
  labs(x = "KH Groups", y = "CP2 vs. KH Group\n(Euclidean Distance)") +
  theme(axis.title.y=element_text(vjust=1)) +
  theme(plot.margin=unit(c(1.5,1.5,1.5,1.2), "lines"))

## Print the figure
svg("Figure SX - CP_KH.svg", width = 11, height = 10.5)
  grid.arrange(pf_1_2, pf_1_3, pf_2_3, box, ncol=2)
  grid.text("A", x = unit(0.03, "npc"), y = unit(0.98, "npc"), gp=gpar(fontsize=30))
  grid.text("B", x = unit(0.54, "npc"), y = unit(0.98, "npc"), gp=gpar(fontsize=30))
  grid.text("C", x = unit(0.03, "npc"), y = unit(0.49, "npc"), gp=gpar(fontsize=30))
  grid.text("D", x = unit(0.54, "npc"), y = unit(0.49, "npc"), gp=gpar(fontsize=30))
dev.off()


######################################################
########## K-means CLUSTERING TO ID GROUPS ###########
######################################################

#####
## Identify clusters using K-means
#####

## These settings appropriately group the CP3 and CP4 groups with themselves
## And groups the nearby ERR samples in there too.
grp <- find.clusters(pf_gl, n.pca = 4, max.n.clust = 10, choose.n.clust = FALSE, criterion = "diffNgroup")
dapc <- dapc(pf_gl, grp$grp, n.pca = 4, n.da = 20)
scatter(dapc)

## Then cycle through the group memberships and figure out which correspond to CP3 and CP4

pf_gl$ind.names[dapc$grp == 4] # CP4
pf_gl$ind.names[dapc$grp == 6] # CP3

## Now go do some calculations in popgenome_pi_calc.r

