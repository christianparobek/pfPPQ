## Started 25 September 2015
## Updated 28 June 2016
## To plot AFS for my collected and simulated data.
## Need to mod input files in the following ways:
##    1) Remove all lines with dashes
##    2) For mono/cp2 & pf/pv all comparisons, lengthen shorter one with tabbed zeros

#####################################
########### READ IN DATA ############
#####################################

cp1 <- read.table("afsPlots/data/cp1/dadi_output_afs/cp1.afs", header = TRUE)
cp2 <- read.table("afsPlots/data/cp2/dadi_output_afs/cp2.afs", header = TRUE)
cp3 <- read.table("afsPlots/data/cp3/dadi_output_afs/cp3.afs", header = TRUE)
cp4 <- read.table("afsPlots/data/cp4/dadi_output_afs/cp4.afs", header = TRUE)
slimmed84fix <- read.table("afsPlots/data/slimmed84fix/dadi_output_afs/slimmed84fix.afs", header = TRUE)


#####################################
######### PF & PV ALL DATA ##########
#####################################

plotter <- function(df, name) {
  barplot(t(df[,1:4]), beside = TRUE, col=c("#d7191c","#fdae61","#abd9e9","#2c7bb6"),
          ylab=name, xlab="", border = NA, axes = FALSE)
  axis(2, las = 1)
  axis(1, at = 1:nrow(df)*5-2, labels = 1:nrow(df), tick = FALSE, cex.axis = 0.5, line = -0.5)
  lines(1:36*5-2, df$ms[1:36], lwd=2, col = "black", lty=2) # showing sims as lines here
}

svg(filename = "subgroup_afs_sans3.svg", width = 8, height = 7.2)
par(mfrow = c(4, 1),
    oma = c(5,2,0,0),
    mar = c(0,4,2,0))
plotter(cp1, "CP1")
plotter(cp2, "CP2")
#plotter(cp3, "CP3")
plotter(cp4, "CP4")
plotter(slimmed84fix, "Entire Population")
mtext("Minor Allele Frequency", side = 1, line = 1.8, cex = 0.7)
legend(x = 160, y = 0.25, 
       legend = c("Synonymous", "Nonsynonymous", "Genic", "Intergenic", "Simulation"), 
       col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6", "black"), 
       lty=c(NA, NA, NA, NA, 2), 
       pch=c(15, 15, 15, 15, NA),
       pt.cex=c(1.8, 1.8, 1.8, 1.8, NA), lwd=2, bty = "n", cex=1)
dev.off()

