## An R script to plot the output of Gubbins
## Started 15 July 2016

################################################
########## LOAD THE LIBRARIES WE NEED ##########
################################################

library(ape)
library(pegas)


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

gl <- genlight.maker("slimmed.recode.84fix_with3D7.vcf")
  # read in VCF

tree <- root(nj(dist(as.matrix(gl))), "Pf3D7", resolve.root = TRUE)
  # make tree
  # root tree

boots <- boot.phylo(tree, as.matrix(gl), function(x) nj(dist(x)), B = 100, block = 5)
  # bootstrap tree

svg("treeroot_something.svg", width = 10, height = 15)
plot.phylo(x=tree, typ="phy", show.tip=FALSE, x.lim = c(0, 100), show.node.label = TRUE)
dev.off()


write.tree(treeroot, file = "treeroot.nwk")



####################################
############ EXAMPLES ##############
####################################


### Are bootstrap values stable?
for (i in 1:5)
  print(boot.phylo(tr, woodmouse, f, quiet = TRUE))
### How many partitions in 100 random trees of 10 labels?...


###########################################
################# Add Color ###############
###########################################

cp1 <- c("OM352", "BB084", "SN003", "SN019", "SN032", "SN043", "SN060", "SN066", "SN072", "SN082", "SN083", "SN093", "SN099", "BB059", "KP054", "BB080", "KP065")
cp2 <- c("BB085", "KP001", "KP004", "KP027", "KP030", "KP059", "KP062", "KP073", "SN064", "SN076", "SN078", "SN079", "SN091", "SN097", "SN103", "SN107", "SN109", "SN111")
cp3 <- c("BB052", "BB082", "BB068", "BB069", "SN022", "SN035", "SN038", "SN039", "SN044", "SN048", "SN052", "SN057", "SN058", "SN061", "SN085", "SN092", "SN094", "SN095", "SN105", "SN117")
cp4 <- c("SN015", "SN016", "SN030", "SN031", "SN042", "SN046", "SN047", "SN062", "SN063", "SN071", "SN074", "SN075", "SN077", "SN081", "SN086", "SN096", "SN098", "SN101", "SN102", "SN106", "SN108", "SN114", "SN116")

colors <- tree$tip.label

colors[colors %in% cp1] <- "deepskyblue"
colors[colors %in% cp2] <- "brown1"
colors[colors %in% cp3] <- "darkolivegreen3"
colors[colors %in% cp4] <- "darkgoldenrod1"
colors[colors == "Pf3D7"] <- "black"


###########################################
################# MY TRY ##################
###########################################

func <- function(xx) root(nj(dist(xx)), "Pf3D7", resolve.root = TRUE)

bstrees <- boot.phylo(tree, as.matrix(gl), func, trees = TRUE)$trees
clad <- prop.clades(tree, bstrees, rooted = TRUE)

clad[clad < 90] <- NA

svg("treeroot_6.svg", width = 4, height = 7)
plot.phylo(tree, tip.color = colors, no.margin = TRUE, label.offset = 1.2, show.tip.label = TRUE, cex = 0.5)
#nodelabels(clad, bg = "white", adj = c(1.1, -0.5), font = 3, frame = "none", cex = 0.5)
nodelabels(clad, bg = "white", adj = c(1.1, -0.5), font = 3, frame = "none", cex = 0.4, col = "brown1")
tiplabels(pch=21, col="black", bg=colors, lwd=0.8, cex=1, adj = c(1.2, 0.5))
legend(8,13, legend = c("CP1", "CP2", "CP3", "CP4", "Pf3D7"), pch = 21, col = "black",
       pt.bg = c("deepskyblue", "brown1", "darkolivegreen3", "darkgoldenrod1", "black"), pt.cex = 1.5, cex = 0.7)

start <- 2
length <- 2
segments(2,1,4,1)
text(start+length*0.5, 0, 
     labels = signif(length * length(gl$loc.names) / 23292622, digits = 1), 
     cex = 0.5)

dev.off()

##########################
## NOTE ABOUT SCALE BAR ##
##########################
  ## add.scale.bar(cex=0.7, length = 1)
  ## adds an incorrect scale bar because assumes genome size is ~7K
  ## And returns a very high substitution rate
  ## So instead I took the length of the bar (eg. 2), representing 2 subst/site
  ## multiplied it by the number of sites gl length (~7K)
  ## Then divided it by the true genome length (~23M)




