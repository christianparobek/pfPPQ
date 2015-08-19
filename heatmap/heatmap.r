## An ape script to analyze Jess / Jon Parr's 
## Modified from the Acineto pheatmap.r
## Begun 18 August, 2015
## Example: http://stackoverflow.com/questions/15153202/how-to-create-a-heatmap-with-a-fixed-external-hierarchical-cluster

################################################
########## LOAD THE LIBRARIES WE NEED ##########
################################################

#library(ape)
#library(phangorn)
library(adegenet)
library(pegas)

################################################
######### DEFINE SOME USEFUL FUNCTIONS #########
################################################

## Function to create genlight from VCF.
genlight.maker <- function(infile) {
  loci <- read.vcf(infile)
  genlight <- new("genlight", loci) # convert data frame into genlight object
  ploidy(genlight) <- 1 # add back population information
  return(genlight)
}


diff.paster <- function(x){
  paste(rownames(ord_mat)[x], " (", snps[x], ")", sep="")
} # add the num of snp diffs to the names

################################################
################# READ IN DATA #################
################################################

gl <- genlight.maker("/run/user/1001/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/pf/variants/our_goods_UG.pass.vcf") # make genlight


matrix <- as.matrix(gl)
heatmap(matrix)


mat <- matrix(rnorm(80*8),nrow=80) #Some random data to plot
heatmap(mat, Rowv=NA, Colv=NA, margins=c(5,10))

################################################
###### COUNT DIFFS BETWEEN A03 AND OTHERS ######
################################################

good42_matrix <- as.matrix(good42_gl) # convert gl to matrix
ord_mat <- good42_matrix[ultra$tip,] # ord matrix by phylo tree

snps <- apply(ord_mat, 1, a03.comparer) # calculate nuc diffs bt A03 and others
new_names <- unlist(lapply(1:nrow(ord_mat), diff.paster)) # make list of updated names
rownames(ord_mat) <- new_names # assign new names to ord_mat

################################################
############### MAKE THE HEATMAP ###############
################################################

heatmap(ord_mat, Rowv=dend, Colv=NA, labCol="", col=topo.colors(4), margins=c(1,20)) # draw pheatmap
text(0.29, -8, "35551 SNVs")

# ML example: http://cran.r-project.org/web/packages/phangorn/vignettes/Trees.pdf
# Need to run the treeNJ phylo object through the ML Process
# Use the above link for an example of how to do this
fit <- pml(treeNJ, data=as.phyDat(data))




library(reshape)
library(ggplot2)

mat <- matrix(rnorm(80*2),nrow=80) #Some random data to plot
y <- melt(mat)
y <- melt(jp_subset_ordered[,1:2])
p <- ggplot(y, aes(y=X1,x=X2))
p + geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="white", high="steelblue") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
    )

q <- p + geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="white", high="red") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )











jp_subset_ordered





#########################################
############## IMPORT DATA ##############
#########################################

jp_data <- read.table("MRP2_Aug18.txt", header=TRUE, sep = "\t")


#########################################
############### CLEAN DATA ##############
#########################################

subset(jp_data, )




newdata <- subset(mydata, age >= 20 | age < 10,
                  select=c(ID, Weight))




jp_subset <- cbind(jp_data$ppq_ic50, jp_data$mq_ic50, jp_data$pfmdr1_copynum, jp_data$MIIg_present)
colnames(jp_subset) <- c("ppq_ic50", "mq_ic50", "pfmdr1_copynum", "MIIg_present")
rownames(jp_subset) <- jp_data$Screening_ID
heatmap(jp_subset)



jp_subset_ordered <- jp_subset[order(jp_subset[,1]),]

heatmap(jp_subset_ordered, Rowv=NA, Colv=NA)




jp_subset_ordered[1:80,]


library(pheatmap)
pheatmap(log(jp_subset_ordered+1), cluster_rows=FALSE, cluster_cols=FALSE)


jp_subset_ordered <- jp_subset_ordered[-103,]

log(jp_subset+1)


jp_subset_ordered[,1] <- ifelse(jp_subset_ordered[,1] > 1000, NA)
jp_
