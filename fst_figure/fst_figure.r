## Making a figure that calculates FST across genome
## Particularly coloring genes that are mutant in CP4
## but not in CP3

## Christian Parobek
## Started 19 April 2016

##########
### INITIALIZE LIBRARIES
##########

library(PopGenome)
library(stringr)


##########
### USEFUL FUNCTIONS
##########

## Break a GENOME object by gene and name it
get.pf.genes <- function(GENOME){
  genes <- splitting.data(GENOME, subsites="gene", whole.data = FALSE)
  ## Specify whole.data = FALSE because I don't want to concatenate regions
  ## Extracting gene IDs from the gff INFO field
  ## Need to use stringr's str_extract function to actually get the gene ID out of the big long string of INFO it returns
  id_str <- "PF3D7_[0-9, a-zA-Z, \\., -]+"
  chr01.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr01.gff", chr="Pf3D7_01_v3", extract.gene.names=TRUE), id_str)
  chr02.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr02.gff", chr="Pf3D7_02_v3", extract.gene.names=TRUE), id_str)
  chr03.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr03.gff", chr="Pf3D7_03_v3", extract.gene.names=TRUE), id_str)
  chr04.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr04.gff", chr="Pf3D7_04_v3", extract.gene.names=TRUE), id_str)
  chr05.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr05.gff", chr="Pf3D7_05_v3", extract.gene.names=TRUE), id_str)
  chr06.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr06.gff", chr="Pf3D7_06_v3", extract.gene.names=TRUE), id_str)
  chr07.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr07.gff", chr="Pf3D7_07_v3", extract.gene.names=TRUE), id_str)
  chr08.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr08.gff", chr="Pf3D7_08_v3", extract.gene.names=TRUE), id_str)
  chr09.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr09.gff", chr="Pf3D7_09_v3", extract.gene.names=TRUE), id_str)
  chr10.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr10.gff", chr="Pf3D7_10_v3", extract.gene.names=TRUE), id_str)
  chr11.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr11.gff", chr="Pf3D7_11_v3", extract.gene.names=TRUE), id_str)
  chr12.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr12.gff", chr="Pf3D7_12_v3", extract.gene.names=TRUE), id_str)
  chr13.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr13.gff", chr="Pf3D7_13_v3", extract.gene.names=TRUE), id_str)
  chr14.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr14.gff", chr="Pf3D7_14_v3", extract.gene.names=TRUE), id_str)
  ## Make one big list out of all the little lists of gene IDs
  ## Is there a more elegant way to do all this?
  chr.all.ids <- unlist(list(chr01.ids, chr02.ids, chr03.ids, chr04.ids, chr05.ids, chr06.ids, chr07.ids, chr08.ids, chr09.ids, chr10.ids, chr11.ids, chr12.ids, chr13.ids, chr14.ids))
  genes@region.names <- chr.all.ids
  return(genes)
}


##########
### READ IN DATA
##########

## read in genome data
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/pfPPQ/fst_figure/pf_whole_genome/")
pf <- readData("vcf/", format="VCF", gffpath = "gff/")

## read in groupings
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/pfPPQ/dadi/data/cp_groups/")
cp1 <- read.table("cp1.txt")
cp2 <- read.table("cp2.txt")
cp3 <- read.table("cp3.txt")
cp4 <- read.table("cp4.txt")
pfPops <- set.populations(pf, list(as.character(cp1$V1), as.character(cp2$V1), as.character(cp3$V1), as.character(cp4$V1)))

## read in metadata
## these are snps that are missense or nonsense in CP4
## but not in CP3
appendix <- read.table("appendixS2", header = TRUE, sep = "\t")
genedata <- read.table("genes_to_mark", header = TRUE, sep = "\t")


## To get gene names
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/pfPPQ/fst_figure/")       
pfPops_genenames <- get.pf.genes(pfPops)

##########
### CALCULATE FST
##########

## Split data into genes
pfPops_genes <- splitting.data(pfPops, subsites="gene", whole.data = FALSE)
pfPops_genes <- F_ST.stats(pfPops_genes)
# Specify whole.data = FALSE because I don't want to concatenate regions


##########
### HOUSEKEEPING
##########

fsts <- pfPops_genes@nuc.F_ST.pairwise[6,]

chrs <- str_extract(pfPops_genes@region.names, "chr..")
  # get the chromosome names
starts <- as.numeric(str_trim(str_extract(pfPops_genes@region.names, " \\d+ ")))
  # get the start positions of each gene
  # trim off the leading and trailing whitespaces

genes <- pfPops_genenames@region.names
  # get the gene names for
  # every gene in the genome

mut_genes <- unique(appendix$GeneID)
  # get unique genes from appendix
  # these are mutants in cp4 but not cp3
  # missense and nonsense


##########
### PREPARE DF FOR PLOTTING
##########

all <- as.data.frame(cbind(starts, fsts, chrs, genes))
  # make a df of relevant data
complete <- all[!all$fsts == "NaN",]
  # remove the NaN, missing fst values
complete$starts <- as.numeric(as.character(complete$starts))
  # fix encoding
complete$fsts <- as.numeric(as.character(complete$fsts))
  # fix encoding
complete$cp4_defining <- as.integer(complete$genes %in% mut_genes)
  # add a column


##########
### DO THE PLOTTING
##########

fst.plotter <- function(chr, name) {
  
  subset <- complete[complete$chrs == chr,]
  drivers <- subset[subset$cp4_defining == 1,]
  riders <- subset[subset$cp4_defining == 0,]
  chr_max <- max(complete[complete$chrs == chr,]$starts)
  
  plot(subset$fsts ~ subset$starts, 
       xlim = c(0,3290953), 
       ylim = c(0,1.1), 
       axes = FALSE, type = "n",
       xlab = "", ylab = "")
  polygon(c(0, 0, chr_max, chr_max), c(0,1,1,0), col = "gray95", border = NA)
  points(riders$fsts ~ riders$starts, pch = 19, cex = 0.5, col = "gray40")
  points(drivers$fsts ~ drivers$starts, pch = 19, cex = 0.5, col = "brown1")
  axis(1, at = c(1, chr_max), line = 0.5, labels = FALSE, tick = FALSE)
  axis(2, at = c(0,1), line = -1, las = 2, col = "grey40", col.axis = "grey40")
  mtext(name, side = 2, outer = FALSE, line = 1.5, cex = 0.8, col = "grey40")
  mtext(expression(italic(F)["ST"]), side = 2, outer = FALSE, line = 0.5, cex = 0.5, col = "grey40")
  
  if (nrow(genedata[genedata$chr == chr,]) > 0){
    
    print("happy")
    
    segments(x0 = genedata[genedata$chr == chr,]$pos, y0 = 0, 
             x1 = genedata[genedata$chr == chr,]$pos, y1 = 1, 
             col = "brown1", lwd = 1)
    
    mtext(genedata[genedata$chr == chr,]$name, side = 1, 
          at = genedata[genedata$chr == chr,]$pos, 
          col = "brown1", font = 4, cex = 0.5, line = -0.25)
    
  }
  

}

svg("fst_by_pc.svg", width = 7, height = 9)

par(mfrow = c(14,1), mar = c(0.75,3,0,0))

fst.plotter("chr01", "Chr 1")
fst.plotter("chr02", "Chr 2")
fst.plotter("chr03", "Chr 3")
fst.plotter("chr04", "Chr 4")
fst.plotter("chr05", "Chr 5")
fst.plotter("chr06", "Chr 6")
fst.plotter("chr07", "Chr 7")
fst.plotter("chr08", "Chr 8")
fst.plotter("chr09", "Chr 9")
fst.plotter("chr10", "Chr 10")
fst.plotter("chr11", "Chr 11")
fst.plotter("chr12", "Chr 12")
fst.plotter("chr13", "Chr 13")
fst.plotter("chr14", "Chr 14")

dev.off()

