## Want to calculate pi for my CP3 and CP4
## As compared to the underlying miotto populations
## Updated 19 April 2016


###################################################
################# Load Libraries ##################
###################################################

library(PopGenome)
library(stringr)


#####################################################
################### Import Data #####################
#####################################################

cp3_mine <- readData("cp3_mine/vcf/", format="VCF", gffpath = "gff/")
cp3_miotto <- readData("cp3_miotto/vcf/", format="VCF", gffpath = "gff/")
cp4_mine <- readData("cp4_mine/vcf/", format="VCF", gffpath = "gff/")
cp4_miotto <- readData("cp4_miotto/vcf/", format="VCF", gffpath = "gff/")
  ## group membership was determined using pca_explorer.r
  ## in the k-means section at the bottom

#####################################################
################### Calculate Pi ####################
#####################################################

cp3_mine <- diversity.stats(cp3_mine, pi = TRUE)
cp3_miotto <- diversity.stats(cp3_miotto, pi = TRUE)
cp4_mine <- diversity.stats(cp4_mine, pi = TRUE)
cp4_miotto <- diversity.stats(cp4_miotto, pi = TRUE)

cp3_mine@Pi
cp3_miotto@Pi
cp4_mine@Pi
cp4_miotto@Pi

## Calculate the per-nucleotide diversity
## they all need the same denominator
## Since popgenome n.sites is supposed to represent length of chr
## But can differ from dataset to dataset
cp3_mine_pi <- cp3_mine@Pi/cp4_miotto@n.sites 
cp3_mioto_pi <- cp3_miotto@Pi/cp4_miotto@n.sites
cp4_mine_pi <- cp4_mine@Pi/cp4_miotto@n.sites
cp4_miotto_pi <- cp4_miotto@Pi/cp4_miotto@n.sites

table <- as.data.frame(cbind(cp3_mine_pi, cp3_mioto_pi, cp4_mine_pi, cp4_miotto_pi))
names(table) <- c("cp3_mine", "cp3_miotto", "cp4_mine", "cp4_miotto")

write.table(table, file = "pi.txt", row.names = FALSE, quote = FALSE, sep = "\t")
