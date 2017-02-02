######
# Take JP's SAS data and try out heatmaps in ggplot
# Started 23 July 2015
######

################################################################
###### ggplot2 #################################################
################################################################

##############################
####### LOAD LIBRARIES #######
##############################
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)

##############################
##### READ AND CLEAN DATA ####
##############################

data <- read.table("../pca_vs_ic50/data/PPQ_res_2017-01-30.txt", 
                   header=TRUE, sep = "\t", na.strings = ".")
# read in data
# make the periods NAs

data <- data[!is.na(data$WGS_ID),]
# remove non-WGS samples

data <- data[order(data$pca_group),]
# sort the data frame by pca group

data$ppq_ic50[data$ppq_ic50 > 200] <- 200 
# knock down all high ppq values to 150


##############################
######### SUBSET DATA ########
##############################

## by ic50s
ic50.subsetter <- function(data, group, order_by) {
  data <- subset(data, select=c("WGS_ID", "ppq_ic90", "mq_ic50"))[data$pca_group == group,] # subset by group
  data[is.na(data)] <- -500
  data$WGS_ID <- with(data, eval(parse(text = paste("reorder(WGS_ID,", order_by, ")")))) # order
    # need the eval() expression to pass the "order_by" string as part of the reorder() expression
  #data <- na.omit(data) # remove rows with missing ic50s
  data <- melt(data) # for ggplot
  return(data)
}

cp1 <- ic50.subsetter(data, 1, "mq_ic50")
cp2 <- ic50.subsetter(data, 2, "mq_ic50")
cp3 <- ic50.subsetter(data, 3, "ppq_ic50")
cp4 <- ic50.subsetter(data, 4, "ppq_ic50")

## by ppq genotypes
ppq.subsetter <- function(data, group, order_by) {
  data <- subset(data, select=c("WGS_ID", "exo_E415G", "Plas_2_3_dup", "mIIg_present", "X5r_cn_high", "CRT_350", "cviet_present", "MDR6_repeats_gt6", "ppq_ic50", "mq_ic50"))[data$pca_group == group,] # subset by group
  data[,c('ppq_ic50','mq_ic50')][is.na(data[,c('ppq_ic50','mq_ic50')])] <- -50 # replace mq and ppq ic50 NA values with an integer, -50
  data$WGS_ID <- with(data, eval(parse(text = paste("reorder(WGS_ID,", order_by, ")")))) # order
  #data <- data[!with(data,is.na(ppq_ic50)),] # exclude rows w/NA in ppq col
  #data <- data[!with(data,is.na(mq_ic50)),] # exclude rows w/NA in mq col
  data <- subset(data, select=-c(ppq_ic50, mq_ic50)) # remove the ppq and mq ic50 cols
  data <- data.frame(lapply(data, as.factor)) # convert all cols to factors
  data <- melt(data, id.vars = "WGS_ID")
  return(data)
}

cp1_ppq <- ppq.subsetter(data, 1, "mq_ic50")
cp2_ppq <- ppq.subsetter(data, 2, "mq_ic50")
cp3_ppq <- ppq.subsetter(data, 3, "ppq_ic50")
cp4_ppq <- ppq.subsetter(data, 4, "ppq_ic50")

## by mq genotypes ###### Turn this into names
mq.subsetter <- function(data, group, order_by) {
  data <- subset(data, select=c("WGS_ID", "mdrcn_high", "ppq_ic50", "mq_ic50"))[data$pca_group == group,] # subset by group
  data[,c('ppq_ic50','mq_ic50')][is.na(data[,c('ppq_ic50','mq_ic50')])] <- -50 # replace mq and ppq ic50 NA values with an integer, -50
  data$WGS_ID <- with(data, eval(parse(text = paste("reorder(WGS_ID,", order_by, ")")))) # order
  #data <- data[!with(data,is.na(ppq_ic50)),] # exclude rows w/NA in ppq col
  #data <- data[!with(data,is.na(mq_ic50)),] # exclude rows w/NA in mq col
  data <- subset(data, select=-c(ppq_ic50, mq_ic50)) # remove the ppq and mq ic50 cols
  data <- data.frame(lapply(data, as.factor)) # convert all cols to factors
  data <- melt(data, id.vars = "WGS_ID")
  return(data)
}

cp1_mq <- mq.subsetter(data, 1, "mq_ic50")
cp2_mq <- mq.subsetter(data, 2, "mq_ic50")
cp3_mq <- mq.subsetter(data, 3, "ppq_ic50")
cp4_mq <- mq.subsetter(data, 4, "ppq_ic50")

## by art genotypes (Miotto's backbone)
art.subsetter <- function(data, group, order_by) {
  data <- subset(data, select=c("WGS_ID", "miotto_mdr2", "miotto_arps10", "miotto_crt326", "miotto_crt356", "miotto_fd", "miotto_pph", "ppq_ic50", "mq_ic50"))[data$pca_group == group,]
  data[,c('ppq_ic50','mq_ic50')][is.na(data[,c('ppq_ic50','mq_ic50')])] <- -50 # replace mq and ppq ic50 NA values with an integer, -50
  data$WGS_ID <- with(data, eval(parse(text = paste("reorder(WGS_ID,", order_by, ")")))) # order
  #data <- data[!with(data,is.na(ppq_ic50)),] # exclude rows w/NA in ppq col
  #data <- data[!with(data,is.na(mq_ic50)),] # exclude rows w/NA in mq col
  data <- subset(data, select=-c(ppq_ic50, mq_ic50)) # remove the ppq and mq ic50 cols
  data <- data.frame(lapply(data, as.factor)) # convert all cols to factors
  data <- melt(data, id.vars = "WGS_ID")
  return(data)
}

cp1_art <- art.subsetter(data, 1, "mq_ic50")
cp2_art <- art.subsetter(data, 2, "mq_ic50")
cp3_art <- art.subsetter(data, 3, "ppq_ic50")
cp4_art <- art.subsetter(data, 4, "ppq_ic50")

## by K13 genotypes
k13.subsetter <- function(data, group, order_by) {
  data <- subset(data, select=c("WGS_ID", "k13_r539t_present", "k13_c580y_present", "k13_other_present", "ppq_ic50", "mq_ic50"))[data$pca_group == group,]
  data[,c('ppq_ic50','mq_ic50')][is.na(data[,c('ppq_ic50','mq_ic50')])] <- -50 # replace mq and ppq ic50 NA values with an integer, -50
  data$WGS_ID <- with(data, eval(parse(text = paste("reorder(WGS_ID,", order_by, ")")))) # order
  #data <- data[!with(data,is.na(ppq_ic50)),] # exclude rows w/NA in ppq col
  #data <- data[!with(data,is.na(mq_ic50)),] # exclude rows w/NA in mq col
  data <- subset(data, select=-c(ppq_ic50, mq_ic50)) # remove the ppq and mq ic50 cols
  data <- data.frame(lapply(data, as.factor)) # convert all cols to factors
  data <- melt(data, id.vars = "WGS_ID")
  return(data)
}

cp1_k13 <- k13.subsetter(data, 1, "mq_ic50")
cp2_k13 <- k13.subsetter(data, 2, "mq_ic50")
cp3_k13 <- k13.subsetter(data, 3, "ppq_ic50")
cp4_k13 <- k13.subsetter(data, 4, "ppq_ic50")

##############################
##### DEFINE PLOT PARAMS #####
##############################

## DEFINE THEMES ##
theme <- theme(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  legend.position = "none",
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  plot.margin=unit(c(0,0,0,0),"cm"))

bottom_row_theme <- theme(axis.text = element_text(size = rel(0.5)), axis.text.y = element_blank()) 

left_row_theme <- theme(axis.title = element_text(size = rel(0.8)), plot.margin=unit(c(0,0,0,0),"cm"), axis.title.x = element_blank())

top_row_theme <- theme(plot.margin=unit(c(0.75,0,-0.25,-0.5),"cm"), plot.title = element_text(size = rel(0.8)))

## DEFINE GEOMETRY ##
geometry <- geom_tile(aes(fill = value), colour = "white", alpha = 0.8)

## DEFINE PLOT ##
aesthetics <- aes(variable, WGS_ID)

## Colorbrewer recommended colors
gradient <- scale_fill_gradient(low = "grey90", high = "black", limits = c(-50,200), na.value = "#ffffbf")
bw <- scale_fill_manual(values = c("0" = "#91bfdb", "1" = "#fc8d59"), na.value = "#ffffbf")


#################################
########## MAKE GROBS ###########
#################################

## PLOT IC50 COLUMN ##
p_cp1 <- ggplot(cp1, aesthetics) + geometry + gradient + theme
p_cp2 <- ggplot(cp2, aesthetics) + geometry + gradient + theme
p_cp3 <- ggplot(cp3, aesthetics) + geometry + gradient + theme
p_cp4 <- ggplot(cp4, aesthetics) + geometry + gradient + theme

## PLOT PPQ COLUMN ##
p_cp1_ppq <- ggplot(cp1_ppq, aesthetics) + geometry + theme + bw
p_cp2_ppq <- ggplot(cp2_ppq, aesthetics) + geometry + theme + bw
p_cp3_ppq <- ggplot(cp3_ppq, aesthetics) + geometry + theme + bw
p_cp4_ppq <- ggplot(cp4_ppq, aesthetics) + geometry + theme + bw

## PLOT MQ COLUMN ##
p_cp1_mq <- ggplot(cp1_mq, aesthetics) + geometry + theme + bw
p_cp2_mq <- ggplot(cp2_mq, aesthetics) + geometry + theme + bw
p_cp3_mq <- ggplot(cp3_mq, aesthetics) + geometry + theme + bw
p_cp4_mq <- ggplot(cp4_mq, aesthetics) + geometry + theme + bw

## PLOT K13 COLUMN ##
p_cp1_k13 <- ggplot(cp1_k13, aesthetics) + geometry + theme + bw
p_cp2_k13 <- ggplot(cp2_k13, aesthetics) + geometry + theme + bw
p_cp3_k13 <- ggplot(cp3_k13, aesthetics) + geometry + theme + bw
p_cp4_k13 <- ggplot(cp4_k13, aesthetics) + geometry + theme + bw

## PLOT ART COLUMN ##
p_cp1_art <- ggplot(cp1_art, aesthetics) + geometry + theme + bw
p_cp2_art <- ggplot(cp2_art, aesthetics) + geometry + theme + bw
p_cp3_art <- ggplot(cp3_art, aesthetics) + geometry + theme + bw
p_cp4_art <- ggplot(cp4_art, aesthetics) + geometry + theme + bw


###############################
########## PLOT IT ############
###############################

svg("heatmap_ic90.svg", width=5, height=8)

grid.arrange(p_cp1_art, p_cp1_k13, p_cp1_mq, p_cp1_ppq, p_cp1, p_cp2_art, p_cp2_k13, p_cp2_mq, p_cp2_ppq, p_cp2, 
             p_cp3_art, p_cp3_k13, p_cp3_mq, p_cp3_ppq, p_cp3, p_cp4_art, p_cp4_k13, p_cp4_mq, p_cp4_ppq, p_cp4, 
             ncol = 5, widths = c(6/19,3/19,1/19,7/19,2/19), heights = c(17/78, 18/78, 20/78, 23/78),
             top = "", left = "", bottom = "", right = "", padding = unit(1, "line"))

## Place column-group labels
grid.text("Background", x = 0.21, y = 0.985, gp=gpar(fontsize=10))
grid.text("K13", x = 0.415, y = 0.985, gp=gpar(fontsize=10))
grid.text("MQ", x = 0.505, y = 0.985, gp=gpar(fontsize=10))
grid.text("PPQ", x = 0.69, y = 0.985, gp=gpar(fontsize=10))
grid.text("IC90/50", x = 0.895, y = 0.985, gp=gpar(fontsize=10))

## Add the sample IDs
grid.text(rev(data[data$pca_group=="1",3]), x=0.97, y = ((0:16)/87)+0.770, gp=gpar(fontsize=6))
grid.text(rev(data[data$pca_group=="2",3]), x=0.97, y = ((0:17)/87)+0.558, gp=gpar(fontsize=6))
grid.text(rev(data[data$pca_group=="3",3]), x=0.97, y = ((0:19)/87)+0.323, gp=gpar(fontsize=6))
grid.text(rev(data[data$pca_group=="4",3]), x=0.97, y = ((0:22)/87)+0.053, gp=gpar(fontsize=6))

## Place line segments over columns
grid.segments(x0 = 0.08, y0 = 0.97, x1 = 0.335, y1 = 0.97)
grid.segments(x0 = 0.3525, y0 = 0.97, x1 = 0.47, y1 = 0.97)
grid.segments(x0 = 0.4875, y0 = 0.97, x1 = 0.52, y1 = 0.97)
grid.segments(x0 = 0.535, y0 = 0.97, x1 = 0.8375, y1 = 0.97)
grid.segments(x0 = 0.855, y0 = 0.97, x1 = 0.93, y1 = 0.97)

## Name the CP groups
grid.text("CP1", rot = 90, x = 0.04, y = 0.86)
grid.text("CP2", rot = 90, x = 0.04, y = 0.66)
grid.text("CP3", rot = 90, x = 0.04, y = 0.43)
grid.text("CP4", rot = 90, x = 0.04, y = 0.18)

## Add the gene names
grid.text(c("mdr2", "arps10", "crt326", "crt356", "fd", "pph"), rot = 45, x = ((0:5)/23)+0.10, y = 0.04, gp=gpar(fontsize=6), just = "right")
grid.text(c("R539T", "C580Y", "other"), rot = 45, x = ((0:2)/23)+0.375, y = 0.04, gp=gpar(fontsize=6), just = "right")
grid.text(c("exo415", "plas2/3", "mrp2", "X5r", "crt350", "CVIET", "mdr6"), rot = 45, x = ((0:6)/23)+0.56, y = 0.04, gp=gpar(fontsize=6), just = "right")
grid.text(c("mdr1"), rot = 45, x = 0.51, y = 0.04, gp=gpar(fontsize=6), just = "right")
grid.text(c("PPQ", "MQ"), rot = 45, x = ((0:1)/23)+0.875, y = 0.04, gp=gpar(fontsize=6), just = "right")

dev.off()

