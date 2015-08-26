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
library(reshape)
library(ggplot2)
library(grid)
library(gridExtra)

##############################
##### READ AND CLEAN DATA ####
##############################

data <- read.table("../pca_vs_ic50/PPQ_res_Aug20_v2.txt", 
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
  data <- data[data$pca_group == group, c(11,22,25)] # subset by group
  data$WGS_ID <- with(data, eval(parse(text = paste("reorder(WGS_ID,", order_by, ")")))) # order
    # need the eval() expression to pass the "order_by" string as part of the reorder() expression
  data <- na.omit(data) # remove rows with missing ic50s
  data <- melt(data) # for ggplot
  return(data)
}

cp1 <- ic50.subsetter(data, 1, "mq_ic50")
cp2 <- ic50.subsetter(data, 2, "mq_ic50")
cp3 <- ic50.subsetter(data, 3, "ppq_ic50")
cp4 <- ic50.subsetter(data, 4, "ppq_ic50")

## by ppq genotypes
ppq.subsetter <- function(data, group, order_by) {
  data <- data[data$pca_group == group, c(11,5,6,12,14,22,25)] # subset by group
  data$WGS_ID <- with(data, eval(parse(text = paste("reorder(WGS_ID,", order_by, ")")))) # order
  data <- data[!with(data,is.na(ppq_ic50)),] # exclude rows w/NA in ppq col
  data <- data[!with(data,is.na(mq_ic50)),] # exclude rows w/NA in mq col
  data <- subset(data, select=-c(ppq_ic50, mq_ic50)) # remove the ppq and mq ic50 cols
  data <- data.frame(lapply(data, as.factor)) # convert all cols to factors
  data <- melt(data, id.vars = "WGS_ID")
  return(data)
}

cp1_ppq <- ppq.subsetter(data, 1, "mq_ic50")
cp2_ppq <- ppq.subsetter(data, 2, "mq_ic50")
cp3_ppq <- ppq.subsetter(data, 3, "ppq_ic50")
cp4_ppq <- ppq.subsetter(data, 4, "ppq_ic50")

## by mq genotypes
mq.subsetter <- function(data, group, order_by) {
  data <- data[data$pca_group == group, c(11,28,22,25)] # subset by group
  data$WGS_ID <- with(data, eval(parse(text = paste("reorder(WGS_ID,", order_by, ")")))) # order
  data <- data[!with(data,is.na(ppq_ic50)),] # exclude rows w/NA in ppq col
  data <- data[!with(data,is.na(mq_ic50)),] # exclude rows w/NA in mq col
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
  data <- data[data$pca_group == group, c(11,16:21,22,25)] # subset by group
  data$WGS_ID <- with(data, eval(parse(text = paste("reorder(WGS_ID,", order_by, ")")))) # order
  data <- data[!with(data,is.na(ppq_ic50)),] # exclude rows w/NA in ppq col
  data <- data[!with(data,is.na(mq_ic50)),] # exclude rows w/NA in mq col
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
  data <- data[data$pca_group == group, c(11,32,31,22,25)] # subset by group
  data$WGS_ID <- with(data, eval(parse(text = paste("reorder(WGS_ID,", order_by, ")")))) # order
  data <- data[!with(data,is.na(ppq_ic50)),] # exclude rows w/NA in ppq col
  data <- data[!with(data,is.na(mq_ic50)),] # exclude rows w/NA in mq col
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
  plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))

bottom_row_theme <- theme(axis.text = element_text(angle = 45, hjust = 1, size = rel(0.8)), axis.text.y = element_blank())

left_row_theme <- theme(axis.title = element_text(size = rel(1)), plot.margin=unit(c(0,0,-0.25,0),"cm"), axis.title.x = element_blank())

top_row_theme <- theme(plot.margin=unit(c(0.5,0,-0.25,-0.5),"cm"), plot.title = element_text(size = rel(1)))

## DEFINE GENOMETRY ##
geometry <- geom_tile(aes(fill = value), colour = "white", alpha = 0.8)

## DEFINE PLOT ##
aesthetics <- aes(variable, WGS_ID)

## DEFINE SCALES ##
gradient <- scale_fill_gradient(low = "white", high = " black", limits = c(0,200))

## DEFINE COLOR SCHEMES ##
yellow <- scale_fill_manual(values = c("0" = "darkgoldenrod1", "1" = "darkgoldenrod"), na.value = "gray90")
green <- scale_fill_manual(values = c("0" = "darkseagreen1","1" = "darkseagreen3"), na.value = "gray90")
blue <- scale_fill_manual(values = c("0" = "steelblue1","1" = "steelblue3"), na.value = "gray90")
purple <- scale_fill_manual(values = c("0" = "mediumorchid1","1" = "mediumorchid4"), na.value = "gray90")

bw <- scale_fill_manual(values = c("0" = "gray", "1" = "black"), na.value = "white")

cp1_color <- bw
cp2_color <- bw
cp3_color <- bw
cp4_color <- bw

##############################
########## PLOT IT ###########
##############################

## PLOT IC50 COLUMN ##
p_cp1 <- ggplot(cp1, aesthetics) + geometry + gradient + theme + top_row_theme + labs(title = "IC50s")
p_cp2 <- ggplot(cp2, aesthetics) + geometry + gradient + theme
p_cp3 <- ggplot(cp3, aesthetics) + geometry + gradient + theme
p_cp4 <- ggplot(cp4, aesthetics) + geometry + gradient + theme + bottom_row_theme +
  scale_x_discrete(labels = c("ppq_ic50" = "PPQ", "mq_ic50" = "MQ"))

## PLOT PPQ COLUMN ##
p_cp1_ppq <- ggplot(cp1_ppq, aesthetics) + geometry + theme + cp1_color + top_row_theme + labs(title = "PPQ-R")
p_cp2_ppq <- ggplot(cp2_ppq, aesthetics) + geometry + theme + cp2_color
p_cp3_ppq <- ggplot(cp3_ppq, aesthetics) + geometry + theme + cp3_color
p_cp4_ppq <- ggplot(cp4_ppq, aesthetics) + geometry + theme + cp4_color + bottom_row_theme +
  scale_x_discrete(labels = c("CRT_350" = "crt350", "MDR6_repeats_gt6" = "mdr6", "MDR1_1246" = "mdr1", "MIIg_present" = "MIIg"))

## PLOT MQ COLUMN ##
p_cp1_mq <- ggplot(cp1_mq, aesthetics) + geometry + theme + cp1_color + top_row_theme + labs(title = "MQ-R")
p_cp2_mq <- ggplot(cp2_mq, aesthetics) + geometry + theme + cp2_color
p_cp3_mq <- ggplot(cp3_mq, aesthetics) + geometry + theme + cp3_color
p_cp4_mq <- ggplot(cp4_mq, aesthetics) + geometry + theme + cp4_color + bottom_row_theme +
  scale_x_discrete(labels = c("mdrcn_high" = "mdr CN"))

## PLOT K13 COLUMN ##
p_cp1_k13 <- ggplot(cp1_k13, aesthetics) + geometry + theme + cp1_color + top_row_theme + labs(title = "K13")
p_cp2_k13 <- ggplot(cp2_k13, aesthetics) + geometry + theme + cp2_color
p_cp3_k13 <- ggplot(cp3_k13, aesthetics) + geometry + theme + cp3_color
p_cp4_k13 <- ggplot(cp4_k13, aesthetics) + geometry + theme + cp4_color + bottom_row_theme +
  scale_x_discrete(labels = c( "k13_c580y_present" = "C580Y", "k13_r539t_present" = "R539T"))


## PLOT ART COLUMN ##
p_cp1_art <- ggplot(cp1_art, aesthetics) + geometry + theme + cp1_color + left_row_theme + top_row_theme + theme(plot.margin=unit(c(0.5,0,-0.25,0),"cm")) + labs(title = "Background") + ylab("CP1")
p_cp2_art <- ggplot(cp2_art, aesthetics) + geometry + theme + cp2_color + left_row_theme + ylab("CP2")
p_cp3_art <- ggplot(cp3_art, aesthetics) + geometry + theme + cp3_color + left_row_theme + ylab("CP3")
p_cp4_art <- ggplot(cp4_art, aesthetics) + geometry + theme + cp4_color + left_row_theme + bottom_row_theme + ylab("CP4") +
  scale_x_discrete(labels = c("miotto_mdr2" = "mdr2", "miotto_arps10" = "arps10", "miotto_crt326" = "crt326", "miotto_crt356" = "crt356", "miotto_fd" = "fd", "miotto_pph" = "pph"))

grid.arrange(p_cp1_art, p_cp1_k13, p_cp1_mq, p_cp1_ppq, p_cp1, 
             p_cp2_art, p_cp2_k13, p_cp2_mq, p_cp2_ppq, p_cp2, 
             p_cp3_art, p_cp3_k13, p_cp3_mq, p_cp3_ppq, p_cp3, 
             p_cp4_art, p_cp4_k13, p_cp4_mq, p_cp4_ppq, p_cp4, 
             ncol = 5, heights=c(0.9,.588,1,1.25), widths=c(3,1,0.5,2,1))
grid.text("ART-R", x = unit(0.35, "npc"), y = unit(0.95, "npc"))
