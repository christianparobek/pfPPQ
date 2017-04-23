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

data <- read.table("PPQ_res_2017-01-30_modified_for_heatmap_SRRs.txt", 
                   header=TRUE, sep = "\t", na.strings = ".")
# read in data
# make the periods NAs

data <- data[!data$SRS == 0,]
# remove non-WGS samples

data$ppq_ic90[data$ppq_ic90 > 3000] <- 3000 # knock down the really high values
data$ppq_ic90 <- as.numeric(scale(data$ppq_ic90))
data$mq_ic50 <- as.numeric(scale(data$mq_ic50))
# scale the data


##############################
######### SUBSET DATA ########
##############################

## by ic50s
ic50.subsetter <- function(data, group, order_by) {
  data <- subset(data, select=c("SRS", "ppq_ic90", "mq_ic50"))[data$pca_group == group,] # subset by group
  data[is.na(data)] <- -500
  data$SRS <- with(data, eval(parse(text = paste("reorder(SRS,", order_by, ")")))) # order
  # need the eval() expression to pass the "order_by" string as part of the reorder() expression
  #data <- na.omit(data) # remove rows with missing ic50s
  data <- melt(data) # for ggplot
  return(data)
}

cp1 <- ic50.subsetter(data, 1, "mq_ic50")
cp2 <- ic50.subsetter(data, 2, "mq_ic50")
cp3 <- ic50.subsetter(data, 3, "ppq_ic90")
cp4 <- ic50.subsetter(data, 4, "ppq_ic90")

## by generic genotype subsetter
subsetter <- function(data, group, order_by, selection) {
  data <- subset(data, select=selection)[data$pca_group == group,] # subset by group
  data[,c('ppq_ic90','mq_ic50')][is.na(data[,c('ppq_ic90','mq_ic50')])] <- -50 # replace mq and ppq ic50 NA values with an integer, -50
  data$SRS <- with(data, eval(parse(text = paste("reorder(SRS,", order_by, ")")))) # order
  data <- subset(data, select=-c(ppq_ic90, mq_ic50)) # remove the ppq and mq ic50 cols
  data <- data.frame(lapply(data, as.factor)) # convert all cols to factors
  data <- melt(data, id.vars = "SRS")
  return(data)
}

## subset ppq blocks
ppq_cols <- c("SRS", "exo_E415G", "Plas_2_3_dup", 
              "ppq_ic90", "mq_ic50")
cp1_ppq <- subsetter(data, 1, "mq_ic50", ppq_cols)
cp2_ppq <- subsetter(data, 2, "mq_ic50", ppq_cols)
cp3_ppq <- subsetter(data, 3, "ppq_ic90", ppq_cols)
cp4_ppq <- subsetter(data, 4, "ppq_ic90", ppq_cols)

## by mq genotypes
mq_cols <- c("SRS", "mdrcn_high", "ppq_ic90", "mq_ic50")
cp1_mq <- subsetter(data, 1, "mq_ic50", mq_cols)
cp2_mq <- subsetter(data, 2, "mq_ic50", mq_cols)
cp3_mq <- subsetter(data, 3, "ppq_ic90", mq_cols)
cp4_mq <- subsetter(data, 4, "ppq_ic90", mq_cols)

## by art genotypes (Miotto's backbone)
miotto_cols <- c("SRS", "miotto_mdr2", "miotto_arps10", 
                 "miotto_crt326", "miotto_crt356", "miotto_fd", 
                 "miotto_pph", "ppq_ic90", "mq_ic50")
cp1_art <- subsetter(data, 1, "mq_ic50", miotto_cols)
cp2_art <- subsetter(data, 2, "mq_ic50", miotto_cols)
cp3_art <- subsetter(data, 3, "ppq_ic90", miotto_cols)
cp4_art <- subsetter(data, 4, "ppq_ic90", miotto_cols)

## by K13 genotypes
k13_cols <- c("SRS", "k13_r539t_present", "k13_c580y_present", 
              "k13_other_present", "ppq_ic90", "mq_ic50")
cp1_k13 <- subsetter(data, 1, "mq_ic50", k13_cols)
cp2_k13 <- subsetter(data, 2, "mq_ic50", k13_cols)
cp3_k13 <- subsetter(data, 3, "ppq_ic90", k13_cols)
cp4_k13 <- subsetter(data, 4, "ppq_ic90", k13_cols)

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
aesthetics <- aes(as.factor(variable), SRS, group=SRS)

## Setting colors for gradient and binary 
min <- min(c(data$ppq_ic90, data$mq_ic50), na.rm = TRUE)
max <- max(c(data$ppq_ic90, data$mq_ic50), na.rm = TRUE)
gradient <- scale_fill_gradient(low = "grey90", high = "black", limits = c(min, max), na.value = "#ffffbf")
bw <- scale_fill_manual(values = c("0" = "#91bfdb", "1" = "#fc8d59"), na.value = "#ffffbf") # original way


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

## Plot Label plots
cp1_lab <- ggplot(cp1[cp1$variable == "ppq_ic90",], aesthetics) + geom_text(aes(label=SRS), size = 2) + theme
cp2_lab <- ggplot(cp2[cp2$variable == "ppq_ic90",], aesthetics) + geom_text(aes(label=SRS), size = 2) + theme
cp3_lab <- ggplot(cp3[cp3$variable == "ppq_ic90",], aesthetics) + geom_text(aes(label=SRS), size = 2) + theme
cp4_lab <- ggplot(cp4[cp4$variable == "ppq_ic90",], aesthetics) + geom_text(aes(label=SRS), size = 2) + theme


###############################
########## PLOT IT ############
###############################

svg("heatmap_ic90_new.svg", width=4, height=8)

grid.arrange(p_cp1_art, p_cp1_k13, p_cp1_mq, p_cp1_ppq, p_cp1, cp1_lab, p_cp2_art, p_cp2_k13, p_cp2_mq, p_cp2_ppq, p_cp2, cp2_lab,
             p_cp3_art, p_cp3_k13, p_cp3_mq, p_cp3_ppq, p_cp3, cp3_lab, p_cp4_art, p_cp4_k13, p_cp4_mq, p_cp4_ppq, p_cp4, cp4_lab,
             ncol = 6, widths = c(5.75/14,3/14,1/14,2/14,2/14,2.5/14), heights = c(17/78, 18/78, 20/78, 23/78),
             top = "", left = "", bottom = "", padding = unit(1, "line"))

## Place column-group labels
grid.text("Background", x = 0.25, y = 0.985, gp=gpar(fontsize=9))
grid.text("K13", x = 0.50, y = 0.985, gp=gpar(fontsize=9))
grid.text("MQ", x = 0.61, y = 0.985, gp=gpar(fontsize=9))
grid.text("PPQ", x = 0.695, y = 0.985, gp=gpar(fontsize=9))
grid.text("IC90/50", x = 0.81, y = 0.985, gp=gpar(fontsize=9))

## Place line segments over columns
grid.segments(x0 = 0.10, y0 = 0.97, x1 = 0.40, y1 = 0.97)
grid.segments(x0 = 0.425, y0 = 0.97, x1 = 0.57, y1 = 0.97)
grid.segments(x0 = 0.59, y0 = 0.97, x1 = 0.63, y1 = 0.97)
grid.segments(x0 = 0.65, y0 = 0.97, x1 = 0.74, y1 = 0.97)
grid.segments(x0 = 0.762, y0 = 0.97, x1 = 0.853, y1 = 0.97)

## Name the CP groups
grid.text("CP1", rot = 90, x = 0.04, y = 0.86, gp=gpar(fontsize=9))
grid.text("CP2", rot = 90, x = 0.04, y = 0.66, gp=gpar(fontsize=9))
grid.text("CP3", rot = 90, x = 0.04, y = 0.43, gp=gpar(fontsize=9))
grid.text("CP4", rot = 90, x = 0.04, y = 0.18, gp=gpar(fontsize=9))

## Add the gene names
grid.text(c("mdr2", "arps10", "crt326", "crt356", "fd", "pph"), rot = 45, x = ((0:5)/20)+0.13, y = 0.04, gp=gpar(fontsize=6), just = "right")
grid.text(c("R539T", "C580Y", "other"), rot = 45, x = ((0:2)/20)+0.45, y = 0.04, gp=gpar(fontsize=6), just = "right")
grid.text(c("mdr1"), rot = 45, x = 0.615, y = 0.04, gp=gpar(fontsize=6), just = "right")
grid.text(c("exo415", "plas2/3"), rot = 45, x = ((0:1)/20)+0.675, y = 0.04, gp=gpar(fontsize=6), just = "right")
grid.text(c("PPQ", "MQ"), rot = 45, x = ((0:1)/20)+0.79, y = 0.04, gp=gpar(fontsize=6), just = "right")

dev.off()


