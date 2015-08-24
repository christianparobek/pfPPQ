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
  data <- data[data$pca_group == group, c(11,12,6,14,22,25)] # subset by group
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
  data <- data[data$pca_group == group, c(11,28,5,22,25)] # subset by group
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

## by art genotypes
art.subsetter <- function(data, group, order_by) {
  data <- data[data$pca_group == group, c(11,16:21,31,32,22,25)] # subset by group
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

##############################
########### PLOT IT ##########
##############################

# define theme
new_theme <- theme(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  legend.position = "none",
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank()
)

# plot the ic50s
p_cp1 <- ggplot(cp1, aes(variable, WGS_ID)) + geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0,200)) + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp2 <- ggplot(cp2, aes(variable, WGS_ID)) + geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0,200)) + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp3 <- ggplot(cp3, aes(variable, WGS_ID)) + geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0,200)) + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp4 <- ggplot(cp4, aes(variable, WGS_ID)) + geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0,200)) + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))


# plot the ppq stuff
p_cp1_ppq <- ggplot(cp1_ppq, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp2_ppq <- ggplot(cp2_ppq, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp3_ppq <- ggplot(cp3_ppq, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp4_ppq <- ggplot(cp4_ppq, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))

# ppq up close
p_ppq <- ggplot(cp1_ppq, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p_ppq

# plot the mq stuff
p_cp1_mq <- ggplot(cp1_mq, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp2_mq <- ggplot(cp2_mq, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp3_mq <- ggplot(cp3_mq, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp4_mq <- ggplot(cp4_mq, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))

# mq up close
p_mq <- ggplot(cp1_mq, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p_mq

# plot the art stuff
p_cp1_art <- ggplot(cp1_art, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp2_art <- ggplot(cp2_art, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp3_art <- ggplot(cp3_art, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))
p_cp4_art <- ggplot(cp4_art, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") + 
  new_theme + theme(plot.margin=unit(c(0,0,-0.25,-0.5),"cm"))

# mq up close
p_art <- ggplot(cp1_art, aes(variable, WGS_ID)) + geom_tile(aes(fill = value), colour = "white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p_art

grid.arrange(p_cp1_art, p_cp1_mq, p_cp1_ppq, p_cp1, 
             p_cp2_art, p_cp2_mq, p_cp2_ppq, p_cp2, 
             p_cp3_art, p_cp3_mq, p_cp3_ppq, p_cp3, 
             p_cp4_art, p_cp4_mq, p_cp4_ppq, p_cp4, 
             ncol = 4, heights=c(.588,.588,1,1), widths=c(4,1,1.5,1))




gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}


n = 12
cols = gg_color_hue(4)

dev.new(width=4, height=4)
plot(1:n, pch=16, cex=2, col=cols)
