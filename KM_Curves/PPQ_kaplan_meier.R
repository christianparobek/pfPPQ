## Kaplan-Meier curves for PPQ/MQ resistance project
## Originally built by Jon Parr
## Adapted for many submissions
## Most recent update 27 Jan 2017 for GBE resubmission
## Incorporates Jon Parr's edits sent on that date
## Uses file exported by SAS containing recrudescence data from Jess Manning


#####################################
######### IMPORT LIBRARIES ##########
#####################################

library(ggplot2)
library(survival)
library(grid)
library(gridExtra)


#####################################
###### DEFINE USEFUL VARIABLES ######
#####################################

palette <- c("deepskyblue", "brown1", "darkolivegreen3", "darkgoldenrod1")
  # color consistent with the rest of the paper
#fig1_colors <- c("#5a5b5f", "#5ca8dc", "#cd3a34", "#2cb673", "#f5ee79")
#fig1_colors_wo <- c("#5ca8dc", "#cd3a34", "#2cb673", "#f5ee79")
#custom_lines <- c("dotdash", "dashed", "dotted", "solid")
custom_lines2 <- c("solid", "solid", "solid", "solid")
cp_names <- c("CP1", "CP2", "CP3", "CP4")
ppqvar_names <- c("Plas and/or Exo Absent", "Plas and Exo Present")
ppqvar_colors <- c("#91bfdb", "#fc8d59")
  # for WT and mutant samples


#####################################
########### READ IN DATA ############
#####################################

df <- read.table("PPQ_recrud.txt", header=TRUE, sep = "\t")
df <- read.table("PPQ_recrud_jan27_2017_withPCT.txt", header=TRUE, sep = "\t")
df <- read.table("PPQ_recrud_jun08_2017_withPCT.txt", header=TRUE, sep = "\t")

#Unsuccessfully tried to create reciprocal of censor, to allow for cumulative incidence graph
#df$invcensor[censor == 0] <- 1
#df$invcensor[censor == 1] <- 0

###########################################################################
#First, imports a function called ggsurv found here:
# http://www.r-statistics.com/2013/07/creating-good-looking-survival-curves-the-ggsurv-function/
# - which I edited for consistency with other figures (edits are annotated with #)
###########################################################################

ggsurv <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                   cens.col = 'red', lty.est = 1, lty.ci = 2,
                   cens.shape = 3, back.white = F, xlab = 'Time',
                   ylab = 'Survival', main = ''){
  
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = ''){
    
    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)
    
    col <- ifelse(surv.col == 'gg.def', 'black', surv.col)
    
    pl <- ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est, size = 1.2)
    
    pl <- if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci, size = 1.2) +
        geom_step(aes(y = low), color = col, lty = lty.ci, size = 1.2)
    } else (pl)
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col, size = 3)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  
  ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = '') {
    n <- s$strata
    
    groups <- factor(unlist(strsplit(names
                                     (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), '='))[1]
    gr.df <- vector('list', strata)
    ind <- vector('list', strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]
    
    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1))
    }
    
    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)
    
    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group), size = 1.2)
    
    col <- if(length(surv.col == 1)){
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
    
    pl <- if(surv.col[1] != 'gg.def'){
      pl + col
    } else {pl + scale_colour_discrete(name = gr.name)}
    
    line <- if(length(lty.est) == 1){
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {scale_linetype_manual(name = gr.name, values = lty.est)}
    
    pl <- pl + line
    
    pl <- if(CI == T) {
      if(length(surv.col) > 1 && length(lty.est) > 1){
        stop('Either surv.col or lty.est should be of length 1 in order
             to plot 95% CI with multiple strata')
      }else if((length(surv.col) > 1 | surv.col == 'gg.def')[1]){
        pl + geom_step(aes(y = up, color = group), lty = lty.ci, size = 1.2) +
          geom_step(aes(y = low, color = group), lty = lty.ci, size = 1.2)
      } else{pl +  geom_step(aes(y = up, lty = group), col = surv.col, size = 1.2) +
          geom_step(aes(y = low,lty = group), col = surv.col, size = 1.2)}
    } else {pl}
    
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col, size = 3)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  pl <- if(strata == 1) {ggsurv.s(s, CI , plot.cens, surv.col ,
                                  cens.col, lty.est, lty.ci,
                                  cens.shape, back.white, xlab,
                                  ylab, main)
  } else {ggsurv.m(s, CI, plot.cens, surv.col ,
                   cens.col, lty.est, lty.ci,
                   cens.shape, back.white, xlab,
                   ylab, main)}
  pl
}



#####################################################################
###### SURVIVAL CURVE FOR ORIGINAL CP GROUPS (ORIGINAL FIGURE) ######
#####################################################################

recrud.serv <- survfit(Surv(df$Persontime,df$Censor) ~ df$pca_group, data = df)
km1 <- ggsurv(recrud.serv)

km2 <- km1 +
  theme_bw(18) +
  ggtitle(NULL) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.justification=c(0,0), legend.position=c(0.05,0.05),
    legend.background = element_rect(color="black", fill = "grey90")
    ) +
  labs(
    x = "Days", 
    y = "") +  
  scale_color_manual(labels=cp_names, values=palette) +
  scale_linetype_manual(labels=cp_names, values=custom_lines2) + 
  ylim(0, 1)


#####################################################################
##### SURVIVAL CURVE FOR BINARY PPQ VARIANTS (FOR GBE REVISION) #####
#####################################################################

recrud.serv_ppqvar <- survfit(Surv(df$Persontime,df$Censor) ~ df$both_exo_plas_binary, data = df)
km1_ppqvar <- ggsurv(recrud.serv_ppqvar)

km2_ppqvar <- km1_ppqvar +
  theme_bw(18) +
  ggtitle(NULL) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 15),
    legend.title = element_blank(),
    legend.justification=c(0,0), legend.position=c(0.05,0.05),
    legend.background = element_rect(color="black", fill = "grey90")
  ) +
  labs(
    x = "Days", 
    y = "") +  
  scale_color_manual(labels=ppqvar_names, values=ppqvar_colors) +
  scale_linetype_manual(labels=ppqvar_names, values=custom_lines2) + 
  ylim(0, 1)


#####################################
######## MAKE COMBINED PLOT #########
#####################################

svg("kaplans.svg", width = 6, height = 6)
  grid.arrange(km2, km2_ppqvar, ncol = 1, heights = c(1, 1.2))
  grid.text("Proportion of Recrudescence", hjust = 0.35, 
            x = unit(0.025, "npc"), gp=gpar(fontsize=15), rot = 90)
dev.off()


