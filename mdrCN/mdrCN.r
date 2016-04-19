## Plot pfmdr1 copynumber plotter
## 18 April 2016


## Read in Library
library(ggplot2)


## Read in data
data <- read.table("../pca_vs_ic50/data/PPQ_res_Mar20.txt", header = TRUE, sep = "\t", na.strings = "999")

## Exclude CP0 data
new <- data[!data$pca_group == 0,]

## Look at some stats
## Difference between groups?
kruskal.test(pfmdr1_copynum ~ pca_group, data = new)

## Plot the stuff
svg("pfmdrCN.svg", width = 6, height = 5)

p <- ggplot(new, aes(factor(pca_group), pfmdr1_copynum))
p + geom_boxplot() + theme_bw() + + #scale_y_continuous(limits = c(1, 4)) + 
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14)) + 
  labs(x = "CP Groups", y = expression(paste(italic("pfmdr1"), "CN")))

dev.off()