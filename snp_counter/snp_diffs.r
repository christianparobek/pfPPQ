## To plot the pairwise differennces between groups

library(ggplot2)


cp1 <- read.table("CP1.out")
cp1 <- cp1[upper.tri(cp1, diag = FALSE)]

cp2 <- read.table("CP2.out")
cp2 <- cp2[upper.tri(cp2, diag = FALSE)]

cp3 <- read.table("CP3.out")
cp3 <- cp3[upper.tri(cp3, diag = FALSE)]

cp4 <- read.table("CP4.out")
cp4 <- cp4[upper.tri(cp4, diag = FALSE)]

cp1 <- cbind(cp1, 1)
cp2 <- cbind(cp2, 2)
cp3 <- cbind(cp3, 3)
cp4 <- cbind(cp4, 4)

df <- as.data.frame(rbind(cp1, cp2, cp3, cp4))

names(df) <- c("snps", "cp")

plot(df[,1] ~ df[,2])

svg("diffs.svg", width = 4, height = 3)
ggplot(df, aes(cp, snps)) + 
  geom_jitter() + 
  theme_bw() + 
  xlab("CP Group") + 
  ylab("Pairwise SNP Differences")
dev.off()
