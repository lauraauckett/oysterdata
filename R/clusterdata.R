library(hierfstat)
library(adegenet)
library(plotly)
library(pegas)
library(ade4)

setwd("~/Desktop/oysterdata")
clusterdata <- read.genepop("ByCluster.GEN", ncode = 2)

head(clusterdata)
summary(clusterdata)
names(clusterdata)
clusterdata$pop

clustersummary <- basic.stats(clusterdata)

clustersummary$Ho
Ho <- colMeans(clustersummary$Ho)
Hs <- colMeans(clustersummary$Hs)
Fis <- colMeans(clustersummary$Fis, na.rm = TRUE)
means <- cbind(Ho, Hs, Fis)
means

plot(clustersummary$Ho, clustersummary$Fis)


HW <- hw.test(clusterdata, B = 1000)
plot(HW, col="Coral2", pch=18)
head(HW)
summary(HW)
HW

PCA2 <- indpca(clusterdata) 
plot(PCA2, cex = 0.7, col="Coral2")
