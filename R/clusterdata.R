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
names(clustersummary)
clustersummary <- basic.stats(clusterdata)

clustersummary$Ho
Ho <- colMeans(clustersummary$Ho)
Hs <- colMeans(clustersummary$Hs)
Fis <- colMeans(clustersummary$Fis, na.rm = TRUE)
means <- cbind(Ho, Hs, Fis)
means

plot(clustersummary$Ho, clustersummary$Fis)

HWalldata <- hw.test(MyData, B = 1000)
plot(HW, col="Coral2", pch=18)

HWcluster <- hw.test(clusterdata, B = 1000)
plot(HWcluster, pch=18, col = "aquamarine") 
head(HW)
summary(HW)
HW

# PCA

PCA2 <- indpca(clusterdata) 
plot(PCA2, cex = 0.7, col="darkmagenta")

PCA <- indpca(MyData) 
plot(PCA, cex = 0.7, col="aquamarine")

# DAPC

kclust <- find.clusters(clusterdata, max.n.clust=50)

dapc2 <- dapc(clusterdata, kclust$kclust)
dapc2
scatter(dapc2)
