library(usethis)
library(hierfstat)
library(ggplot2)
library(diveRsity)
library(adegenet)
library(plotly)
library(pegas)
library(ade4)
library(HardyWeinberg)

setwd("~/Desktop/oysterdata")
MyData <- read.genepop("oysterdata.GEN", ncode = 2)
clusterdata <- read.genepop("ByCluster.GEN", ncode = 2)

# My Data
MyData <- na.omit(MyData)
head(MyData)
tab(MyData)
div <- summary(MyData)

MySummary<-basic.stats(MyData)

# HW

HWalldata <- hw.test(MyData, B = 1000)
plot(HWalldata, col="Coral2", pch=18)

# means for table all data

alldatano <- colMeans(MySummary$n.ind.samp)
alldataHo <- colMeans(MySummary$Ho)
alldataHs <- colMeans(MySummary$Hs)
alldataFis <- colMeans(MySummary$Fis, na.rm = TRUE)
alldataHW <- colMeans(HWalldata)
alldatameans <- cbind(alldatano, alldataHo, alldataHs, alldataFis, alldataHW)
MEANS1 <- alldatameans

# means for cluster data

clustersummary <- basic.stats(clusterdata)

Ho <- colMeans(clustersummary$Ho)
Hs <- colMeans(clustersummary$Hs)
Fis <- colMeans(clustersummary$Fis, na.rm = TRUE)
means <- cbind(Ho, Hs, Fis)
MEANS2 <- means


# plots 

plot(MySummary$Ho, MySummary$Hs, pch=20, cex=0.7)
abline(0,1,lty=2)

pop(MyData)
barplot(table(pop(MyData)), col=funky(13), las=2,
        xlab="Population", ylab="Sample size")
temp <- locNames(MyData, withAlleles=TRUE)

plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus", col="coral2", pch=20)


# PCA ALL DATA

PCA1 <- indpca(MyData) 
plot(PCA1, cex = 0.7, col="Coral2")

# PCA CLUSTERS

PCA2 <- indpca(clusterdata) 
plot(PCA2, cex = 0.7, col="Coral2")

# PCA SCATTER ALL DATA

library(poppr)
poppr(
  MyData,
  total = TRUE,
  sublist = "ALL",
  exclude = NULL,
  blacklist = NULL,
  sample = 0,
  method = 1,
  missing = "ignore",
  cutoff = 0.05,
  quiet = FALSE,
  clonecorrect = FALSE,
  strata = 1,
  keep = 1,
  plot = TRUE,
  hist = TRUE,
  index = "rbarD",
  minsamp = 10,
  legend = FALSE)

names(div$pop.n.all)

dev.off()

sum(is.na(MyData$tab))
X <- tab(MyData, freq = TRUE, NA.method = "mean")
class(X)
dim(X)
X[1:5,1:5]
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

s.label(pca1$li)
title("PCA of Oyster dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)
pca1

# We can use s.class to represent both
# the genotypes and inertia ellipses for populations.

s.class(pca1$li, pop(MyData))
title("PCA of Oyster dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)


s.class(pca1$li,pop(MyData),xax=1,yax=3,sub="PCA 1-3",csub=2)
title("PCA of Oyster dataset\naxes 1-3")
add.scatter.eig(pca1$eig[1:20],nf=3,xax=1,yax=3)

#  we can remove the grid, choose different colors for the groups, use larger dots
# and transparency to better assess the density of points, and remove internal segments of the
# ellipses:

col <- funky(15)
s.class(pca1$li, pop(MyData),xax=1,yax=3, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE)

# Colors are based on the first three PCs of the PCA, recoded respectively on the red, green,
# and blue channel. In this figure, the genetic diversity is represented in two complementary
# ways: by the distances (further away = more genetically different), and by the colors (more
# different colors = more genetically different).

colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA of Oyster dataset\naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

# We can represent the diversity on the third axis similarly

colorplot(pca1$li[c(1,3)], pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
title("PCA of Oyster dataset\naxes 1-3")
abline(v=0,h=0,col="grey", lty=2)

# DAPC ALL DATA

grp <- find.clusters(MyData, max.n.clust=50)
names(grp)
head(grp$Kstat, 8)
grp$stat
head(grp$grp, 10)
table(pop(MyData), grp$grp)
table.value(table(pop(MyData), grp$grp), col.lab=paste("inf", 1:6),
            row.lab=paste("ori", 1:6))
dapc1 <- dapc(MyData, grp$grp)
scatter(dapc1)
# Colour and style
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)
myCol <- c("darkblue","purple","green","orange","red","blue")
scatter(dapc1, posi.da="bottomright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft")
# Add legend
scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:6))

# can also add a minimum spanning tree based on the (squared) distances between
# populations within the entire space. This allows one to bear in mind the actual proximities
# between populations inside the entire space, which are not always well represented in susbsets
# of discriminant functions of lesser rank. We also indicate the centre of each group with
# crosses. Lastly, we remove the DAPC eigenvalues, not very useful in this case, and replace
# them manually by a graph of PCA eigenva

scatter(dapc1, ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=myCol, solid=.4, cex=3, clab=0,
        mstree=TRUE, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, txt.leg=paste("Cluster",1:6))
par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="black")
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=myCol)
myInset <- function(){
  temp <- dapc1$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc1$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated variance (%)",
       cex=1, pch=20, type="h", lwd=2)
}
add.scatter(myInset(), posi="bottomright",
            inset=c(-0.03,-0.01), ratio=.28,
            bg=transp("white"))
scatter(dapc1,1,1, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)

