# Library
install.packages("xtable")
library(usethis)
library(hierfstat)
library(ggplot2)
library(diveRsity)
library(adegenet)
library(plotly)
library(pegas)
library(ade4)
# install.packages("HardyWeinberg")
library(HardyWeinberg)
library(captioner)


setwd("~/Desktop/oysterdata")
MyData <- read.genepop("oysterdata.GEN", ncode = 2)
MyData <- na.omit(MyData)
head(MyData)
tab(MyData)
summary(MyData)
div <- summary(MyData)
names(div)
head(div)
div

MySummary<-basic.stats(MyData)
MySummary$n.ind.samp
MySummary$Ho
mean(MySummary$Ho)

# HW 

HWalldata <- hw.test(MyData, B = 1000)
plot(HW, col="Coral2", pch=18)
head(HW)
summary(HW)
HW

# means for table

alldatano <- colMeans(MySummary$n.ind.samp)
alldataHo <- colMeans(MySummary$Ho)
alldataHs <- colMeans(MySummary$Hs)
alldataFis <- colMeans(MySummary$Fis, na.rm = TRUE)
alldataHW <- colMeans(HWalldata)
alldatameans <- cbind(alldatano, alldataHo, alldataHs, alldataFis, alldataHW)
alldatameans


# seppop
?seppop
splits <- seppop(MyData,pop=MyData@pop,truenames=TRUE,res.type=c("genind","matrix"),
                 drop=FALSE, treatOther=TRUE, quiet=TRUE)


PSHB2 <- splits$PSHB2
PSHB2summary <- basic.stats(PSHB2)

pshb2ho <- (PSHB2summary$Ho)
pshb2ho <- as.data.frame(pshb2ho)
SummaryHo <- MySummary$Ho
pshb2ho$pop <- c(rep("PSHB2", nrow(pshb2ho)))
PSHB2 <- as.data.frame
SummaryHo <- as.data.frame(SummaryHo)
SummaryHo$pop <- c(rep("ALL", nrow(SummaryHo)))


my_summaries <- rbind(SummaryHo, pshb2ho)
library(dplyr)

SummaryHo <- SummaryHo[,-9] # dropped the 9th column
colMeans(SummaryHo)

SummaryHs <- MySummary$Hs

SummaryHs <- SummaryHs[,-9] # dropped the 9th column
colMeans(SummaryHs)
meansHs

SummaryFis <- MySummary$Fis
SummaryFis <- SummaryFis[,-9] # dropped the 9th column
colMeans(SummaryFis,na.rm = TRUE)

my_summaries <- as.data.frame(my_summaries)

HoMeans <- group_by(my_summaries, pop) %>%
  summarise(
    count = n(), mean = mean(Ho, na.rm = TRUE),
    sd = sd(Ho, na.rm = TRUE)
  )

# Pull out each stat seperately for each pop (Ho, He, etc). 
# Create c(rep("B2", nrow(B2_Ho))

par(mfrow=c(2,2))
barplot(div$Hexp-div$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")
barplot(div$n.by.pop, main="Sample sizes per population",
        ylab="Number of genotypes",las=3)
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus", col="coral2", pch=20)

plot(div$Hobs, div$Hexp, xlab="Observed Heterozygosity", ylab="Expected Heterozygosity", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")

plot((div$Hexp - div$Hobs)/(div$Hexp), xlab="Observed Heterozygosity", ylab="Expected Heterozygosity", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")
abline

testhw <- (div$Hexp - div$Hobs)/(div$Hexp)
summary(testhw)

figs(name="Figure 1","Expected heterozygosity as a function of observed heterozygosity per locus")


###
plot((div$Hexp - div$Hobs)/(div$Hexp))
test <- (div$Hexp - div$Hobs)
boot.ppfst(dat=MyData,nboot=100,quant=c(0.025,0.975),diploid=TRUE)

pop <- pop(MyData)
basicstat <- basic.stats(MyData, diploid = TRUE, digits = 2)
div
###

names(basicstat)
head(basicstat)
boot.ppfis(MyData)

PCA <- indpca(MyData) 
plot(PCA, cex = 0.7, col="Coral2")



MyClusters <- find.clusters(MyData, max.n.clust = 50)
head(MyClusters)
MyClusters

barplot(div$Hexp-div$Hobs, main="Heterozygosity: expected-observed", ylab="Hexp - Hobs")
   
barplot(div$n.by.pop, main="Sample sizes per population",
        ylab="Number of genotypes",las=3)

t.test(div$Hexp,div$Hobs,pair=T,var.equal=TRUE,alter="greater")


oyster.hwt <- hw.test(MyData, B=0)
oyster.hwt

Fst(as.loci(MyData))

MyData@pop

B2 <- seppop(MyData)$PSHB2
B3 <- seppop(MyData)$WBB2 
B4 <- seppop(MyData)$QBB2

temp2 <- inbreeding(B2, N=100)
class(temp2)
head(names(temp2))
head(temp2[[1]],20)
Fbar2 <- sapply(temp2, mean)
hist(Fbar2, col="firebrick", main="Average inbreeding in PSHB2")

temp3 <- inbreeding(B3, N=100)
class(temp3)
head(names(temp3))
head(temp3[[1]],20)
Fbar3 <- sapply(temp3, mean)
hist(Fbar3, col="firebrick", main="Average inbreeding in WBB2")

temp4 <- inbreeding(B4, N=100)
class(temp4)
head(names(temp4))
head(temp4[[1]],20)
Fbar4 <- sapply(temp4, mean)
hist(Fbar4,breaks=3, col=rgb(1,0,0,0.5), xlab="Average inbreeding coefficient", 
     ylab="Frequency", main="Average inbreeding in QBB2")

hist(Fbar3, breaks=3, col=rgb(0,0,1,0.5), add=T)
?par
par(
  mfrow=c(1,2),
  mar=c(4,4,1,0)
)

hist(Fbar4,breaks=3, col=rgb(1,0,0,0.5), xlab="Average inbreeding coefficient", 
     ylab="Frequency")
hist(Fbar3, breaks=3, col=rgb(0,0,1,0.5), xlab="Average inbreeding coefficient", ylab = "")
legend("topright", legend=c("Quibray Bay","Woolooware Bay"), col=c(rgb(1,0,0,0.5), 
                                                      rgb(0,0,1,0.5)), pt.cex=2, pch=10 )

inbreed <- rbind(Fbar2, Fbar3, Fbar4)



class(inbreed)
INBData <- as.data.frame(inbreed)

head(INBData)
INBData$`55`

dev.off()
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

# sPCA - couldn't work out! 


# Discriminant Analysis of Principal Components (DAPC)

MyData
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



names(div)
div$Hobs
lm <- lm(div$Hexp~div$Hobs)
plot(lm)
lm2 <- lm(div$n.by.pop~div$pop.n.all)
plot(lm2)

qqnorm((div$Hexp - div$Hobs)/div$Hexp)
qqline((div$Hexp - div$Hobs)/div$Hexp)
t.test((div$Hexp - div$Hobs)/div$Hexp)

qqnorm(div$Hexp - div$Hobs)

hist(div$Hexp - div$Hobs)
plot(div$Hexp ~ div$Hobs)
qqnorm(div$Hexp - div$Hobs)
qqline(div$Hexp - div$Hobs)



genpop <- genind2genpop(
  MyData,
  pop = NULL,
  quiet = FALSE,
  process.other = FALSE,
  other.action = mean
)

install.packages("gtsummary")
library(gtsummary)

trial2 <- data1 %>% select(Hobs, pop.n.all, Hexp)

install.packages("captioner")
library(captioner)
figs <- captioner(prefix="Figure")
tbls <- captioner(prefix="Table")



# kable in knitr for tables 

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

dapc3 <- dapc(MyData, grp$grp)
scatter(dapc3)

