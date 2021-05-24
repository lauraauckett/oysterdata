# Library
library(usethis)
library(hierfstat)
library(ggplot2)
library(diveRsity)
library(adegenet)
library(plotly)
library(pegas)
library(ade4)
install.packages("HardyWeinberg")
library(HardyWeinberg)

MyData <- read.genepop("oysterdata.GEN", ncode = 2)
head(MyData)
tab(MyData)
summary(MyData)
div <- summary(MyData)
div
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus", col="coral2", pch=20)
plot(div$Hobs, div$Hexp, xlab="Observed Heterozygosity", ylab="Expected Heterozygosity", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")

bartlett.test(list(div$Hexp, div$Hobs))
basic.stats(MyData[,-1])
boot.ppfst(dat=MyData,nboot=100,quant=c(0.025,0.975),diploid=TRUE)

pop(MyData)
basicstat <- basic.stats(MyData, diploid = TRUE, digits = 2)

names(basicstat)
head(basicstat)
boot.ppfis(MyData)
PCA <- indpca(MyData) 
plot(PCA, cex = 0.7, col="Coral2")
HW <- hw.test(MyData, B = 1000)
plot(HW, col="Coral2", pch=18)
head(HW)
MyClusters <- find.clusters(MyData, max.n.clust = 50)
head(MyClusters)
MyClusters

barplot(div$Hexp-div$Hobs, main="Heterozygosity: expected-observed", ylab="Hexp - Hobs")
   
barplot(div$n.by.pop, main="Sample sizes per population",
        ylab="Number of genotypes",las=3)

bartlett.test(list(div$Hexp,div$Hobs))
t.test(div$Hexp,div$Hobs,pair=T,var.equal=TRUE,alter="greater")

t.test.lm = (value ~ names, data=div) 
div2 <- as.data.frame(div) 

oyster.hwt <- hw.test(MyData, B=0)
oyster.hwt

Fst(as.loci(MyData))


fstat(MyData, pop=NULL, fstonly=FALSE)
pairwise.fst(x, pop=NULL, res.type=c("dist","matrix"), truenames=TRUE)

B2 <- seppop(MyData)$PSHB 
B2
B3 <- seppop(MyData)$WBB2 

temp2 <- inbreeding(B2, N=100)
class(temp2)
head(names(temp2))
head(temp2[[1]],20)
Fbar2 <- sapply(temp2, mean)
hist(Fbar2, col="firebrick", main="Average inbreeding in Hatchery B2")

temp3 <- inbreeding(B3, N=100)
class(temp3)
head(names(temp3))
head(temp3[[1]],20)
Fbar3 <- sapply(temp3, mean)
hist(Fbar3, col="firebrick", main="Average inbreeding in Woolooware Bay B2")


B4 <- seppop(MyData)$QBB2
B4
temp4 <- inbreeding(B4, N=100)
class(temp4)
head(names(temp4))
head(temp4[[1]],20)
Fbar4 <- sapply(temp4, mean)
hist(Fbar4,breaks=3, col=rgb(1,0,0,0.5), xlab="Average inbreeding coefficient", 
     ylab="Frequency", main="Average inbreeding in Woolooware Bay B2")
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

basic.stats(MyData,diploid=TRUE,digits=4)
dev.off()

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

poppr.amova(
  MyData,
  hier = NULL,
  clonecorrect = FALSE,
  within = TRUE,
  dist = NULL,
  squared = TRUE,
  freq = TRUE,
  correction = "quasieuclid",
  sep = "_",
  filter = FALSE,
  threshold = 0,
  algorithm = "farthest_neighbor",
  threads = 1L,
  missing = "loci",
  cutoff = 0.05,
  quiet = FALSE,
  method = c("ade4", "pegas"),
  nperm = 0
)

?strata
library(graph4lg)
genpop <- genind_to_genepop(MyData, output = "data.frame")
genpop
diffCalc(infile = 'genpop', outfile = 'MyOut', fst = FALSE, pairwise = FALSE,
         bs_locus = FALSE, bs_pairwise = FALSE, boots = NULL,
         ci_type = "individuals", alpha = 0.05, para = FALSE)


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

