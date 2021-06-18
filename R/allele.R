install.packages("PopGenReport")
library(PopGenReport)

allelrich <- allel.rich(MyData)


colMeans(allelrich$pop.sizes)


allelrich_cluster <- allel.rich(clusterdata)


meanrichmydata <- allelrich$mean.richness
meaanrichcluster <- allelrich_cluster$mean.richness

meanrichmydata
meaanrichcluster

# F STATS
library(devtools)
install_github("jgx65/hierfstat")
library(hierfstat)

pop.freq(MyData)

pairwise.neifst(MyData, pop = NULL, res.type = c("dist", "matrix"))
?hierfstat

dataframe <- as.data.frame(MyData)

Fst(as.loci(MyData))
Gtest <- gstat.randtest(MyData,nsim=99)

require(devtools)
install.packages("devtools")
install_version("hierfstat", version = "0.04-22", repos = "http://cran.us.r-project.org")


fstat(MyData)
fstat(MyData, fstonly=TRUE)
mydata.gtest <- gstat.randtest(MyData)
mydata.gtest
plot(mydata.gtest)


mydata.matFst <- pairwise.fst(MyData,res.type="matrix")
mydata.matFst[1:4,1:4]
mydata.tree <- nj(mydata.matFst)
plot(mydata.tree, type="unr", tip.col=funky(nPop(MyData)), font=2)
annot <- round(mydata.tree$edge.length,2)
edgelabels(annot[annot>0], which(annot>0), frame="n")
add.scale.bar()

mydata.matFst

table.paint(mydata.matFst, col.labels=1:16, csize = 0.7, clegend = 0.9, 
            clabel.col=0.7, clabel.row = 0.7)
?table.paint
data_box <- mydata.matFst
diag(data_box) <- NA
boxplot(data_box, col=funky(nPop(MyData)),las=2,cex.lab = 0.9, cex.axis=0.7,
        xlab="Population", ylab="Fst")



pca1$eig[1]
pc1 <- pca1$li[,1]
var(pc1)
var(pc1)*89/90
mean(pc1^2)
n <- length(pc1)
0.5*mean(dist(pc1)^2)*((n-1)/n)
eig.perc <- 100*pca1$eig/sum(pca1$eig)
head(eig.perc)
s.arrow(pca1$c1)
loadingplot(pca1$c1^2)
contrib <- loadingplot(pca1$c1, axis=2,
                       thres=.065, lab.jitter=1)
contrib


install.packages("PopGenReport")
library(PopGenReport)
popgenreport(MyData, mk.counts = TRUE, mk.map = FALSE,
             maptype = "satellite", mapdotcolor = "blue", mapdotsize = 1,
             mapdotalpha = 0.4, mapdottype = 19, mapzoom = NULL,
             mk.locihz = FALSE, mk.hwe = FALSE, mk.fst = FALSE,
             mk.gd.smouse = FALSE, mk.gd.kosman = FALSE, mk.pcoa = FALSE,
             mk.spautocor = FALSE, mk.allele.dist = TRUE, mk.null.all = FALSE,
             mk.allel.rich = TRUE, mk.differ.stats = FALSE, mk.custom = FALSE,
             fname = "PopGenReport", foldername = "results", path.pgr = NULL,
             mk.Rcode = FALSE, mk.complete = FALSE, mk.pdf = TRUE)

R.Version()
