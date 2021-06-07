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



