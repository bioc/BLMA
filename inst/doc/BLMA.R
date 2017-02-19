### R code from vignette source 'BLMA.Rnw'

###################################################
### code chunk number 1: BLMA.Rnw:114-142
###################################################
# one-sample tests
library(BLMA)
set.seed(1)
x <- rnorm(10, mean = 0)
# one-sample left-tailed t-test
t.test(x, mu=1, alternative = "less")$p.value
# one-sample left-tailed intra-experiment analysis with t-test
intraAnalysisClassic(x, func=t.test, mu=1, alternative = "less")
# one-sample right-tailed t-test
t.test(x, mu=1, alternative = "greater")$p.value
# one-sample right-tailed intra-experiment analysis with t-test
intraAnalysisClassic(x, func=t.test, mu=1, alternative = "greater")
# one-sample two-tailed t-test
t.test(x, mu=1)$p.value
# one-sample two-tailed intra-experiment analysis with t-test
intraAnalysisClassic(x, func=t.test, mu=1)
# one-sample left-tailed Wilcoxon test
wilcox.test(x, mu=1, alternative = "less")$p.value
# one-sample left-tailed intra-experiment analysis with Wilcoxon test
intraAnalysisClassic(x, func=wilcox.test, mu=1, alternative = "less")
# one-sample right-tailed Wilcoxon test
wilcox.test(x, mu=1, alternative = "greater")$p.value
# one-sample right-tailed intra-experiment analysis with Wilcoxon test
intraAnalysisClassic(x, func=wilcox.test, mu=1, alternative = "greater")
# one-sample two-tailed Wilcoxon test
wilcox.test(x, mu=1)$p.value
# one-sample two-tailed intra-experiment analysis with Wilcoxon test
intraAnalysisClassic(x, func=wilcox.test, mu=1)


###################################################
### code chunk number 2: BLMA.Rnw:147-162
###################################################
# two-sample tests
set.seed(1)
x <- rnorm(20, mean=0); y=rnorm(20, mean=1)
# two-sample left-tailed t-test
t.test(x,y,alternative="less")$p.value
# two-sample left-tailed intra-experiment analysis with t-test
intraAnalysisClassic(x, y, func=t.test, alternative = "less")
# two-sample right-tailed t-test
t.test(x,y,alternative="greater")$p.value
# two-sample right-tailed intra-experiment analysis with t-test
intraAnalysisClassic(x, y, func=t.test, alternative = "greater")
# two-sample two-tailed t-test
t.test(x,y)$p.value
# two-sample two-tailed intra-experiment analysis with t-test
intraAnalysisClassic(x, y, func=t.test)


###################################################
### code chunk number 3: BLMA.Rnw:168-189
###################################################
# one-sample tests
set.seed(1)
l1 <- lapply(as.list(seq(3)),FUN=function (x) rnorm(n=10, mean=1))
l0 <- lapply(as.list(seq(3)),FUN=function (x) rnorm(n=10, mean=0))
# one-sample right-tailed t-test
lapply(l1, FUN=function(x) t.test(x, alternative="greater")$p.value)
# combining the p-values of one-sample t-test:
addCLT(unlist(lapply(l1, FUN=function(x) 
    t.test(x, alternative="greater")$p.value)))
#Bi-level meta-analysis with one-sample right-tailed t-test
bilevelAnalysisClassic(x=l1, func=t.test, alternative="greater")
# two-sample left-tailed t-test
lapply(seq(l1), FUN=function(i,l1,l0) 
    t.test(l1[[i]], l0[[i]], alternative="greater")$p.value, l1, l0)
# combining the p-values of one-sample t-test:
addCLT(unlist(lapply(seq(l1), FUN=function(i,l1,l0) 
    t.test(l1[[i]], l0[[i]], alternative="greater")$p.value, l1, l0)))
#Bi-level meta-analysis with two-sample right-tailed t-test
bilevelAnalysisClassic(x=l1, y=l0, func=t.test, alternative="greater")
#Bi-level meta-analysis with two-sample left-tailed t-test
bilevelAnalysisClassic(x=l1, y=l0, func=t.test, alternative="less")


###################################################
### code chunk number 4: BLMA.Rnw:232-253
###################################################
library(BLMA)
# load KEGG pathways and create genesets
x=loadKEGGPathways()
gslist <- lapply(x$kpg,FUN=function(y){return (nodes(y));})
gs.names <- x$kpn[names(gslist)]

# load the 4 AML datasets
dataSets <- c("GSE17054", "GSE57194", "GSE33223", "GSE42140")
data(list=dataSets, package="BLMA")

# prepare dataList and groupList
names(dataSets) <- dataSets
dataList <- lapply(dataSets, function(dataset) get(paste0("data_", dataset)))
groupList <- lapply(dataSets, function(dataset) get(paste0("group_", dataset)))

# perform bi-level meta-analysis in conjunction with ORA
system.time(ORAComb <- bilevelAnalysisGeneset(gslist, gs.names, dataList, 
                                              groupList, enrichment = "ORA"))
#print the results
options(digits=2)
head(ORAComb[, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")])


###################################################
### code chunk number 5: BLMA.Rnw:268-275
###################################################
# perform bi-level meta-analysis in conjunction with GSA
system.time(GSAComb <- bilevelAnalysisGeneset(gslist, gs.names, dataList, 
                                            groupList, enrichment = "GSA", 
                                            nperms=200, random.seed = 1))
#print the results
options(digits=2)
head(GSAComb[, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")])


###################################################
### code chunk number 6: BLMA.Rnw:288-294
###################################################
set.seed(1)
system.time(PADOGComb <- bilevelAnalysisGeneset(gslist, gs.names, dataList, 
                                    groupList, enrichment = "PADOG", NI=200))
#print the results
options(digits=2)
head(PADOGComb[, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")])


###################################################
### code chunk number 7: BLMA.Rnw:304-309
###################################################
x <- loadKEGGPathways()  
system.time(IAComb <- bilevelAnalysisPathway(x$kpg, x$kpn, dataList, groupList))
#print the results
options(digits=2)
head(IAComb[, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")])


###################################################
### code chunk number 8: BLMA.Rnw:337-355
###################################################
#perform intra-experiment analysis of the dataset GSE33223 using addCLT
library(BLMA)
data(GSE33223)
system.time(X <- intraAnalysisGene(data_GSE33223, group_GSE33223))
X <- X[order(X$pIntra), ]
# top 10 genes
head(X)
# bottom 10 genes
tail(X)

#perform intra-experiment analysis of GSE33223 using Fisher's method
system.time(Y <- intraAnalysisGene(data_GSE33223, group_GSE33223, 
                      metaMethod=fisherMethod))
Y = Y[order(Y$pIntra), ]
# top 10 genes
head(Y)
# bottom 10 genes
tail(Y)


###################################################
### code chunk number 9: BLMA.Rnw:367-373
###################################################
system.time(Z <- bilevelAnalysisGene(dataList = dataList, groupList = groupList))

# top 10 genes
head(Z)
# bottom 10 genes
tail(Z)


