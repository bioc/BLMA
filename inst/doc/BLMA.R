### R code from vignette source 'BLMA.Rnw'

###################################################
### code chunk number 1: BLMA.Rnw:110-138
###################################################
# one-sample tests
library(BLMA)
set.seed(1)
x=rnorm(10, mean = 0)
# one-sample left-tailed t-test
t.test(x, mu=1, alternative = "less")$p.value
# one-sample left-tailed intra-experiment analysis with t-test
IntraAnalysisClassic(x, func=t.test, mu=1, alternative = "less")
# one-sample right-tailed t-test
t.test(x, mu=1, alternative = "greater")$p.value
# one-sample right-tailed intra-experiment analysis with t-test
IntraAnalysisClassic(x, func=t.test, mu=1, alternative = "greater")
# one-sample two-tailed t-test
t.test(x, mu=1)$p.value
# one-sample two-tailed intra-experiment analysis with t-test
IntraAnalysisClassic(x, func=t.test, mu=1)
# one-sample left-tailed Wilcoxon test
wilcox.test(x, mu=1, alternative = "less")$p.value
# one-sample left-tailed intra-experiment analysis with Wilcoxon test
IntraAnalysisClassic(x, func=wilcox.test, mu=1, alternative = "less")
# one-sample right-tailed Wilcoxon test
wilcox.test(x, mu=1, alternative = "greater")$p.value
# one-sample right-tailed intra-experiment analysis with Wilcoxon test
IntraAnalysisClassic(x, func=wilcox.test, mu=1, alternative = "greater")
# one-sample two-tailed Wilcoxon test
wilcox.test(x, mu=1)$p.value
# one-sample two-tailed intra-experiment analysis with Wilcoxon test
IntraAnalysisClassic(x, func=wilcox.test, mu=1)


###################################################
### code chunk number 2: BLMA.Rnw:143-158
###################################################
# two-sample tests
set.seed(1)
x=rnorm(20, mean=0); y=rnorm(20, mean=1)
# two-sample left-tailed t-test
t.test(x,y,alternative="less")$p.value
# two-sample left-tailed intra-experiment analysis with t-test
IntraAnalysisClassic(x, y, func=t.test, alternative = "less")
# two-sample right-tailed t-test
t.test(x,y,alternative="greater")$p.value
# two-sample right-tailed intra-experiment analysis with t-test
IntraAnalysisClassic(x, y, func=t.test, alternative = "greater")
# two-sample two-tailed t-test
t.test(x,y)$p.value
# two-sample two-tailed intra-experiment analysis with t-test
IntraAnalysisClassic(x, y, func=t.test)


###################################################
### code chunk number 3: BLMA.Rnw:164-185
###################################################
# one-sample tests
set.seed(1)
l1 = lapply(as.list(seq(3)),FUN=function (x) rnorm(n=10, mean=1))
l0 = lapply(as.list(seq(3)),FUN=function (x) rnorm(n=10, mean=0))
# one-sample right-tailed t-test
lapply(l1, FUN=function(x) t.test(x, alternative="greater")$p.value)
# combining the p-values of one-sample t-test:
addCLT(unlist(lapply(l1, FUN=function(x) 
    t.test(x, alternative="greater")$p.value)))
#Bi-level meta-analysis with one-sample right-tailed t-test
BilevelAnalysisClassic(x=l1, func=t.test, alternative="greater")
# two-sample left-tailed t-test
lapply(seq(l1), FUN=function(i,l1,l0) 
    t.test(l1[[i]], l0[[i]], alternative="greater")$p.value, l1, l0)
# combining the p-values of one-sample t-test:
addCLT(unlist(lapply(seq(l1), FUN=function(i,l1,l0) 
    t.test(l1[[i]], l0[[i]], alternative="greater")$p.value, l1, l0)))
#Bi-level meta-analysis with two-sample right-tailed t-test
BilevelAnalysisClassic(x=l1, y=l0, func=t.test, alternative="greater")
#Bi-level meta-analysis with two-sample left-tailed t-test
BilevelAnalysisClassic(x=l1, y=l0, func=t.test, alternative="less")


###################################################
### code chunk number 4: BLMA.Rnw:226-259
###################################################
library(BLMA)
# load KEGG pathways and create genesets
x=loadKEGGPathways()
gslist=lapply(x$kpg,FUN=function(y){return (nodes(y));})
gs.names=x$kpn[names(gslist)]

# load the 6 AML datasets
dataSets=c("GSE14924_CD4", "GSE17054", "GSE57194", "GSE33223", 
           "GSE42140", "GSE8023")
data(list=dataSets, package="BLMA")

# prepare dataList and groupList
dataList <- list()
groupList <- list()
for (i in 1:length(dataSets)) {
   dataset=dataSets[i]
   group <- get(paste("group_",dataset,sep=""))
   data=get(paste("data_",dataset,sep=""))
   dataList[[i]] = data
   groupList[[i]] = group
}
names(dataList)=names(groupList)=dataSets

# perform bi-level meta-analysis in conjunction with ORA
t1=Sys.time()
ORAComb=BilevelAnalysisGeneset(gslist = gslist, gs.names = gs.names, 
            dataList = dataList, groupList = groupList, enrichment = "ORA")
t2=Sys.time(); 
# running time
t2-t1
#print the results
options(digits=3)
ORAComb[1:10, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")]


###################################################
### code chunk number 5: BLMA.Rnw:274-285
###################################################
# perform bi-level meta-analysis in conjunction with GSA
t1=Sys.time()
GSAComb=BilevelAnalysisGeneset(gslist = gslist, gs.names = gs.names, 
            dataList = dataList, groupList = groupList, enrichment = "GSA", 
            mc.cores=2, nperms=200, random.seed = 1)
t2=Sys.time(); 
# running time
t2-t1
#print the results
options(digits=3)
GSAComb[1:10, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")]


###################################################
### code chunk number 6: BLMA.Rnw:300-311
###################################################
set.seed(1)
t1=Sys.time()
PADOGComb=BilevelAnalysisGeneset(gslist = gslist, gs.names = gs.names, 
            dataList = dataList, groupList = groupList, enrichment = "PADOG", 
            mc.cores=2, NI=200)
t2=Sys.time(); 
# running time
t2-t1
#print the results
options(digits=3)
PADOGComb[1:10, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")]


###################################################
### code chunk number 7: BLMA.Rnw:321-331
###################################################
x=loadKEGGPathways()  
t1=Sys.time()
IAComb=BilevelAnalysisPathway(kpg = x$kpg, kpn = x$kpn, dataList = dataList, 
                              groupList = groupList, mc.cores = 2)
t2=Sys.time(); 
# running time
t2-t1
#print the results
options(digits=3)
IAComb[1:10, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")]


###################################################
### code chunk number 8: BLMA.Rnw:359-381
###################################################
#perform intra-experiment analysis of the dataset GSE14924_CD4 using addCLT
library(BLMA)
data(GSE14924_CD4)
t1=Sys.time()
X = IntraAnalysisGene(data_GSE14924_CD4, group_GSE14924_CD4)
Sys.time()-t1
X = X[order(X$pIntra), ]
# top 10 genes
X[1:10,]
# bottom 10 genes
X[(nrow(X)-10):nrow(X),]

#perform intra-experiment analysis of GSE14924_CD4 using Fisher's method
t1=Sys.time()
Y = IntraAnalysisGene(data_GSE14924_CD4, group_GSE14924_CD4, 
                      metaMethod=fishersMethod)
Sys.time()-t1
Y = Y[order(Y$pIntra), ]
# top 10 genes
Y[1:10,]
# bottom 10 genes
Y[(nrow(Y)-10):nrow(Y),]


###################################################
### code chunk number 9: BLMA.Rnw:393-416
###################################################
dataSets=c("GSE14924_CD4", "GSE17054", "GSE57194", "GSE33223", "GSE42140", 
           "GSE8023")
data(list=dataSets, package="BLMA")
dataList <- list()
groupList <- list()
for (i in 1:length(dataSets)) {
   dataset=dataSets[i]
   group <- get(paste("group_",dataset,sep=""))
   data=get(paste("data_",dataset,sep=""))
   dataList[[i]] = data
   groupList[[i]] = group
}
names(dataList)=names(groupList)=dataSets

# running time
t1=Sys.time()
Z=BilevelAnalysisGene(dataList = dataList, groupList = groupList)
Sys.time()-t1

# top 10 genes
Z[1:10,]
# bottom 10 genes
Z[(nrow(Z)-10):nrow(Z),]


