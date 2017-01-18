test_IntraAnalysisClassic <- function() {
    x=rnorm(10, mean = 0)
    a=IntraAnalysisClassic(x, func=t.test, mu=1, alternative = "less")
    l=splitS(x)
    b=addCLT(unlist(lapply(l, FUN=function(x) t.test(x, mu=1, alternative="less")$p.value)))
    checkEquals(a,b)
}

test_BilevelAnalysisClassic <- function() {
    l1 = lapply(as.list(seq(3)),FUN=function (x) rnorm(n=10, mean=1))
    a=addCLT(unlist(lapply(l1, FUN=IntraAnalysisClassic, alternative="greater")))
    b=BilevelAnalysisClassic(x=l1, alternative="greater")
    checkEquals(a,b)
}

test_IntraAnalysisGene <- function() {
    data(GSE17054)
    data = data_GSE17054
    group=group_GSE17054
    X = IntraAnalysisGene(data, group)
    checkEquals(rownames(data),rownames(X))
}

test_BilevelAnalysisGene <- function() {
    data(GSE17054)
    data = data_GSE17054
    group=group_GSE17054
    dataList=list(data)
    groupList=list(group)
    names(dataList)=names(groupList)="GSE17054"
    a=IntraAnalysisGene(data, group)
    a=a[order(a$pIntra),]
    b=BilevelAnalysisGene(dataList = dataList, groupList = groupList)
    checkEqualsNumeric(a$pIntra,b$pBilevel)
}