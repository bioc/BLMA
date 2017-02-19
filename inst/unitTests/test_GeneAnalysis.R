test_intraAnalysisClassic <- function() {
    x <- rnorm(10, mean = 0)
    a <- intraAnalysisClassic(x, func=t.test, mu=1, alternative = "less")
    l=list(x[1:5],x[6:10])
    b=addCLT(unlist(lapply(l, FUN=function(x) t.test(x, mu=1, alternative="less")$p.value)))
    checkEquals(a,b)
}

test_bilevelAnalysisClassic <- function() {
    l1 = lapply(as.list(seq(3)),FUN=function (x) rnorm(n=10, mean=1))
    a=addCLT(unlist(lapply(l1, FUN=intraAnalysisClassic, alternative="greater")))
    b=bilevelAnalysisClassic(x=l1, alternative="greater")
    checkEquals(a,b)
}

test_intraAnalysisGene <- function() {
    data(GSE17054)
    data = data_GSE17054
    group=group_GSE17054
    X = intraAnalysisGene(data, group)
    checkEquals(rownames(data),rownames(X))
}

test_bilevelAnalysisGene <- function() {
    data(GSE17054)
    data = data_GSE17054
    group=group_GSE17054
    dataList=list(data)
    groupList=list(group)
    names(dataList)=names(groupList)="GSE17054"
    a=intraAnalysisGene(data, group)
    a=a[order(a$pIntra),]
    b=bilevelAnalysisGene(dataList = dataList, groupList = groupList)
    checkEqualsNumeric(a$pIntra,b$pBilevel)
}