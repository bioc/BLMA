test_loadKEGGPathways <- function() {
    x=loadKEGGPathways()
    checkIdentical(names(x$kpg), names(x$kpn))
}

test_bilevelAnalysisGeneset <- function() {
    dataSets=c("GSE17054")
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
    
    # load KEGG pathways and create gene sets
    x=loadKEGGPathways()
    gslist=lapply(x$kpg,FUN=function(y){return (nodes(y));})
    gs.names=x$kpn[names(gslist)]
    
    ORAComb=bilevelAnalysisGeneset(gslist = gslist, gs.names = gs.names, 
            dataList = dataList, groupList = groupList, enrichment = "ORA")
    
    GSAComb=bilevelAnalysisGeneset(gslist = gslist, gs.names = gs.names, 
            dataList = dataList, groupList = groupList, enrichment = "GSA", 
            mc.cores=1, nperms=200, random.seed = 1)
    
    PADOGComb=bilevelAnalysisGeneset(gslist = gslist, gs.names = gs.names, 
            dataList = dataList, groupList = groupList, enrichment = "PADOG", 
            mc.cores=1, NI=200)
    
    checkEqualsNumeric(dim(ORAComb),dim(GSAComb))
    checkEqualsNumeric(dim(ORAComb),dim(PADOGComb))
}

test_bilevelAnalysisPathway <- function() {
    dataSets=c("GSE17054")
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
    
    x=loadKEGGPathways() 
    IAComb=bilevelAnalysisPathway(kpg = x$kpg, kpn = x$kpn, dataList = dataList, 
                                  groupList = groupList, mc.cores = 1)
    
    checkTrue(setequal(rownames(IAComb),names(x$kpg)))
}

