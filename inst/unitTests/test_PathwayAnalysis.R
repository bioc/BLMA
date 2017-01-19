test_loadKEGGPathways <- function() {
    x=loadKEGGPathways()
    checkIdentical(names(x$kpg), names(x$kpn))
}

test_BilevelAnalysisGeneset <- function() {
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
    
    ORAComb=BilevelAnalysisGeneset(gslist = gslist, gs.names = gs.names, 
            dataList = dataList, groupList = groupList, enrichment = "ORA")
    
    GSAComb=BilevelAnalysisGeneset(gslist = gslist, gs.names = gs.names, 
            dataList = dataList, groupList = groupList, enrichment = "GSA", 
            mc.cores=1, nperms=200, random.seed = 1)
    
    PADOGComb=BilevelAnalysisGeneset(gslist = gslist, gs.names = gs.names, 
            dataList = dataList, groupList = groupList, enrichment = "PADOG", 
            mc.cores=1, NI=200)
    
    checkEqualsNumeric(dim(ORAComb),dim(GSAComb))
    checkEqualsNumeric(dim(ORAComb),dim(PADOGComb))
}

test_BilevelAnalysisPathway <- function() {
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
    IAComb=BilevelAnalysisPathway(kpg = x$kpg, kpn = x$kpn, dataList = dataList, 
                                  groupList = groupList, mc.cores = 1)
    
    checkTrue(setequal(rownames(IAComb),names(x$kpg)))
}

