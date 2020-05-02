#' @title Intergrative genes statistic
#' @description Calculate genes summary statistic across multiple datasets
#' 
#' @param allGenes Vector of all genes names for the analysis.
#' @param dataList A list of expression matrices, in which rows are genes and columns are samples.
#' @param groupList A list of vectors indicating sample group corresponding with expression matrices in dataList.
#' @param ncores Number of core to use in prallel processing.
#' @param method Function for combining p-values. It must accept one input which is a vector of p-values and return a combined p-value. Three methods are embeded in this package are addCLT, fisherMethod, and stoufferMethod.
#' 
#' @details 
#' 
#' To estimate the effect sizes of genes across all studies, first standardized mean difference 
#' for each gene in individual studies is compute. Next, the overall efect size and standard error are estimated using
#' the random-efects model. This overall efect size represents the gene's expression change under the efect of
#' the condition. The, z-scores and p-values of observing such efect sizes are computed. The p-values is obtained 
#' from classical hypothesis testing. By default, linear model and empirical Bayesian testing \(limma\) are used 
#' to compute the p-values for diferential expression. The two-tailed p-values are converted to one-tailed p-values (lef- and right-tailed). 
#' For each gene, the one-tailed p-values across all datasets are then combined using the addCLT, stouffer or fisher method.
#' These p-values represent how likely the diferential expression is observed by chance.
#' 
#' @return A data.frame of gene statistics with following columns:
#' 
#' \describe{
#'   \item{pTwoTails}{Two-tailed p-values}
#'   \item{pTwoTails.fdr}{Two-tailed p-values with false discovery rate correction}
#'   \item{pLeft}{left-tailed p-values}
#'   \item{pLeft.fdr}{left-tailed p-values with false discovery rate correction}
#'   \item{pRight.fdr}{right-tailed p-values with false discovery rate correction}
#'   \item{pRight}{right-tailed p-values}
#'   \item{ES}{Effect size}
#'   \item{ES.pTwoTails}{Two-tailed p-values for effect size}
#'   \item{ES.pTwoTails.fdr}{Two-tailed p-values for effect size with false discovery rate correction}
#'   \item{ES.pLeft}{Left-tailed p-values for effect size}
#'   \item{ES.pLeft.fdr}{Left-tailed p-values for effect size with false discovery rate correction}
#'   \item{ES.pRight}{Right-tailed p-values for effect size}
#'   \item{ES.pRight.fdr}{Right-tailed p-values for effect size with false discovery rate correction}
#' }
#' 
#' @author Tin Nguyen, Hung Nguyen, and Sorin Draghici
#' 
#' @references
#' 
#' Nguyen, T., Shafi, A., Nguyen, T. M., Schissler, A. G., & Draghici, S. (2020). NBIA: a network-based integrative analysis framework-applied to pathway analysis. Scientific reports, 10(1), 1-11.
#' Nguyen, T., Tagett, R., Donato, M., Mitrea, C., & Draghici, S. (2016). A novel bi-level meta-analysis approach: applied to biological pathway analysis. Bioinformatics, 32(3), 409-416.
#' Smyth, G. K. (2005). Limma: linear models for microarray data. In Bioinformatics and computational biology solutions using R and Bioconductor (pp. 397-420). Springer, New York, NY.
#' 
#' @seealso \code{\link{addCLT}}
#' 
#' @examples
#' 
#' datasets <- c("GSE17054", "GSE57194", "GSE33223", "GSE42140")
#' data(list = datasets, package = "BLMA")
#' dataList <- lapply(datasets, function(dataset) {
#'     get(paste0("data_", dataset))
#' })
#' groupList <- lapply(datasets, function(dataset) {
#'     get(paste0("group_", dataset))
#' })
#' names(dataList) <- datasets
#' names(groupList) <- datasets
#' 
#' allGenes <- Reduce(intersect, lapply(dataList, rownames))
#' 
#' geneStat <- getStatistics(allGenes, dataList, groupList)
#' head(geneStat)
#' 
#' 
#' @import limma
#' @import metafor
#' @export 
getStatistics <- function (allGenes, dataList, groupList, ncores = 1, method=addCLT) {
    # calculate left and right p-values
    pValue.Left <- pValue.Right <- data.frame(row.names=allGenes)
    for (i in 1:length(dataList)){
        dataset=names(dataList)[i]
        
        data <- dataList[[i]]
        group <- groupList[[i]]
        
        data=data[rownames(data)%in%allGenes,]
        
        gset=ExpressionSet(as.matrix(data))
        fl <- as.factor(group)
        gset$description <- fl
        design <- model.matrix(~ description + 0, gset)
        colnames(design) <- levels(fl)
        fit <- lmFit(gset, design)
        cont.matrix <- makeContrasts(d-c, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
        tT$P.Left = ifelse(tT$logFC <0, tT$P.Value/2, 1-tT$P.Value/2)
        tT$P.Right = ifelse(tT$logFC >0, tT$P.Value/2, 1-tT$P.Value/2)
        pValue.Left[rownames(tT),dataset] <- tT[,"P.Left"]
        pValue.Right[rownames(tT),dataset] <- tT[,"P.Right"]
    }
    pValue.Left$pAddCLT <- apply(pValue.Left[,1:length(dataList)], FUN = function(x) {method(na.omit(x))}, MARGIN = 1)
    pValue.Right$pAddCLT <- apply(pValue.Right[,1:length(dataList)], FUN = function(x) {method(na.omit(x))}, MARGIN = 1)
    
    
    controlSize <- controlMean <- controlStd <- data.frame(row.names=allGenes)
    diseaseSize <- diseaseMean <- diseaseStd <- data.frame(row.names=allGenes)
    for (i in 1:length(dataList)){
        dataset=names(dataList)[i]
        
        data <- dataList[[i]]
        group <- groupList[[i]]
        
        data=data[rownames(data)%in%allGenes,]
        
        controlDat = data[,names(group)[group=="c"]]
        diseaseDat = data[,names(group)[group=="d"]]
        
        controlSize[rownames(controlDat), dataset] <- ncol(controlDat)
        diseaseSize[rownames(diseaseDat), dataset] <- ncol(diseaseDat)
        
        controlMean[rownames(controlDat), dataset] <- apply(controlDat, FUN=mean, MARGIN = 1)
        diseaseMean[rownames(diseaseDat), dataset] <- apply(diseaseDat, FUN=mean, MARGIN = 1)
        
        controlStd[rownames(controlDat), dataset] <- apply(controlDat, FUN=sd, MARGIN = 1)
        diseaseStd[rownames(diseaseDat), dataset] <- apply(diseaseDat, FUN=sd, MARGIN = 1)
    }
    
    fc=mclapply(allGenes, FUN=calculateFC, 
                controlSize=controlSize, controlMean=controlMean, controlStd=controlStd,
                diseaseSize=diseaseSize, diseaseMean=diseaseMean, diseaseStd=diseaseStd,
                mc.cores=ncores)
    names(fc) <- allGenes
    effectSizes <- unlist(lapply(fc,FUN=function(x){x$beta}))
    pValues <- unlist(lapply(fc,FUN=function(x){x$pval}))
    
    mRNAStats = data.frame(row.names = allGenes, pLeft=pValue.Left[allGenes,]$pAddCLT, 
                           pRight=pValue.Right[allGenes,]$pAddCLT, 
                           ES=effectSizes[allGenes], ES.pTwoTails=pValues[allGenes])
    mRNAStats$pTwoTails=2*apply(mRNAStats[,c("pLeft","pRight")], FUN=min, MARGIN = 1)
    mRNAStats$ES.pLeft = ifelse(mRNAStats$ES <0, mRNAStats$ES.pTwoTails/2, 1-mRNAStats$ES.pTwoTails/2)
    mRNAStats$ES.pRight = ifelse(mRNAStats$ES >0, mRNAStats$ES.pTwoTails/2, 1-mRNAStats$ES.pTwoTails/2)
    mRNAStats$pLeft.fdr=p.adjust(mRNAStats$pLeft,method="fdr")
    mRNAStats$pRight.fdr=p.adjust(mRNAStats$pRight,method="fdr")
    mRNAStats$pTwoTails.fdr=p.adjust(mRNAStats$pTwoTails,method="fdr")
    mRNAStats$ES.pLeft.fdr=p.adjust(mRNAStats$ES.pLeft,method="fdr")
    mRNAStats$ES.pRight.fdr=p.adjust(mRNAStats$ES.pRight,method="fdr")
    mRNAStats$ES.pTwoTails.fdr=p.adjust(mRNAStats$ES.pTwoTails,method="fdr")
    
    mRNAStats
}


getSimilarityFromGrouping <- function(g) {
    N=length(g)
    S=matrix(0,N,N)
    colnames(S)=rownames(S)=names(g)
    for (j in unique(g)){
        X=rep(0,N);
        X[which(g==j)]=1
        S=S+X%*%t(X)
    }
    S
}

hierClustering <- function(data, method=km, minSize=40) {
    group <- rep(1, nrow(data));names(group) <- rownames(data)
    elements <- rownames(data)
    myGap <- clusGap(data[elements,],FUNcluster=method, K.max=2, B=50)
    split= length(elements)>minSize & maxSE(myGap$Tab[,"gap"], myGap$Tab[,"SE.sim"], method="globalSEmax")>1
    while (split==TRUE) {
        g <- method(data[elements,], 2)$cluster + max(group)
        group[elements] <- g
        
        split <- FALSE
        for (i in sort(unique(group))) {
            elements = names(group)[which(group==i)]
            if (length(elements)>minSize) {
                myGap <- clusGap(data[elements,],FUNcluster=method, K.max=2, B=50)
                if (maxSE(myGap$Tab[,"gap"], myGap$Tab[,"SE.sim"], method="globalSEmax")>1) {
                    split <- TRUE
                    break
                }
            }
        }
        
    }
    
    group <- as.numeric(as.factor(group))
    names(group) <- rownames(data)
    
    group
}

getCommonGenes <- function(dataList, minGeneCount=3) {
    allGenes=NULL
    for (i in 1:length(dataList)) {
        data <- dataList[[i]]
        allGenes=union(allGenes, rownames(data))
    }
    
    dsCounts <- rep(0, length(allGenes))
    names(dsCounts) <- allGenes
    for (i in 1:length(dataList)) {
        data <- dataList[[i]]
        dsCounts[rownames(data)] <- dsCounts[rownames(data)] +1
    } 
    
    allGenes <- names(dsCounts)[dsCounts>=minGeneCount]
    
    allGenes
}

pORACalc=function(DEGeneNames,measuredGenes,geneSet,minSize) {
    universe=measuredGenes
    white=DEGeneNames
    black=setdiff(measuredGenes,DEGeneNames)
    sample=intersect(measuredGenes,geneSet)
    sampleWhite=intersect(sample,DEGeneNames)
    if (length(sampleWhite)<minSize) {res=NA}
    else {
        res=phyper(q=length(sampleWhite)-1,m=length(white),n=length(black),k=length(sample), lower.tail=FALSE)
    }
    res
}

calculateFC = function(gInd, controlSize, controlMean, controlStd, diseaseSize, diseaseMean, diseaseStd) {
    # get dataset names without NA
    ds=controlSize[gInd,]
    ds=names(ds)[!is.na(ds)]
    
    tmp=data.frame(Study=seq(length(ds)),Source=ds, 
                   n1i=as.numeric(diseaseSize[gInd,ds]), m1i=as.numeric(diseaseMean[gInd,ds]), sd1i=as.numeric(diseaseStd[gInd,ds]), 
                   n2i=as.numeric(controlSize[gInd,ds]), m2i=as.numeric(controlMean[gInd,ds]), sd2i=as.numeric(controlStd[gInd,ds]))
    
    tmp1=escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i,m2i=m2i, sd2i=sd2i, n2i=n2i, data=tmp)
    
    res <- try(rma.uni(yi, vi, data=tmp1, control=list(stepadj=0.5,maxiter=10000), method="REML"),silent=TRUE)
    res
}
