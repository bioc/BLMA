#' @title Intra-experiment analysis in conjunction with classical 
#' hypothesis tests
#' @description Perform an intra-experiment analysis in conjunction with any 
#' of the classical hypothesis testing methods, such as t-test, 
#' Wilcoxon test, etc.
#' @param x a numeric vector of data values
#' @param y an optional numeric vector of values
#' @param splitSize the minimum number of size in each split sample. 
#' splitSize should be at least 3. By default, splitSize=5
#' @param metaMethod the method used to combine p-values. This should be one 
#' of addCLT (additive method [1]), fishersMethod (Fisher's method [5]), 
#' stoufferMethod (Stouffer's method [6]), max (maxP method [7]), 
#' or min (minP method [8])
#' @param func the name of the hypothesis test. By default func=t.test
#' @param p.value the component that returns the p-value after performing 
#' the test provided by the \emph{func} parameter. For example, the function 
#' t-test returns the class "htest" where the component "p.value" is the 
#' p-value of the test. By default, p.value="p.value"
#' @param ... additional parameters for \emph{func}
#' @details This function performs an intra-experiment analysis for the given 
#' sample(s) [1]. Given x as the numeric vector, this function first splits x 
#' into smaller samples with size \emph{splitSize}, performs hypothesis 
#' testing using \emph{func}, and then combines the p-values 
#' using \emph{metaMethod}
#' @return
#' intra-experiment p-value
#' @author
#' Tin Nguyen and Sorin Draghici
#' @references
#' [1] T. Nguyen, R. Tagett, M. Donato, C. Mitrea, and S. Draghici. A novel 
#' bi-level meta-analysis approach -- applied to biological pathway analysis. 
#' Bioinformatics, 32(3):409-416, 2016.
#' 
#' @seealso \code{\link{BilevelAnalysisClassic}}, \code{\link{IntraAnalysisGene}}, \code{\link{BilevelAnalysisGene}}
#' @examples
#' set.seed(1)
#' x=rnorm(10, mean = 0)

#' # p-value obtained from a one-sample t-test
#' t.test(x, mu=1, alternative = "less")$p.value
#' # p-value obtained from an intra-experiment analysis
#' IntraAnalysisClassic(x, func=t.test, mu=1, alternative = "less")
#' 
#' # p-value obtained from a one-sample wilcoxon test
#' wilcox.test(x, mu=1, alternative = "less")$p.value
#' # p-value obtained from an intra-experiment analysis
#' IntraAnalysisClassic(x, func=wilcox.test, mu=1, alternative = "less")
#' 
#' set.seed(1)
#' x=rnorm(20, mean=0); y=rnorm(20, mean=1)

#' # p-value obtained from a two-sample t-test
#' t.test(x,y,alternative="less")$p.value
#' # p-value obtained from an intra-experiment analysis
#' IntraAnalysisClassic(x, y, func=t.test, alternative = "less")

#' # p-value obtained from a two-sample wilcoxon test
#' wilcox.test(x,y,alternative="less")$p.value
#' # p-value obtained from an intra-experiment analysis
#' IntraAnalysisClassic(x, y, func=wilcox.test, alternative = "less")
#' 
#' @export
IntraAnalysisClassic <- function(x, y=NULL, splitSize=5, metaMethod=addCLT, func=t.test, p.value="p.value", ...) {
    if (splitSize < 3) {
        stop("splitSize should be at least 3")
    }
    
    ret=NULL
    
    l = splitS(x, splitSize)
    
    retList <- lapply(l, 
                      FUN=function(z, y, func) {
                          func(z , y, ...)[p.value]
                      }, y=y, func=func)
    
    metaMethod(unlist(retList))
}

#' @title Bi-level meta-analysis in conjunction with a classical 
#' hypothesis testing method
#' @description Perform a bi-level meta-analysis in conjunction with any of 
#' the classical hypothesis testing methods, such as t-test, Wilcoxon test, etc.
#' @param x a list of numeric vectors
#' @param y an optional list of numeric vectors
#' @param splitSize the minimum number of size in each split sample. 
#' splitSize should be at least 3. By default, splitSize=5
#' @param metaMethod the method used to combine p-values. This should be one 
#' of addCLT (additive method [1]), fishersMethod (Fisher's method [5]), 
#' stoufferMethod (Stouffer's method [6]), max (maxP method [7]), 
#' or min (minP method [8])
#' @param func the name of the hypothesis test. By default func=t.test
#' @param p.value the component that returns the p-value after performing the 
#' test provided by the \emph{func} parameter. For example, the function 
#' t-test returns the class "htest" where the component "p.value" is the 
#' p-value of the test. By default, p.value="p.value"
#' @param ... additional parameters for \emph{func}
#' @details This function performs a bi-level meta-analysis for the lists 
#' of samples [1]. It performs intra-experiment analyses to compare the 
#' vectors in x agains the corresponding vectors in y using the function 
#' \code{\link{IntraAnalysisClassic}} in conjunction with the test provided 
#' in \emph{func}. For example, it compares the first vector in x with the 
#' first vector in y, the second vector in x with the second vector in y, etc. 
#' When y is null, then the comparisons are reduced to one-sample tests. After 
#' these comparisons, we have a list of p-values, one for each comparision. 
#' The function then combines these p-values to obtain a single p-value 
#' using \emph{metaMethod}.
#' @return
#' the combined p-value
#' @author
#' Tin Nguyen and Sorin Draghici
#' @references
#' [1] T. Nguyen, R. Tagett, M. Donato, C. Mitrea, and S. Draghici. A novel 
#' bi-level meta-analysis approach -- applied to biological pathway analysis. 
#' Bioinformatics, 32(3):409-416, 2016.
#' 
#' @seealso \code{\link{IntraAnalysisClassic}}, \code{\link{IntraAnalysisGene}}, \code{\link{BilevelAnalysisGene}}
#' @examples
#' set.seed(1)
#' l1 = lapply(as.list(seq(3)),FUN=function (x) rnorm(n=10, mean=1))
#' l1
#' # one-sample t-test
#' lapply(l1, FUN=function(x) t.test(x, alternative="greater")$p.value)
#' # combining the p-values of one-sample t-tests:
#' addCLT(unlist(lapply(l1, FUN=function(x) t.test(x, alter="g")$p.value)))
#' #Bi-level meta-analysis
#' BilevelAnalysisClassic(x=l1, alternative="greater")
#' @export
BilevelAnalysisClassic <- function(x, y=NULL, splitSize=5, metaMethod=addCLT, func=t.test, p.value="p.value", ...) {
    if (splitSize < 3) {
        stop("splitSize should be at least 3")
    }
    
    metaMethod(sapply(seq(x),FUN=function(i) IntraAnalysisClassic(x[[i]], 
                y[[i]], splitSize, metaMethod, func, p.value, ...)))
}



#' @title Intra-experiment analysis of an expression dataset at the gene-level
#' @description perform an intra-experiment analysis in conjunction with 
#' the moderated t-test (limma package) for the purpose of differential 
#' expression analysis of a gene expression dataset
#' @param data a data frame where the rows are the gene IDs and the 
#' columns are the samples
#' @param group sample grouping. The elements of \emph{group} are 
#' either 'c' (control) or 'd' (disease). names(group) should be 
#' identical to colnames(data)
#' @param splitSize the minimum number of disease samples in each split 
#' dataset. splitSize should be at least 3. By default, splitSize=5
#' @param metaMethod the method used to combine p-values. This should be one 
#' of addCLT (additive method [1]), fishersMethod (Fisher's method [5]), 
#' stoufferMethod (Stouffer's method [6]), max (maxP method [7]), 
#' or min (minP method [8])
#' @details This function performs an intra-experiment analysis [1] for 
#' individual genes of the given dataset. The function first splits the 
#' dataset into smaller datasets, performs a moderated t-test (limma package) 
#' for the genes of the split datasets, 
#' and then combines the p-values for individual genes using \emph{metaMethod}
#' @return
#' A data frame (rownames are gene IDs) that consists of the 
#' following information:
#' \itemize{
#' \item \emph{logFC:} log foldchange (diseases versus controls)
#' \item \emph{pLimma:} p-value obtained from limma without spliting
#' \item \emph{pLimma.fdr:} FDR-corrected p-values of pLimma
#' \item \emph{pIntra:} p-value obtained from intra-experiment analysis
#' \item \emph{pIntra.fdr:} FDR-corrected p-values of pIntra
#' }
#' @author
#' Tin Nguyen and Sorin Draghici
#' @references
#' [1] T. Nguyen, R. Tagett, M. Donato, C. Mitrea, and S. Draghici. A novel 
#' bi-level meta-analysis approach -- applied to biological pathway analysis. 
#' Bioinformatics, 32(3):409-416, 2016.
#' 
#' @seealso \code{\link{BilevelAnalysisGene}}, \code{\link{IntraAnalysisClassic}}, \code{link{BilevelAnalysisClassic}}
#' @examples
#' data(GSE14924_CD4)
#' data = data_GSE14924_CD4
#' group=group_GSE14924_CD4
#' X = IntraAnalysisGene(data, group)
#' View(X)
#' 
#' @import limma
#' @export
IntraAnalysisGene <- function (data, group, splitSize=5, metaMethod=addCLT) {
    if (splitSize < 3) {
        stop("splitSize should be at least 3")
    }
        
    group <- group[order(group)]
    data <- data[,names(group)]
    
    diseases=paste(names(group)[which(group=="d")])
    
    l = splitS(diseases, splitSize)
    
    resList=lapply(l, 
                    FUN= function (sample, group, data) {
                        ret=rep(NA, nrow(data))
                        names(ret)=rownames(data)
                            
                        s=c(names(group)[which(group=="c")], sample)
                        g=group[s]
                        d <- data[,names(g)]
                                     
                        gset=ExpressionSet(as.matrix(d))
                        fl <- as.factor(g)
                        gset$description <- fl
                        design <- model.matrix(~ description + 0, gset)
                        colnames(design) <- levels(fl)
                        fit <- lmFit(gset, design)
                        cont.matrix <- makeContrasts(d-c, levels=design)
                        fit2 <- contrasts.fit(fit, cont.matrix)
                        fit2 <- eBayes(fit2)
                        tT <- topTable(fit2, adjust="none", sort.by="none", number=Inf)
                                 
                        ret[rownames(tT)] = tT$P.Value
                         }
                        , group=group, data=data)
    
    result=data.frame(row.names=rownames(data))
    result=cbind(result,do.call(cbind, resList))
    result$pIntra=apply(result, 1, FUN = metaMethod)
    
    gset=ExpressionSet(as.matrix(data))
    fl <- as.factor(group)
    gset$description <- fl
    design <- model.matrix(~ description + 0, gset)
    colnames(design) <- levels(fl)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(d-c, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    tT <- topTable(fit2, adjust="none", sort.by="none", number=Inf)

    FinalResult = data.frame(row.names=rownames(data))
    FinalResult$logFC = tT$logFC
    FinalResult$pLimma = tT$P.Value
    FinalResult$pLimma.fdr=p.adjust(FinalResult$pLimma,method = "fdr")
    FinalResult$pIntra = result$pIntra
    FinalResult$pIntra.fdr=p.adjust(FinalResult$pIntra,method = "fdr")
    FinalResult
}

#' @title Bi-level meta-analysis of multiple expression datasets at the 
#' gene-level
#' @description Perform a bi-level meta-analysis in conjunction with the 
#' moderate t-test (limma package) for the purpose of 
#' differential expression analysis of multiple gene expression datasets
#' @param dataList a list of datasets. Each dataset is a data frame where the 
#' rows are the gene IDs and the columns are the samples
#' @param groupList a list of vectors. Each vector represents the phenotypes 
#' of the corresponding dataset in dataList, 
#' which are either 'c' (control) or 'd' (disease). 
#' @param splitSize the minimum number of disease samples in each split 
#' dataset. splitSize should be at least 3. By default, splitSize=5
#' @param metaMethod the method used to combine p-values. This should be one 
#' of addCLT (additive method [1]), 
#' fishersMethod (Fisher's method [5]), stoufferMethod (Stouffer's method [6]), 
#' max (maxP method [7]), or min (minP method [8])
#' @details The bi-level framework combines the datasets at two levels: an 
#' intra- experiment analysis, and an inter-experiment analysis [1]. At the 
#' intra-experiment analysis, the framework splits a dataset into 
#' smaller datasets, performs a moderated t-test (limma package) on split 
#' datasets, and then combines p-values of individual genes 
#' using \emph{metaMethod}. At the inter-experiment analysis, the p-values 
#' obtained for each individual datasets are combined using \emph{metaMethod}
#' @return
#' A data frame containing the following components:
#' \itemize{
#' \item \emph{rownames:} gene IDs that are common in all datasets
#' \item \emph{pLimma:} the p-values obtained by combining pLimma 
#' values of individual datasets
#' \item \emph{pLimma.fdr:} FDR-corrected p-values of pLimma
#' \item \emph{pBilevel:} the p-values obtained from combining pIntra 
#' values of individual datasets
#' \item \emph{pBilevel.fdr:} FDR-corrected p-values of pBilevel
#' }
#' @author
#' Tin Nguyen and Sorin Draghici
#' @references
#' [1] T. Nguyen, R. Tagett, M. Donato, C. Mitrea, and S. Draghici. A novel 
#' bi-level meta-analysis approach -- applied to biological pathway analysis. 
#' Bioinformatics, 32(3):409-416, 2016.
#' 
#' @seealso \code{\link{BilevelAnalysisGene}}, \code{\link{IntraAnalysisClassic}}
#' @examples
#' dataSets=c("GSE14924_CD4", "GSE17054", "GSE57194", "GSE33223", "GSE42140", "GSE8023")
#' data(list=dataSets, package="BLMA")
#' dataList <- list()
#' groupList <- list()
#' for (i in 1:length(dataSets)) {
#'     dataset=dataSets[i]
#'     group <- get(paste("group_",dataset,sep=""))
#'     data=get(paste("data_",dataset,sep=""))
#'     dataList[[i]] = data
#'     groupList[[i]] = group
#' }
#' names(dataList)=names(groupList)=dataSets
#' 
#' Z=BilevelAnalysisGene(dataList = dataList, groupList = groupList)
#' View(Z)
#' @export
BilevelAnalysisGene <- function (dataList, groupList, splitSize=5, 
                                 metaMethod=addCLT) {
    retList = lapply(as.list(seq(length(dataList))), 
                FUN = function (i, dataList, groupList, metaMethod) {
                    cat("Working on dataset ", names(dataList)[[i]], ", " , ncol(dataList[[i]]), " samples \n", sep="")
                    
                    IntraAnalysisGene(dataList[[i]], groupList[[i]], splitSize=splitSize, metaMethod=metaMethod)
                }, dataList=dataList, groupList=groupList, metaMethod=metaMethod)
    
    names(retList)=names(dataList)
    
    genes=Reduce(intersect, lapply(dataList, FUN=rownames))
    pTableLimma=data.frame(row.names=genes)
    for (i in 1:length(dataList)) {
        pTableLimma[,i]=retList[[i]][genes,"pLimma"]
    }
    
    pTableIntra=data.frame(row.names=genes)
    for (i in 1:length(dataList)) {
        pTableIntra[,i]=retList[[i]][genes,"pIntra"]
    }
    
    Result=data.frame(row.names=genes, pLimma=apply(pTableLimma, 1, FUN = metaMethod), 
                        pBilevel=apply(pTableIntra, 1, FUN = metaMethod))
    Result$pLimma.fdr=p.adjust(Result$pLimma, method="fdr")
    Result$pBilevel.fdr=p.adjust(Result$pBilevel, method="fdr")
    Result <- Result[c(1,3,2,4)]
    Result=Result[order(Result$pBilevel),]
    Result
}

exampleGene <- function() {
    data(GSE14924_CD4)
    data = data_GSE14924_CD4
    group=group_GSE14924_CD4
    X = IntraAnalysisGene(data, group)
    View(X)
    
    dataSets=c("GSE14924_CD4", "GSE17054", "GSE57194", "GSE33223", "GSE42140", "GSE8023")
    data(list=dataset, package="BLMA")
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
    
    X=BilevelAnalysisGene(dataList,groupList)
    View(X$Bilevel)
    names(X$IntraLevel)
    View(X$IntraLevel[[1]])
}
