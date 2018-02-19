#' @title Load KEGG pathways and names
#' @description Load KEGG pathways and names
#' @param organism organism code. Default value is "hsa" (human)
#' @param updateCache re-download KEGG pathways. Default value is FALSE
#' @return 
#' A list of the following components
#' \itemize{
#' \item \emph{kpg} a list of \code{\link{graphNEL}} objects 
#' encoding the pathway information. 
#' \item \emph{kpn} a named vector of pathway tiles. 
#' The names of the vector are the pathway KEGG IDs.
#' }
#' @author
#' Tin Nguyen and Sorin Draghici
#' @seealso \code{\link{keggPathwayGraphs}}, \code{\link{keggPathwayNames}}
#' @examples
#' x <- loadKEGGPathways()
#' @import ROntoTools 
#' @import graph
#' @export
loadKEGGPathways <- function (organism="hsa", updateCache=FALSE) {
    # pathwayFile <- paste(system.file("data", package="BLMA"), "/KEGGPathways_", organism, ".RData", sep = "") 
    # #pathwayFile <- paste("./inst/data/KEGGPathways_", organism, ".RData", sep = "") 
    # if (updateCache | !file.exists(pathwayFile)) {
    #     kpg <- keggPathwayGraphs(organism = organism, updateCache=updateCache)
    #     kpg <- setEdgeWeights(kpg)
    #     kpg <- setNodeWeights(kpg, defaultWeight = 1)
    #     suppressMessages(kpn <- keggPathwayNames(organism, updateCache = FALSE, verbose = FALSE))
    #     kpn <- kpn[names(kpg)]
    #     #dbInfo=keggInfo("pathway")
    #     #message(paste("Database info:\n", dbInfo, sep =""))
    #     save(kpg, kpn, file=pathwayFile)
    # } else {
    #     load(file=pathwayFile)
    #     #message(paste("Database info:\n", dbInfo, sep =""))
    # }

    kpg <- keggPathwayGraphs(organism = organism, updateCache=updateCache)
    kpg <- setEdgeWeights(kpg)
    kpg <- setNodeWeights(kpg, defaultWeight = 1)
    suppressMessages(kpn <- keggPathwayNames(organism, updateCache = FALSE, verbose = FALSE))
    kpn <- kpn[names(kpg)]

    list(kpg=kpg, kpn=kpn)
}

#' @import stats
pORACalc <- function(geneSet, DEGenes,measuredGenes, minSize =0) {
    #minSize is the minimum number of DE genes in the geneSet
    universe <- measuredGenes
    white <- DEGenes
    black <- setdiff(measuredGenes,DEGenes)
    sample <- intersect(measuredGenes,geneSet)
    sampleWhite <- intersect(sample,DEGenes)
    if (length(sampleWhite) < minSize) {
        res <- NA
    } else {
        res <- phyper(q=length(sampleWhite)-1,m=length(white),n=length(black), k=length(sample), lower.tail=FALSE)
    }
    res
}

#' @title Bi-level meta-analysis -- applied to pathway analysis 
#' 
#' @description Perform a bi-level meta-analysis conjunction with 
#' Impact Analysis to integrate multiple gene expression datasets
#' 
#' @param kpg list of pathway graphs as objects of type graph 
#' (e.g., \code{\link{graphNEL}})
#' @param kpn names of the pathways.
#' @param dataList a list of datasets to be combined. Each dataset is a data 
#' frame where the rows are the gene IDs and the columns are the samples.
#' @param groupList a list of vectors. Each vector represents the phenotypes 
#' of the corresponding dataset in dataList, which are either 'c' (control) 
#' or 'd' (disease). 
#' @param splitSize the minimum number of disease samples in each split 
#' dataset. splitSize should be at least 3. By default, splitSize=5
#' @param metaMethod the method used to combine p-values. This should be one 
#' of addCLT (additive method [1]), fisherMethod (Fisher's method [5]), 
#' stoufferMethod (Stouffer's method [6]), max (maxP method [7]), or 
#' min (minP method [8])
#' @param pCutoff cutoff p-value used to identify differentially 
#' expressed (DE) genes. This parameter is used only when the enrichment 
#' method is "ORA". By default, pCutoff=0.05 (five percent)
#' @param percent percentage of genes with highest foldchange to be considered 
#' as differentially expressed (DE). This parameter is used when the enrichment 
#' method is "ORA". By default percent=0.05 (five percent). Please note that 
#' only genes with p-value less than pCutoff will be considered
#' @param nboot number of bootstrap iterations. By default, nboot=200
#' @param seed seed. By default, seed=1.
#' @param mc.cores the number of cores to be used in parallel computing. 
#' By default, mc.cores=1
#' 
#' @details The bi-level framework combines the datasets at two levels: an 
#' intra-experiment analysis, and an inter-experiment analysis [1]. At the 
#' intra-level analysis, the framework splits a dataset into smaller datasets, 
#' performs pathway analysis for each split dataset using 
#' Impact Analysis [2,3], and then combines the results of these split datasets 
#' using \emph{metaMethod}. At the inter-level analysis, the results obtained 
#' for individual datasets are combined using \emph{metaMethod}
#' 
#' @return
#' A data frame (rownames are geneset/pathway IDs) that consists of the 
#' following information:
#' \itemize{
#' \item \emph{Name:} name/description of the corresponding pathway/geneset
#' \item Columns that include the pvalues obtained from the intra-experiment 
#' analysis of individual datasets
#' \item \emph{pBLMA:} p-value obtained from the inter-experiment 
#' analysis using addCLT
#' \item \emph{rBLMA:} ranking of the geneset/pathway using addCLT
#' \item \emph{pBLMA.fdr:} FDR-corrected p-values
#' }
#' 
#' @author
#' Tin Nguyen and Sorin Draghici
#' 
#' @references
#' [1] T. Nguyen, R. Tagett, M. Donato, C. Mitrea, and S. Draghici. A novel 
#' bi-level meta-analysis approach -- applied to biological pathway analysis. 
#' Bioinformatics, 32(3):409-416, 2016.
#' 
#' [2] A. L. Tarca, S. Draghici, P. Khatri, S. S. Hassan, P. Mittal, 
#' J.-s. Kim, C. J. Kim, J. P. Kusanovic, and R. Romero. A novel signaling 
#' pathway impact analysis. Bioinformatics, 25(1):75-82, 2009.
#' 
#' [3] S. Draghici, P. Khatri, A. L. Tarca, K. Amin, A. Done, C. Voichita, 
#' C. Georgescu, and R. Romero. A systems biology approach for pathway 
#' level analysis. Genome Research, 17(10):1537-1545, 2007.
#' 
#' [4] R. A. Fisher. Statistical methods for research workers. 
#' Oliver & Boyd, Edinburgh, 1925.
#' 
#' [5] S. Stouffer, E. Suchman, L. DeVinney, S. Star, and J. Williams, RM. 
#' The American Soldier: Adjustment during army life, volume 1. 
#' Princeton University Press, Princeton, 1949.
#' 
#' [6] L. H. C. Tippett. The methods of statistics. 
#' The Methods of Statistics, 1931.
#' 
#' [7] B. Wilkinson. A statistical consideration in psychological research. 
#' Psychological Bulletin, 48(2):156, 1951.
#' @seealso \code{\link{bilevelAnalysisGeneset}}, \code{\link{pe}}, \code{\link{phyper}}
#'
#' @examples
#' # load KEGG pathways
#' x <- loadKEGGPathways()  
#' 
#' # load example data
#' dataSets <- c("GSE17054", "GSE57194", "GSE33223", "GSE42140")
#' data(list=dataSets, package="BLMA")
#' names(dataSets) <- dataSets
#' dataList <- lapply(dataSets, function(dataset) get(paste0("data_", dataset)))
#' groupList <- lapply(dataSets, function(dataset) get(paste0("group_", dataset)))
#' 
#' IAComb <- bilevelAnalysisPathway(x$kpg, x$kpn, dataList, groupList)
#' head(IAComb[, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")])
#' 
#' @import ROntoTools 
#' @import graph
#' @import limma
#' @import parallel
#' @import Biobase
#' @import utils
#' @import stats
#' @export
bilevelAnalysisPathway <- function (kpg, kpn, dataList, groupList, splitSize=5, metaMethod=addCLT, pCutoff=0.05, percent=0.05, mc.cores=1, nboot=200, seed=1) {
    if (splitSize < 3) {
        stop("splitSize should be at least 3")
    }
    
    allGenes <- unique(unlist(lapply(kpg,FUN=function(x){return (x@nodes)})))
    
    Results <- list()
    for (i in 1:length(dataList)) {
        cat("Working on dataset ", names(dataList)[[i]], ", " , ncol(dataList[[i]]), " samples \n", sep="")
        
        data <- dataList[[i]]
        group <- groupList[[i]]
        
        data <- data[rownames(data)%in%allGenes,]
        
        # make sure that we are not dividing by very small number
        if (min(data) < 2) {
            data <- data - min(data) + 2
        }
        
        
        diseases <- names(group)[which(group=="d")]
        
        l <- splitS(diseases, splitSize)

        resList <- mclapply(l, 
                    FUN= function (sample, group, data, kpg) {
                            ret <- rep(1, length(kpg))
                            names(ret) <- names(kpg)
                                 
                            message(paste(sample, collapse=", "))
                                 
                            s <- c(names(group)[which(group=="c")], sample)
                            g <- group[s]
                            d <- data[,names(g)]
                                 
                            gset <- ExpressionSet(as.matrix(d))
                            fl <- as.factor(g)
                            gset$description <- fl
                            design <- model.matrix(~ description + 0, gset)
                            colnames(design) <- levels(fl)
                            fit <- lmFit(gset, design)
                            cont.matrix <- makeContrasts(d-c, levels=design)
                            fit2 <- contrasts.fit(fit, cont.matrix)
                            fit2 <- eBayes(fit2)
                            tT <- topTable(fit2, adjust="none", sort.by="logFC", number=nrow(d)*percent, p.value=pCutoff)
                                 
                            measuredGenes <- rownames(d)
                            geneFC <- tT$logFC; names(geneFC) <- rownames(tT)
                                 
                            a <- capture.output(
                                peRes <- pe(x=geneFC,graphs=kpg, ref=measuredGenes, 
                                        nboot=nboot, seed=seed, verbose=FALSE))
                            res <- Summary(peRes)
                                 
                            if (length(geneFC)<length(measuredGenes)) {
                                ret[rownames(res)] <- as.numeric(res$pComb)
                            } else ret[rownames(res)] <- as.numeric(res$pAcc)
                                 
                            ret
                        }
                        , group=group, data=data, kpg=kpg, mc.cores=mc.cores)
        result <- data.frame(row.names=names(kpg), Names=kpn[names(kpg)])
        result <- cbind(result,do.call(cbind, resList))
        

        if ((1+length(diseases)/splitSize)>=3) {
            result$pIntra <- apply(result[,2:(1+length(diseases)/splitSize)], 1, FUN = metaMethod)
        } else {
            result$pIntra <- result[,2]
        }
        
        Results[[i]] <- result
    }
    
    r <- names(kpn)
    FinalResults <- data.frame(row.names=r, Name=kpn[r])
    for (i in 1:length(dataList)) {
        FinalResults[r, i+1] <- Results[[i]][r,"pIntra"]
    }
    colnames(FinalResults) <- c("Name", names(dataList))
    
    if (length(dataList)>1) {
        FinalResults$pBLMA <- apply(FinalResults[,2:(length(dataList)+1)], 1, FUN = metaMethod) 
    } else {
        FinalResults$pBLMA <- FinalResults[,2]
    }
    FinalResults$rBLMA[order(FinalResults$pBLMA)] <- seq(length(kpg))
    FinalResults$pBLMA.fdr <- p.adjust(FinalResults$pBLMA, method="fdr")
    
    FinalResults <- FinalResults[order(FinalResults$pBLMA),]
    FinalResults
} 

#' @title Bi-level meta-analysis -- applied to geneset enrichment analysis 
#' @description Perform a bi-level meta-analysis in conjunction with 
#' geneset enrichment methods (ORA/GSA/PADOG) to integrate 
#' multiple gene expression datasets.
#' @param gslist a list of gene sets.
#' @param gs.names names of the gene sets.
#' @param dataList a list of datasets to be combined. Each dataset is a 
#' data frame where the rows are the gene IDs and the columns are the samples.
#' @param groupList a list of vectors. Each vector represents the phenotypes 
#' of the corresponding dataset in dataList. 
#' The elements of each vector are either 'c' (control) or 'd' (disease).
#' @param splitSize the minimum number of disease samples in each split 
#' dataset. splitSize should be at least 3. By default, splitSize=5
#' @param metaMethod the method used to combine p-values. This should be one 
#' of addCLT (additive method [1]), fisherMethod (Fisher's method [5]), 
#' stoufferMethod (Stouffer's method [6]), max (maxP method [7]), 
#' or min (minP method [8])
#' @param enrichment the method used for enrichment analysis. This should be 
#' one of "ORA", "GSA", or "PADOG". By default, enrichment is set to "ORA".
#' @param pCutoff cutoff p-value used to identify differentially 
#' expressed (DE) genes. This parameter is used only when the enrichment 
#' method is "ORA". By default, pCutoff=0.05 (five percent)
#' @param percent percentage of genes with highest foldchange to be considered 
#' as differentially expressed (DE). This parameter is used when the 
#' enrichment method is "ORA". By default percent=0.05 (five percent). Please 
#' note that only genes with p-value less than pCutoff will be considered
#' @param mc.cores the number of cores to be used in parallel computing. 
#' By default, mc.cores=1
#' @param ... additional parameters of the GSA/PADOG functions
#' @details The bi-level framework combines the datasets at two levels: 
#' an intra- experiment analysis, and an inter-experiment analysis [1]. At the 
#' intra-level analysis, the framework splits a dataset into smaller datasets, 
#' performs enrichment analysis for each split dataset (using ORA [2], 
#' GSA [3], or PADOG [4]), and then combines the results of these split 
#' datasets using \emph{metaMethod}. At the inter-level analysis, the results 
#' obtained for individual datasets are combined using \emph{metaMethod}
#' @return
#' A data frame (rownames are geneset/pathway IDs) that consists of the 
#' following information:
#' \itemize{
#' \item \emph{Name:} name/description of the corresponding pathway/geneset
#' \item Columns that include the pvalues obtained from the intra-experiment 
#' analysis of individual datasets
#' \item \emph{pBLMA:} p-value obtained from the inter-experiment analysis 
#' using addCLT
#' \item \emph{rBLMA:} ranking of the geneset/pathway using addCLT
#' \item \emph{pBLMA.fdr:} FDR-corrected p-values
#' }
#' @author
#' Tin Nguyen and Sorin Draghici
#' @references
#' [1] T. Nguyen, R. Tagett, M. Donato, C. Mitrea, and S. Draghici. A novel 
#' bi-level meta-analysis approach -- applied to biological pathway analysis. 
#' Bioinformatics, 32(3):409-416, 2016.
#' 
#' [2] S. Draghici, P. Khatri, R. P. Martin, G. C. Ostermeier, and 
#' S. A. Krawetz. Global functional profiling of gene expression. 
#' Genomics, 81(2):98-104, 2003.
#' 
#' [3] B. Efron and R. Tibshirani. On testing the significance of sets of 
#' genes. The Annals of Applied Statistics, 1(1):107-129, 2007.
#' 
#' [4] A. L. Tarca, S. Draghici, G. Bhatti, and R. Romero. Down-weighting 
#' overlapping genes improves gene set analysis. 
#' BMC Bioinformatics, 13(1):136, 2012.
#' 
#' [5] R. A. Fisher. Statistical methods for research workers. 
#' Oliver & Boyd, Edinburgh, 1925.
#' 
#' [6] S. Stouffer, E. Suchman, L. DeVinney, S. Star, and J. Williams, RM. 
#' The American Soldier: Adjustment during army life, volume 1. 
#' Princeton University Press, Princeton, 1949.
#' 
#' [7] L. H. C. Tippett. The methods of statistics. 
#' The Methods of Statistics, 1931.
#' 
#' [8] B. Wilkinson. A statistical consideration in psychological research. 
#' Psychological Bulletin, 48(2):156, 1951.
#' @seealso \code{\link{bilevelAnalysisPathway}}, \code{\link{phyper}}, \code{\link{GSA}}, \code{\link{padog}}
#' @examples 
#' # load KEGG pathways and create gene sets
#' x <- loadKEGGPathways()
#' gslist <- lapply(x$kpg,FUN=function(y){return (nodes(y));})
#' gs.names <- x$kpn[names(gslist)]
#' 
#' # load example data
#' dataSets <- c("GSE17054", "GSE57194", "GSE33223", "GSE42140")
#' data(list=dataSets, package="BLMA")
#' names(dataSets) <- dataSets
#' dataList <- lapply(dataSets, function(dataset) get(paste0("data_", dataset)))
#' groupList <- lapply(dataSets, function(dataset) get(paste0("group_", dataset)))
#' # perform bi-level meta-analysis in conjunction with ORA
#' ORAComb <- bilevelAnalysisGeneset(gslist, gs.names, dataList, groupList, enrichment = "ORA")
#' head(ORAComb[, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")])
#' 
#' # perform bi-level meta-analysis in conjunction with GSA
#' GSAComb <- bilevelAnalysisGeneset(gslist, gs.names, dataList, groupList, enrichment = "GSA", nperms = 200, random.seed = 1)
#' head(GSAComb[, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")])
#' 
#' # perform bi-level meta-analysi in conjunction with PADOG
#' set.seed(1)
#' PADOGComb <- bilevelAnalysisGeneset(gslist, gs.names, dataList, groupList, enrichment = "PADOG", NI=200)
#' head(PADOGComb[, c("Name", "pBLMA", "pBLMA.fdr", "rBLMA")])
#' 
#' @import GSA
#' @import PADOG
#' @import limma
#' @import stats
#' @import parallel
#' @import Biobase
#' @import utils
#' @export
bilevelAnalysisGeneset <- function (gslist, gs.names, dataList, groupList, splitSize=5, metaMethod=addCLT, enrichment = "ORA", pCutoff=0.05, percent=0.05, mc.cores=1, ...) {  
    if (splitSize < 3) {
        stop("splitSize should be at least 3")
    }
    if (!is.element(enrichment, c("GSA","PADOG","ORA"))) {
        stop("enrichment must be GSA, PADOG, or ORA")
    }
    
    allGenes <- unique(unlist(gslist))
    
    Results <- list()
    for (i in 1:length(dataList)) {
        cat("Working on dataset ", names(dataList)[[i]], ", " , ncol(dataList[[i]]), " samples \n", sep="")
        
        data <- dataList[[i]]
        group <- groupList[[i]]
        
        data <- data[rownames(data)%in%allGenes,]
        group <- sort(group)
        data <- data[,names(group)]
        
        diseases <- paste(names(group)[which(group=="d")])
        
        l <- splitS(diseases, splitSize)
        
        resList <- mclapply(l, 
            FUN = function (sample, group, data, gslist) {
                ret <- rep(NA, length(gslist))
                names(ret) <- names(gslist)
                
                message(paste(sample, collapse=", "))
                
                s <- c(names(group)[which(group=="c")], sample)
                g <- group[s]
                d <- data[,names(g)]
                
                if (enrichment == "GSA") {
                    x <- as.character(g); x[x=="c"] <- 1; x[x=="d"] <- 2; x <- as.numeric(x)
                    a <- capture.output(resgsa <- GSA(x = as.matrix(d), y = x, 
                                 genesets = gslist, 
                                 genenames = rownames(data),
                                 resp.type = "Two class unpaired", 
                                 ...))
                    
                    ret <- 2 * apply(cbind(resgsa$pvalues.lo, 
                                          resgsa$pvalues.hi), 1, min)
                } else if (enrichment == "PADOG") {
                    a <- capture.output(respadog <- padog(
                        esetm=as.matrix(d),
                        group=g,
                        gslist=gslist,
                        gs.names=gs.names,
                        ...))
                    
                    ret[rownames(respadog)] <- as.numeric(paste(respadog[,"Ppadog"]))
                } else if (enrichment == "ORA") {
                          
                    gset <- ExpressionSet(as.matrix(d))
                    fl <- as.factor(g)
                    gset$description <- fl
                    design <- model.matrix(~ description + 0, gset)
                    colnames(design) <- levels(fl)
                    fit <- lmFit(gset, design)
                    cont.matrix <- makeContrasts(d-c, levels=design)
                    fit2 <- contrasts.fit(fit, cont.matrix)
                    fit2 <- eBayes(fit2)
                    tT <- topTable(fit2, adjust="none", sort.by="logFC", number=nrow(d)*percent, p.value=pCutoff)
                    
                    measuredGenes <- rownames(d)
                    DEGenes <- rownames(tT[tT$p])
                    
                    resORA <- unlist(lapply(gslist, FUN=pORACalc, DEGenes = DEGenes, measuredGenes = measuredGenes))
                    
                    ret[names(resORA)] <- resORA
                }
            
                ret
            }
            , group=group, data=data, gslist=gslist, mc.cores=mc.cores)
        
        result <- data.frame(row.names=names(gslist), Names=gs.names[names(gslist)])
        result <- cbind(result,do.call(cbind, resList))
        
        if ((1+length(diseases)/splitSize)>=3) {
            result$pIntra <- apply(result[,2:(1+length(diseases)/splitSize)], 1, FUN = metaMethod)
        } else {
            result$pIntra <- result[,2]
        }
        Results[[i]] <- result
    }
    
    r <- names(gslist)
    FinalResults <- data.frame(row.names=r, Name=gs.names[r])
    for (i in 1:length(dataList)) {
        FinalResults[r, i+1] <- Results[[i]][r,"pIntra"]
    }
    colnames(FinalResults) <- c("Name", names(dataList))
    
    if (length(dataList)>1) {
        FinalResults$pBLMA <- apply(FinalResults[,2:(length(dataList)+1)], 1, FUN = metaMethod) 
    } else {
        FinalResults$pBLMA <- FinalResults[,2]
    }
    FinalResults$rBLMA[order(FinalResults$pBLMA)] <- seq(length(gslist))
    FinalResults$pBLMA.fdr <- p.adjust(FinalResults$pBLMA, method="fdr")
    
    FinalResults <- FinalResults[order(FinalResults$pBLMA),]
    FinalResults
}





#' Gene expression dataset GSE17054 from Majeti et al.
#' @name GSE17054
#' @aliases 
#' data_GSE17054
#' group_GSE17054
#' @docType data
#' @description This dataset consists of 5 acute myeloid leukemia and 
#' 4 control samples. The data frame data_GSE17054 includes the expression data 
#' while the vector group_GSE17054 includes the grouping information.
#' @references 
#' Majeti et al. Dysregulated gene expression networks in human acute 
#' myelogenous leukemia stem cells. Proceedings of the National Academy of 
#' Sciences, 106(9):3396-3401, 2009.
#' @source 
#' Obtained from \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17054}
#' 
#' @usage data(GSE17054)
#' @format data_GSE17054 is a data frame with 4738 rows and 9 columns. 
#' The rows represent the genes and the columns represent the samples.
#' 
#' group_GSE17054 is a vector that represents the sample grouping for data_GSE17054. 
#' The elements of \emph{group_GSE17054} are either 'c' (control) or 'd' (disease). 
#' @keywords dataset
NULL



#' Gene expression dataset GSE57194 from Abdul-Nabi et al.
#' @name GSE57194
#' @aliases 
#' data_GSE57194
#' group_GSE57194
#' @docType data
#' @description This dataset consists of 6 acute myeloid leukemia and 
#' 6 control samples. The data frame data_GSE57194 includes the expression data 
#' while the vector group_GSE57194 includes the grouping information.
#' @references 
#' Abdul-Nabi et al. In vitro transformation of primary human CD34+ cells by 
#' AML fusion oncogenes: early gene expression profiling reveals possible drug 
#' target in AML. PLoS One, 5(8):e12464, 2010.
#' @source
#' Obtained from \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57194}
#' @usage data(GSE57194)
#' @format data_GSE57194 is a data frame with 4114 rows and 12 columns. 
#' The rows represent the genes and the columns represent the samples.
#' 
#' group_GSE57194 is a vector that represents the sample grouping for data_GSE57194. 
#' The elements of \emph{group_GSE57194} are either 'c' (control) or 'd' (disease). 
#' @keywords dataset
NULL



#' Gene expression dataset GSE33223 from Bacher et al.
#' @name GSE33223
#' @aliases 
#' data_GSE33223
#' group_GSE33223
#' @docType data
#' @description This dataset consists of 20 acute myeloid leukemia and 
#' 10 control samples. The data frame data_GSE33223 includes the expression data 
#' while the vector group_GSE33223 includes the grouping information.
#' @references 
#' Bacher et al. Multilineage dysplasia does not influence prognosis in 
#' CEBPA-mutated AML, supporting the WHO proposal to classify these patients 
#' as a unique entity. Blood, 119(20):4719-22, 2012.
#' @source 
#' Obtained from \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33223}
#' @usage data(GSE33223)
#' @format data_GSE33223 is a data frame with 4114 rows and 30 columns. 
#' The rows represent the genes and the columns represent the samples.
#' 
#' group_GSE33223 is a vector that represents the sample grouping for 
#' data_GSE33223. The elements of \emph{group_GSE33223} are either 
#' 'c' (control) or 'd' (disease). 
#' @keywords dataset
NULL


#' The gene expression dataset GSE42140 obtained from Gene Expression Omnibus
#' @name GSE42140
#' @aliases 
#' data_GSE42140
#' group_GSE42140
#' @docType data
#' @description This dataset consists of 26 acute myeloid leukemia and 
#' 5 control samples. The data frame data_GSE42140 includes the expression data 
#' while the vector group_GSE42140 includes the grouping information.
#' @references \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42140}
#' @usage data(GSE42140)
#' @format data_GSE42140 is a data frame with 4114 rows and 31 columns. 
#' The rows represent the genes and the columns represent the samples.
#' 
#' group_GSE42140 is a vector that represents the sample grouping for data_GSE42140. 
#' The elements of \emph{group_GSE42140} are either 'c' (control) or 'd' (disease). 
#' @keywords dataset
NULL
