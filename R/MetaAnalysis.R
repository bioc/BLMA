#' @title Fisher's method for meta-analysis
#' @description Combine independent p-values using the minus log product
#' 
#' @param x is an array of independent p-values
#' 
#' @details
#' 
#' Considering a set of \emph{m} independent significance tests, the resulted 
#' p-values are independent and uniformly distributed between \emph{0} and 
#' \emph{1} under the null hypothesis. Fisher's method uses the minus 
#' log product of the p-values as the summary statistic, which follows a 
#' chi-square distribution with \emph{2m} degrees of freedom. 
#' This chi-square distribution is used to calculate the combined p-value.
#' 
#' @return
#' 
#' combined p-value
#' 
#' @author
#' 
#' Tin Nguyen and Sorin Draghici
#' 
#' @references
#' 
#' [1] R. A. Fisher. Statistical methods for research workers. 
#' Oliver & Boyd, Edinburgh, 1925.
#' 
#' @seealso \code{\link{stoufferMethod}}, \code{\link{addCLT}}
#' 
#' @examples
#' 
#' x=rep(0,10)
#' fishersMethod(x)
#' 
#' x=runif(10)
#' fishersMethod(x)
#' 
#' @import stats
#' @export 
fishersMethod = function(x) {
    if(sum(is.na(x))>0) NA
    else pchisq(-2 * sum(log(x)), df=2*length(x), lower=FALSE)
}


#' @title Stouffer's method for meta-analysis
#' @description Combine independent studies using the sum of p-values 
#' transformed into standard normal variables
#' 
#' @param x is an array of independent p-values
#' 
#' @details
#' 
#' Considering a set of \emph{m} independent significance tests, the resulted 
#' p-values are independent and uniformly distributed between \emph{0} and 
#' \emph{1} under the null hypothesis. Stouffer's method is similar to 
#' Fisher's method (\link{fishersMethod}), with the difference is that it 
#' uses the sum of p-values transformed into standard normal variables 
#' instead of the log product.
#' 
#' @return
#' 
#' combined p-value
#' 
#' @author
#' 
#' Tin Nguyen and Sorin Draghici
#' 
#' @references
#' 
#' [1] S. Stouffer, E. Suchman, L. DeVinney, S. Star, and R. M. Williams. 
#' The American Soldier: Adjustment during army life, volume 1. 
#' Princeton University Press, Princeton, 1949.
#' 
#' @seealso \code{\link{fishersMethod}}, \code{\link{addCLT}}
#' 
#' @examples
#' 
#' x=rep(0,10)
#' stoufferMethod(x)
#' 
#' x=runif(10)
#' stoufferMethod(x)
#' 
#' @import stats
#' @export 
stoufferMethod = function(x) {
    if(sum(is.na(x))>0) NA
    else pnorm(sum(qnorm(x)) / sqrt(length(x)))
}

#x is a vector of p-values
IrwinHallDensity = function(x) {
    n=length(x)
    s=sum(x)
    1/factorial(n-1) * sum(sapply(0:floor(s), 
                function(k) (-1)^k * choose(n,k) * (s-k)^(n-1)))
}

#x is a vector of p-values
IrwinHallCumulative = function(x) {
    n=length(x)
    s=sum(x)
    1/factorial(n) * sum(sapply(0:floor(s), 
                function(k) (-1)^k * choose(n,k) * (s-k)^(n)))
}


additiveMethod = function(x) {
    #x is a vector of p-values
    n = length(x)
    if (n <= 20) {
        IrwinHallCumulative(x)
    } else {
        pnorm(sum(x),n/2,sqrt(n/12),lower=TRUE)
    }
}

#x is a vector of p-values
averageDensity = function(x) {
    n=length(x)
    a=mean(x)
    n/factorial(n-1) * sum(sapply(0:floor(n*a), 
                function(k) (-1)^k * choose(n,k) * (n*a-k)^(n-1)))
}

#x is a vector of p-values
averageCumulative = function(x) {
    n=length(x)
    a=mean(x)
    1/factorial(n) * sum(sapply(0:floor(n*a), 
                function(k) (-1)^k * choose(n,k) * (n*a-k)^(n)))
}




#' @title The additive method for meta-analysis
#' @description Combine independent studies using the average of p-values
#' 
#' @param x is an array of independent p-values
#' 
#' @details
#' 
#' This method is based on the fact that sum of independent uniform variables 
#' follow the Irwin-Hall distribution [1a,1b]. When the number of p-values 
#' is small (\emph{n<20}), the distribution of the average of p-values can 
#' be calculated using a linear transformation of the Irwin-Hall distribution.
#' When \emph{n} is large, the distribution is approximated using the 
#' Central Limit Theorem to avoid underflow/overflow problems [2,3,4,5]. 
#' 
#' @return
#' 
#' combined p-value
#' 
#' @author
#' 
#' Tin Nguyen and Sorin Draghici
#' 
#' @references
#' 
#' [1a] P. Hall. The distribution of means for samples of size n drawn from a 
#' population in which the variate takes values between 0 and 1, all such 
#' values being equally probable. Biometrika, 19(3-4):240-244, 1927.
#' 
#' [1b] J. O. Irwin. On the frequency distribution of the means of samples 
#' from a population having any law of frequency with finite moments, with 
#' special reference to Pearson's Type II. Biometrika, 19(3-4):225-239, 1927.
#' 
#' [2] T. Nguyen, R. Tagett, M. Donato, C. Mitrea, and S. Draghici. A novel 
#' bi-level meta-analysis approach -- applied to biological pathway analysis. 
#' Bioinformatics, 32(3):409-416, 2016.
#' 
#' [3] T. Nguyen, C. Mitrea, R. Tagett, and S. Draghici. DANUBE: Data-driven 
#' meta-ANalysis using UnBiased Empirical distributions -- applied to 
#' biological pathway analysis. Proceedings of the IEEE, PP(99):1-20, 2016.
#' 
#' [4] T. Nguyen, D. Diaz, R. Tagett, and S. Draghici. Overcoming the 
#' matched-sample bottleneck: an orthogonal approach to integrate omic data. 
#' Scientific Reports, 6:29251, 2016.
#' 
#' [5] T. Nguyen, D. Diaz, and S. Draghici. TOMAS: A novel TOpology-aware 
#' Meta-Analysis approach applied to System biology. In Proceedings of the 
#' 7th ACM International Conference on Bioin- formatics, Computational Biology, 
#' and Health Informatics, pages 13-22. ACM, 2016.
#' 
#' @seealso \code{\link{fishersMethod}}, \code{\link{stoufferMethod}}
#' 
#' @examples
#' 
#' x=rep(0,10)
#' addCLT(x)
#' 
#' x=runif(10)
#' addCLT(x)
#' 
#' @import stats
#' @export 
addCLT = function(x) {
    if(sum(is.na(x))>0) NA
    else {
        n = length(x)
        if (n <= 20) {
            averageCumulative(x)
        } else {
            pnorm(mean(x),1/2,sqrt(1/(12*n)),lower=TRUE)
        }  
    }
}


#' @title Spliting a vector of strings
#' @description Split a vector of strings into several vectors
#' @param x a vector of string
#' @param splitSize size of split vectors
#' @details
#' The function splits a vector of strings into several vectors
#' @return
#' a list of string vectors
#' @author
#' Tin Nguyen and Sorin Draghici
#' @examples 
#' x=paste("sample", seq(20))
#' splitS((x))
#' @export
splitS <- function(x, splitSize=5) {
    g=ceiling(seq(x)/splitSize)
    g[g==floor(length(x)/splitSize)+1]=floor(length(x)/splitSize)
    split(x, g)
}





