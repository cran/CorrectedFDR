\name{EstimatorsFDR}
\alias{EstimatorsFDR}
\title{Estimators of Local False Discovery Rate (LFDR)}
\usage{
EstimatorsFDR(pvalue)
}
\arguments{
\item{pvalue}{Input numeric vector of pvalues.}
}
\value{
The output returns three lists. It returns the NFDR, CFDR, and RFDR estimators: \tabular{ll}{ NFDR\tab nonlocal FDR\cr CFDR\tab corrected FDR\cr RFDR\tab re-ranked FDR\cr}
}
\description{
EstimatorsFDR is an R function that computes the Nonlocal False Discovery Rate (NFDR) and the estimators of local false discovery rate: Corrected False discovery Rate (CFDR) and Re-ranked False Discovery rate (RFDR).
}
\details{
The input is a list of pvalues. The pvalues can be obtained for example by performing Student's t-test between two datasets. The two groups can be data from healthy and disease states. Let \eqn{i=1, 2, ..., N}, where \eqn{i} represents the \eqn{ith} feature (SNP or gene, for example). Then, for each \eqn{i}, the hypothesis indicator \eqn{A_i} can have two possible values.\cr \eqn{A_i=0}, if the \eqn{ith} null hypothesis is true, or \cr \eqn{A_i=1},  if the \eqn{ith} null hypothesis is not true, \cr
where the null hypothesis is defined by: the \eqn{ith} feature is unaffected by a treatment, unassociated with a disease, etc.
The values for each estimator (NFDR, CFDR, RFDR) indicate the probability that the null hypothesis of the \eqn{ith} feature is true (\eqn{A_i=0}) given the statistics \eqn{T_i}. The alternative hypothesis is true if \eqn{A_i=1}.
For example, in gene expression data analysis, if the null hypothesis is true, this would mean that the genes are not differentially expressed.
}
\examples{
#The examples below are from the "ProData" package.
#In order to use the "Prodata" input you would first need to install the ProData package.
#You will also need the function exprs in this package.
#First, make sure that the ProData package is properly installed:
#source("https://bioconductor.org/biocLite.R")
#biocLite("ProData")
#library(ProData)
#data("f45cbmk")
#q1<- quantile(as(exprs(f45cbmk[, pData(f45cbmk)$GROUP == "B"]), "numeric"), probs = 0.25)
#logish<- function(x){log(x + q1)}
#Vectors of proteins for 20 patients ER/PR-positive and Healthy
#Y<- logish(exprs(f45cbmk[, pData(f45cbmk)$GROUP == "B"])) # Control (Healthy)
#X.ER<- logish(exprs(f45cbmk[, pData(f45cbmk)$GROUP == "C"])) # Case ER/PR-positive
#pvalue<- NULL
#for (i in 1:nrow(X.ER))
#{
# t<-t.test(x=X.ER[i,], y=Y[i,], alternative = "two.sided")
# pvalue[i]<- t$p.value
#}

#The pvalues obtained from the t-test:
pvalue<- c(0.1981, 0.3794, 0.000001443, 0.02325, 0.03264, 0.07263, 0.02965, 0.8016, 0.8888,
          0.9133, 0.2971, 0.4573, 0.2815, 0.0007119, 0.5743, 0.927, 0.369, 0.8478, 0.38, 0.9904)
output<- EstimatorsFDR(pvalue)
#Three lists
output$NFDR
output$CFDR
output$RFDR
}
\references{
Bickel, D. R. (2016). Correcting false discovery rates for their bias toward false positives. Working paper. Retrieved from \url{http://hdl.handle.net/10393/34277}.

Bickel, D. R. (2015). Corrigendum to: Simple estimators of false discovery rates given as few as one or two p-values without strong parametric assumptions. Statistical Applications in Genetics and Molecular Biology, 2015, 14, 225.

Bickel, D. R. (2013). Simple estimators of false discovery rates given as few as one or two p-values without strong parametric assumptions. Statistical Applications in Genetics and Molecular Biology, 2013, 12, 529-543.
}
\author{
Code: Abbas Rahal. \cr Documentation: Anna Akpawu, Justin Chitpin and Abbas Rahal. \cr Maintainer: Abbas Rahal <arahal@uOttawa.ca>
}
\keyword{CFDR}
\keyword{False Discovery Rate}
\keyword{LFDR}
\keyword{NFDR}
\keyword{Null}
\keyword{RFDR}
\keyword{hypothesis}
