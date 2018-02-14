\name{CorrectedFDR-package}
\alias{CorrectedFDR-package}
%\alias{CorrectedFDR}
\docType{package}
\title{
\packageTitle{CorrectedFDR}
}
\description{
\packageDescription{CorrectedFDR}
}
\details{

The DESCRIPTION file:
%\packageDESCRIPTION{CorrectedFDR}
%\packageIndices{CorrectedFDR}
\tabular{ll}{
Package: \tab CorrectedFDR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2018-02-08\cr
License: \tab GPL-3\cr
Depends: \tab R(>= 2.14.2)\cr
Suggests: \tab LFDR.MLE, LFDREmpiricalBayes, ProData\cr
URL: \tab https://davidbickel.com\cr
}
Two functions in \code{CorrectedFDR} package to compute the LFDR estimators. The function \code{EstimatorsFDR} computes the nonlocal false discovery rate (NFDR), the CFDR and the RFDR. The function \code{BlendedLFDR} uses a Benchmark of FDR, and other estimators of LFDR in order to get an estimate of LFDR.
}
\author{
\packageAuthor{CorrectedFDR}

Maintainer: \packageMaintainer{CorrectedFDR}
}
\references{
Bickel, D. R. (2016). Correcting false discovery rates for their bias toward false positives. Working paper. Retrieved from \url{http://hdl.handle.net/10393/34277}.

Bickel, D. R. (2015). Corrigendum to: Simple estimators of false discovery rates given as few as one or two p-values without strong parametric assumptions. Statistical Applications in Genetics and Molecular Biology, 2015, 14, 225.

Bickel, D. R. (2015). Blending Bayesian and frequentist methods according to the precision of prior information with applications to hypothesis testing. Statistical Methods and Applications, 24(4), pp. 523-546.

Bickel, D. R. (2013). Simple estimators of false discovery rates given as few as one or two p-values without strong parametric assumptions. Statistical Applications in Genetics and Molecular Biology, 2013, 12, 529-543.
 }
\keyword{ package }
%\seealso{
  %For more information on how to interpret the outputs, look at the supplementary file in the vignette directory,
  %"Using the CorrectedFDR Package."
%}

