
################################################################################
# Author: Abbas Rahal
# Date: January 5, 2017
# BlendedLFDR is a function used to compute the blended estimator based on a
# benchmark estimator (usually nonlocale false discovery rate estimator (NFDR))
# and a set of estimators of local false discovery rate (LFDR).
# INPUT:
#      -Benchmark: is a vector.
#      -EstLFDR: is a matrix(recommended) or a data.frame.
# OUTPUT:
#      -a vector of Blended estimator of LFDR
# The Blended estimator is an estimator that minimizes the relative entropy of
# set of estimators of LFDR to Benchmark.
################################################################################

BlendedLFDR<- function(Benchmark, EstLFDR){

  EstLFDR[!is.finite(EstLFDR)]<-1
  Benchmark[!is.finite(Benchmark)]<-1
  if (is.vector(EstLFDR)==TRUE) stop("\n The number of columns of EstLFDR should be at least two columns \n")
  if (length(Benchmark)!= nrow(EstLFDR)) stop("\n The lengths of Benchmark and EstLFDR should be equal \n")
  if (min(Benchmark)<0 || max(Benchmark)>1) stop("\n The Benchmark should be between 0 and 1 \n")
  if (min(EstLFDR)<0 || max(EstLFDR)>1) stop("\n The LFDR should be between 0 and 1 \n")

  D<- matrix(NA, nrow = nrow(EstLFDR), ncol = ncol(EstLFDR))

  for (j in 1: ncol(EstLFDR)) {

  for (i in 1: nrow(EstLFDR))
  {
    # KL divergence between two Bernoulli Ber(p) and Ber(q)
    # D<- p*log(p/q) + (1-p)*log((1-p)/(1-q))

    if (EstLFDR[i,j]==0) { D[i,j]<- (1-EstLFDR[i,j])*log((1-EstLFDR[i,j])/(1-Benchmark[i]))
    } else if (EstLFDR[i,j]==1) { D[i,j]<- EstLFDR[i,j]*log(EstLFDR[i,j]/Benchmark[i])
    }else D[i,j]<- EstLFDR[i,j]*log(EstLFDR[i,j]/Benchmark[i]) + (1-EstLFDR[i,j])*log((1-EstLFDR[i,j])/(1- Benchmark[i]))

    }

  }
  # compute the minimun for each row
  minimum=apply( D[,1:ncol(EstLFDR)], 1, min)

  Blended<- NULL
  for (i in 1:nrow(EstLFDR)) {

    index<- which(D[i,]== minimum[i])
    Blended[i]<- EstLFDR[i,index[1]]

  }
return(list(Blended=Blended))

}
