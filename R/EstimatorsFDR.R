
#########################################################################
# Author: Abbas Rahal
# Date: January 5, 2017
# EstimatorFDR is a function to compute the NFDR and the estimators
# of LFDR: CFDR and RFDR.
# INPUT: a vector of p-values.
# OUTPUT: NFDR and estimators of LFDR: CFDR and RFDR.
# Reference: Correcting false discovery rates for their bias toward
#            false positives. David R. Bickel, February 12, 2016
#########################################################################

EstimatorsFDR<- function (pvalue){

  if (is.vector(pvalue)==FALSE || is.numeric(pvalue)==FALSE) stop("\n The entry object should be a numeric vector \n")
  pvalue[!is.finite(pvalue)]<-1
  if (min(pvalue)<0 || max(pvalue)>1) stop("\n The p-value should be between 0 and 1 \n")
  d= length(pvalue)
  # Create a data frame that contains an index and pvalue
  index.pvalue<-data.frame(index=1:d, pvalue)
  # Order the data by pvalue
  index.pvalue<- index.pvalue[order(index.pvalue$pvalue),]
  # Intialiaze those vectors
  NFDR<- CFDR<- RFDR<- NULL

  # Compute NFDR and CFDR
  for (i in 1:d)
    {
        NFDR[i]<- index.pvalue$pvalue[i]*d/i
        if ( NFDR[i] >1 )  {NFDR[i]=1}

        S=0
        for (k in 1:i) { S = S + 1/(i-k+1) }

        CFDR[i] = S*NFDR[i]
        if ( CFDR[i] >1 )  {CFDR[i]=1}
    }

  # Compute RFDR
  for (i in 1:d)
    {
        RFDR[i]<-NFDR[round(i/(1-1/exp(1)))] # round the order of the p-value
        if ( round(i/(1-1/exp(1))) > d )  RFDR[i] = 1

     }

  # Add the NFDR, CFDR and RFDR to index.pvalue data frame
  index.pvalue<- data.frame(index.pvalue,NFDR=NFDR, CFDR=CFDR, RFDR=RFDR)

  # Reorder the data by the index
  index.pvalue<- index.pvalue[order(index.pvalue$index),]

  # Return
  return(list(NFDR=index.pvalue$NFDR, CFDR=index.pvalue$CFDR, RFDR=index.pvalue$RFDR))

  } # End of the function

