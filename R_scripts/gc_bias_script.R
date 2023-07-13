calcDeltaMFromGCBias <- function(gcBias, include.ref=TRUE){
  #' Calculates mutation bias terms $\Delta M$ based on GC content.
  #'
  #' @description Calculates $\Delta M$ based on GC bias.  Calculations
  #  assume freq A = freq T and freq C = freq G.
  #' @param gcBias: frequency of GC within a region or genome.
  #' @param include.ref: indicates whether reference codons should be included in results.
  #' Default value is TRUE
  #'
  #' @return Vector of $\Delta M$ values for each codon. Codon strings are used for entry names.
  #'
  #' @details Reference codons in AnaCoDa are the last alphabetical codon for each
  #' amino acid. By convention, the mutation and selection bias parameters, $\Delta M$ and $\Delta \eta$,
  #' are scaled so that $\Delta M = 0$ and $\Delta \eta = 0$ for each reference codon.
  #' In general the single codon amino acids and the stop codons (M, W, and X) are ignored.
#' These reference codons values are not used by `initializeParameterObject()` when setting
#' values of the `mutation.prior.mean` and `mutation.prior.sd`


  require("AnaCoDa")
  atBias  <-  1- gcBias
  
  ## calculate nt frequencies
  fG <- gcBias/2
  fC <- fG
  fA <- atBias/2
  fT <- fA
  
  fVec = c(fA,fC,fG,fT)
  names(fVec) <- c("A", "C", "G", "T")
  ## calculate M values
  ## note that sum fi = 1
  mVec <- - log(fVec)
  
  
  aaList  <- AnaCoDa::aminoAcids()
  ## Drop stop codons, M, and W
  aaList <- aaList[grep("X|M|W", aaList, invert = TRUE)]
  
  ## Create a vector for mutation terms
  M <- vector()
  deltaM <- vector()
  
  for(aa in aaList){
    ## Get all codons, including reference.
    ## I would expect focal = TRUE to be the correct syntax for this, but it's not
    codons <- AnaCoDa::AAToCodon(aa, focal = FALSE);
    nCodons = length(codons);
    refCodon <- codons[nCodons]
    ##print(c(aa, nCodons, codons))
    
    ## Calculate M for each codon
    for(codon in codons)
    {
      thirdNt <- substring(codon, 3,3)
      M[codon] <- mVec[thirdNt]
    }
    
    ## Scale relative to last codon
    for(codon in codons)
    {
      deltaM[codon] = M[codon]-M[refCodon]
    }
    
  }
  
  if(include.ref == FALSE)
  {
    deltaM <- removeRefCodonsFromVector(deltaM, aaList)
  }
  return(deltaM)
}



calcCodonFreqFromDeltaM <- function(deltaMVec){
#' Calculates expected codon frequencies based on $\Delta M$ based on GC content.
#'
#' @description Calculates expected codon frequencies based on $\Delta M$ values
#'  Does **not** assume freq A = freq T and freq C = freq G.
#'
#' @param deltaMVec: labeled vector containing $\Delta M$ values for each codon (label).
#'
#' @return Vector of of expected codon frequencies. Codon strings are used for entry names.
#' @details Reference codons in AnaCoDa are the last alphabetical codon for each
#' amino acid. By convention, the mutation and selection bias parameters, $\Delta M$ and $\Delta \eta$,
#' are scaled so that $\Delta M = 0$ and $\Delta \eta = 0$ for each reference codon.
#' In general the single codon amino acids and the stop codons (M, W, and X) are ignored.
#' These reference codons values are not used by `initializeParameterObject()` when setting
#' values of the `mutation.prior.mean` and `mutation.prior.sd`


## deltaMVec must include codon identity for its names
namesToKeep <- names(deltaMVec)

require("AnaCoDa")
aaList  <- AnaCoDa::aminoAcids()
## Drop stop codons, M, and W
aaList <- aaList[grep("X|M|W", aaList, invert = TRUE)]

codonFreq <- vector()
expTerms <- vector()


for(aa in aaList){
  codonNames <- names(codonFreq)
  codons <- AnaCoDa::AAToCodon(aa, focal = FALSE); ## Include
  reference codon
  nCodons  <- length(codons)
  refCodon <- codons[nCodons]
  
  ## Check if reference codon is included, if not add it to list
  if( is.na(deltaMVec[refCodon]) ){
    deltaMVec[refCodon] <- 0
  }
  
  for(codon in codons)
  {
    expTerms[codon] = exp(-deltaMVec[codon])
  }
  
  z <- sum(c(expTerms[codons]))
  codonFreq[codons] <- expTerms[codons]/z
  codonNames <- c(codonNames, codons)
  names(codonFreq) <- codonNames
}

## Keep only codons in original list
codonFreq <- codonFreq[names(codonFreq) %in% namesToKeep]



return(codonFreq)
}