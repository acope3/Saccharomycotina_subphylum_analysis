library(AnaCoDa,lib="~/R_dev")
library(argparse)
library(tidyverse)

getCAIweights <- function(referenceGenome,default.weight=0.5)
{
  aa.vec <- aminoAcids()
  aa.vec <- aa.vec[-length(aa.vec)]
  wi.list <- vector(mode = "list", length = length(aa.vec))
  codon.names <- NULL
  for(aa in aa.vec)
  {
  	print(aa)
    codon.names <- c(codon.names, AAToCodon(aa))
    ## create reference table for each codon and gene
    codonCountForAA.ref <- getCodonCountsForAA(aa, genome = referenceGenome)
    fi <- colSums(codonCountForAA.ref)
    fi.max <- max(fi)
    wi.list[[aa]] <- fi / fi.max
  }
  
  wi.vec <- unlist(wi.list)
  names(wi.vec) <- codon.names
  wi.vec[wi.vec == 0.0] <- default.weight
  return(wi.vec)
}

makeRelative <- function(x)
{	
	x$Mean <- unlist(x[nrow(x),"Mean"]) - x$Mean 
	return(x)
}

parser <- ArgumentParser()
parser$add_argument("-i","--input",help="FASTA file with protein-coding sequences",type="character",default="./")
parser$add_argument("--codon_table",help="Codon Table",type="integer",default=1)

args <- parser$parse_args()
input <- args$input
codon_table <- args$codon_table

ribo <- initializeGenomeObject(paste0("RiboProt/",input),codon_table = codon_table)
cai <- getCAIweights(ribo)
cai <- data.frame(AA=unlist(lapply(names(cai),codonToAA)),Codon=names(cai),Mean=unname(cai),stringsAsFactors = F)
print(cai)
cai <- cai %>% filter(!AA %in% c("W","X","M","J"))
cai <- cai %>% group_by(AA) %>% group_map(~ makeRelative(.x)) %>% bind_rows()
cai["AA"] <- unlist(lapply(cai$Codon,codonToAA))
cai <- cai[c("AA","Codon","Mean")]
#cai <- cai %>% filter(Mean!=0)
#write.csv(cai,paste0("CAI_weights_w_reference/",input),sep=",",row.names = F,col.names=T,quote=F)
