library(AnaCoDa,lib.loc="~/R_dev")
library(tidyverse)

splitSeq <- function(gene)
{
  seq <- gene$seq
  seq <- unlist(strsplit(seq, ""))
  seq <- paste(seq[c(T,F,F)], seq[c(F,T,F)], seq[c(F,F,T)], sep="")
  seq <- seq[which(seq %in% codons())]
  aa <- sapply(seq,codonToAA)
  seq.df <- data.frame(Codon=seq,AA=unname(aa))
  seq.df <- seq.df %>% filter(!AA %in% c("M","W","J","X"))
  return(seq.df)
}

candida.deta <- read_csv("../Labella2019_scripts/Results_k_2_selectionShared/hyphopichia_burtonii.max.cds/restart_5/Parameter_est/1_Selection.csv")

candida.deta <- candida.deta %>% separate(Codon,sep=1:3,into=c("NT1","NT2","NT3"),remove = F)

candida.deta <- candida.deta %>% group_by(AA) %>% dplyr::slice(which.min(Mean))

aa.gc <- unlist(candida.deta %>% filter(NT3 == "C" | NT3=="G") %>% dplyr::select(AA))
aa.at <- unlist(candida.deta %>% filter(NT3 == "A" | NT3=="T") %>% dplyr::select(AA))

genome <- initializeGenomeObject("/data2/Labella2019/Genomes/cds_cleaned/hyphopichia_burtonii.max.cds",codon_table = 12)
size <- length(genome)
genes <- lapply(1:size,function(x){genome$getGeneByIndex(x,F)})


seq.df <- lapply(genes,splitSeq)



genome$clear()
for (i in 1:size)
{
  tmp <- seq.df[[i]]
  tmp.gc <- tmp %>% filter(AA %in% aa.gc)
  seq.gc <- paste0(tmp.gc$Codon,collapse="")
  genes[[i]]$seq <- seq.gc
  genome$addGene(genes[[i]],F)
}

genome$writeFasta("hyphopichia_burtonii_gc_favored.fasta",F)
genome$clear()
for (i in 1:size)
{
  tmp <- seq.df[[i]]
  tmp.at <- tmp %>% filter(AA %in% aa.at)
  seq.at <- paste0(tmp.at$Codon,collapse="")
  genes[[i]]$seq <- seq.at
  genome$addGene(genes[[i]],F)
}

genome$writeFasta("hyphopichia_burtonii_at_favored.fasta",F)
genome$clear()