library(tidyverse)
library(Biostrings)
library(parallel)

calculateFop <- function(seq,codon_weights,ctg_ser=F)
{
  seq <- as.character(seq)
  seq <- unlist(str_extract_all(seq,pattern="[:alpha:]{3}"))
  tmp <- data.frame(Codon=seq)
  tmp <- tmp %>% left_join(codon_weights,by="Codon")
  if (ctg_ser)
  {
    tmp <- tmp %>% filter(!Codon %in% c("ATG","TGG","TAG","TAA","TGA","CTG"))
  } else {
    tmp <- tmp %>% filter(!Codon %in% c("ATG","TGG","TAG","TAA","TGA"))
  }
  occ <- tmp %>% count(Optimal) %>% column_to_rownames("Optimal")
  fop <- occ["O","n"]/(occ["O","n"] + occ["N","n"])
  return(fop)

}

fungi <- readLines("../all_fungi.txt")
ctg_ser <- readLines("../ser_fungi.txt")
mclapply(fungi,function(species){
  df <- read.csv(paste0("../CAI_weights_w_reference/",species))
  df <- df %>% group_by(AA) %>% 
        mutate(Minimum=min(Mean),Optimal=if_else(Mean==Minimum,"O","N"))  
  seq <- readDNAStringSet(paste0("/data2/Labella2019/Genomes/cds_codonW/",species))
  gene.ids <- names(seq)
  fop <- lapply(seq,calculateFop,codon_weights=df,ctg_ser = (species %in% ctg_ser))
  fop <- as.data.frame(t(fop %>% bind_rows())) %>% rownames_to_column(var="GeneId")
  fop <- fop %>% dplyr::rename(Fop=V1)
  write.csv(fop,paste0("/data2/Labella2019/Fop/",species),quote = F,row.names = F,col.names = T)
  },mc.cores=16)


# for (species in fungi)
# {
#   df <- read.csv(paste0("../CAI_weights/",species))
#   df <- df %>% group_by(AA) %>% 
#         mutate(Minimum=min(Mean),Optimal=if_else(Mean==Minimum,"O","N"))  
#   seq <- readDNAStringSet(paste0("/data2/Labella2019/Genomes/cds_codonW/",species))
#   gene.ids <- names(seq)
#   fop <- lapply(seq,calculateFop,codon_weights=df,ctg_ser = (species %in% ctg_ser))
#   fop <- as.data.frame(t(fop %>% bind_rows())) %>% rownames_to_column(var="GeneId")
#   fop <- fop %>% rename(Fop=V1)
#   write.csv(fop,paste0("../Fop/",species),quote = F,row.names = F,col.names = T)
# }
