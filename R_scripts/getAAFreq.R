library(AnaCoDa)
library(tidyverse)


all.genomes <- readLines("../Data/all_fungi.txt")
ser.genomes <- readLines("../Data/ser_fungi.txt")
names(all.genomes) <- all.genomes

aa <- aminoAcids()
aa <- aa[which(!aa %in% c("M","X","W","J"))]
names(aa) <- aa


aa.df <- purrr::map(all.genomes,function(species){
  genome <- initializeGenomeObject(file.path("../Data/Genomes/cds_cleaned/",species))
  purrr::map(aa,function(a){
    codon.counts <- getCodonCountsForAA(a,genome)
    if (a == "L" & species %in% ser.genomes)
    {
      aa.c <- AAToCodon(a,F)
      codon.counts <- codon.counts[,which(aa.c != "CTG")]
    }
    sum(codon.counts)
  }) %>% bind_cols()
}) %>% bind_rows(.id="Species")

# aa.df <-aa.df %>% 
#   mutate(across(where(is.numeric),~./rowSums(.,na.rm)))

write_tsv(aa.df,"../Post_analysis/amino_acid_frequency_per_species.tsv")
