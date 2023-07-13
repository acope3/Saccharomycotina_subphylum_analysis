library(tidyverse)

all.species <- readLines("~/Labella2019_scripts/all_fungi.txt")
improved.species <- readLines("~/Labella2019_scripts/improved_species_by_gene_expression_2022_07_09.txt")
nonimproved.species <- setdiff(all.species,improved.species)

k1.phi <- paste0("/data2/Labella2019/Final_runs/Results_k_1/",nonimproved.species,"/run_2/Parameter_est/gene_expression.txt")
k2.phi <- paste0("/data2/Labella2019/Final_runs/Results_k_2_selectionShared/",improved.species,"/run_2/Parameter_est/gene_expression.txt")
improved.species <- str_extract(improved.species,"[A-Za-z0-9_]*[a-z]+_[a-z]+(_JCM[0-9]{4}){0,1}")
nonimproved.species <- str_extract(nonimproved.species,"[A-Za-z0-9_]*[a-z]+_[a-z]+(_JCM[0-9]{4}){0,1}")


names(k1.phi) <- nonimproved.species
names(k2.phi) <- improved.species

phi.files <- c(k1.phi,k2.phi)

phi <- purrr::map(phi.files,read_csv,col_types = cols()) %>% bind_rows(.id="Species")

seqid_index <- read_tsv("/data2/Labella2019/Orthologs_Shen_etal/orthomcl_output/343taxa_protein_IDs_index.txt",col_names = F)
seqid_index <- seqid_index %>% separate(col = X2,sep = "@",into=c("Species","Gene")) %>%
               dplyr::select(-X1) %>% 
               mutate(GeneID=str_extract(X3,"(?<=gene=)\\S+")) %>%
               dplyr::select(GeneID,Species,Gene) %>% 
               mutate(GeneID=ifelse(is.na(GeneID),Gene,GeneID))


cluster <- read_tsv("/data2/Labella2019/Orthologs_Shen_etal/orthomcl_output/343taxa_2408OGs_long_seqIDs.txt",col_names = F)
cluster <- cluster %>% 
  separate(col = X2,sep = "@",into=c("Species","Gene")) %>% 
  dplyr::rename(Orthogroup=X1)

cluster.seq <- cluster %>% inner_join(seqid_index,by=c("Species","Gene"))

cluster.seq <- cluster.seq %>% 
  mutate(Species = purrr::map_chr(Species,function(x){
  tmp <- unlist(str_split(x,"_"))
  if (length(tmp) > 2) {
    if (str_detect(tmp[1],"yH[MABP]{2}"))
    {
      tmp[2] <- tolower(tmp[2])
    } else{
      tmp[1] <- tolower(tmp[1])
    }
  } else{
    tmp[1] <- tolower(tmp[1]) 
  }
  paste(tmp,collapse="_")
}))



cluster.seq.filt <-  cluster.seq %>% 
  filter(Species %in% names(phi.files))

cluster.seq.filt.wide <-  cluster.seq.filt %>% 
                      dplyr::select(Orthogroup,Species,GeneID) %>%
                      pivot_wider(id_cols = Orthogroup,names_from="Species",values_from = "GeneID")
                      

phi <- phi %>%
  mutate(GeneID=str_remove(GeneID,"-mRNA-1"))

seqid_index_phi <- cluster.seq.filt %>%
               left_join(phi,by=c("Species","GeneID")) %>%
               dplyr::select(Orthogroup,Species,GeneID,Gene,Mean)

seqid_index_phi_wide <- seqid_index_phi  %>%
              dplyr::select(Orthogroup,Species,Mean) %>% 
              pivot_wider(id_cols = Orthogroup,names_from="Species",values_from = "Mean")
write_csv(seqid_index_phi_wide,"~/Labella2019_scripts/phi_shen_etal_ortholog_matrix_2022_07_10.csv")