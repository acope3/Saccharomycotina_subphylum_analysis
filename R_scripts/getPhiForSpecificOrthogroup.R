library(tidyverse)


orthomclToDataFrame <- function(orthomcl.string)
{
  orthogroup <- unlist(str_split(orthomcl.string,pattern = ":{0,1} "))
  orthoid <- orthogroup[1]
  orthogroup.genes <- orthogroup[-1]
  orthoid <- rep(orthogroup[1],length(orthogroup.genes))
  df <- data.frame(Orthogroup = orthoid,SeqID = orthogroup.genes)
  return(df)
}

all.species <- readLines("../Data/all_fungi.txt")

all.species.df <- data.frame(Species.full = all.species,
                             Species.cleaned = str_remove(all.species,"_{0,1}[0-9]*\\.max\\.cds"))

improved.fits.species <- readLines("../Data/improved_species_by_gene_expression_25_percent_2023_05_03.txt")


nonimproved.species <- setdiff(all.species,improved.fits.species)

k1.phi <- paste0("../Final_runs/Results_k_1/",nonimproved.species,"/run_2/Parameter_est/gene_expression.txt")
k2.phi <- paste0("../Final_runs/Results_k_2_selectionShared/",improved.fits.species,"/run_2/Parameter_est/gene_expression.txt")
improved.species <- str_extract(improved.fits.species,"[A-Za-z0-9_]*[a-z]+_[a-z]+(_JCM[0-9]{4}){0,1}")
nonimproved.species <- str_extract(nonimproved.species,"[A-Za-z0-9_]*[a-z]+_[a-z]+(_JCM[0-9]{4}){0,1}")


names(k1.phi) <- nonimproved.species
names(k2.phi) <- improved.species

phi.files <- c(k1.phi,k2.phi)

phi <- purrr::map(phi.files,read_csv,col_types = cols()) %>% bind_rows(.id="Species")

phi <- phi %>%
  mutate(GeneID=str_remove(GeneID,"-mRNA-1"))

seqid_index <- read_tsv("../orthomcl_output/343taxa_protein_IDs_index.txt",col_names = F)
seqid_index <- seqid_index %>% 
  separate(col = X2,sep = "@",into=c("Species","Gene")) %>%
  dplyr::select(-X1) %>% 
  mutate(GeneID=str_extract(X3,"(?<=gene=)\\S+")) %>%
  dplyr::select(GeneID,Species,Gene) %>% 
  mutate(GeneID=ifelse(is.na(GeneID),Gene,GeneID))


ortho_seqid <- read_delim("../orthomcl_output/orthomcl_SeqIDs_index.txt",delim = ": ",col_names=F) 
ortho_seqid <- ortho_seqid %>%
  separate(col = X2,sep = "@",into=c("Species","Gene")) %>%
  dplyr::rename(SeqID=X1)



orthogroups <- readLines("../orthomcl_output/orthomcl_clusters.txt")

orthogroups.df <- lapply(orthogroups,orthomclToDataFrame) %>%
  bind_rows()

orthogroups.df <- orthogroups.df %>% 
  left_join(ortho_seqid,by="SeqID") %>%
  left_join(seqid_index,by=c("Species","Gene"))

orthogroups.df <- orthogroups.df %>% 
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

target.ortho <- "OG1606"

tad2 <- orthogroups.df %>% 
  left_join(phi %>% dplyr::select(Species,GeneID,Mean),by=c("Species","GeneID")) %>% 
  filter(Orthogroup == target.ortho)

tad2.gene.count <- tad2 %>% 
  group_by(Species) %>% 
  summarize(Count=n()) %>% 
  arrange(desc(Count))

species.to.include <- tad2.gene.count %>%
  filter(Count == 1) %>%
  dplyr::select(Species) %>%
  deframe()

tad2.cleaned <- tad2 %>% 
  right_join(all.species.df,by=c("Species"="Species.cleaned")) %>%
  mutate(Mean = ifelse(Species %in% species.to.include,Mean,NA),
         Orthogroup = ifelse(is.na(Orthogroup),target.ortho,Orthogroup)) %>%
  distinct(Species,.keep_all = T) %>%
  dplyr::select(-Species) %>%
  dplyr::rename(Species = Species.full) %>%
  relocate(Species)
  

tad2.wide <- tad2.cleaned %>% 
  pivot_wider(id_cols = Orthogroup,names_from=Species,values_from=Mean)

og.wide <- tad2.cleaned %>%
  pivot_wider(id_cols = Orthogroup,names_from=Species,values_from=GeneID)

write_csv(tad2.wide,"../Data/msh2_phi.csv")
write_csv(og.wide,"../Data/shen_etal_msh2_orthogroup.csv")


