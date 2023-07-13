library(Biostrings)
library(tidyverse)


files <- list.files("../Expression/",pattern=".max.cds",recursive = F,)
for (f in files)
{
  fasta <- readDNAStringSet(file.path("/data2/Labella2019/Genomes/cds_cleaned",f))
  seq.names <- unlist(purrr::map(names(fasta) %>% str_split(pattern=" "),1))
  exp <- read_csv(file.path("../","Expression",f))
  colnames(exp)[1] <- "GeneID"
  exp <- exp %>% column_to_rownames("GeneID")
  exp <- exp[seq.names,,drop=F]
  exp <- exp %>%  mutate(
    across(everything(), ~replace_na(.x, 0))
  )
  exp <- exp %>% rownames_to_column(var="GeneID")
  exp <- exp %>% mutate(GeneID=seq.names)
  write_csv(exp,file = file.path("../","Expression_updated",f))
}