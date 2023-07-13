require(tAI)
library(tidyverse)
library(AnaCoDa)
library(ape)
library(parallel)

cleanTree <- function(tree,target.species)
{
  tips.to.drop <- tree$tip.label[which(!(tree$tip.label %in% target.species))]
  cleaned <- drop.tip(tree,tips.to.drop)
  return(cleaned)
}

cleanColumnNames <- function(col_name,list_of_option)
{
  ## find potential column name from list of options
  index <- unlist(purrr::map(col_name, 
                             function(x){
                               ## this species will match "metschnikowia_matae_maris.max.cds" and "metschnikowia_matae.max.cds", so force it to match the latter
                               if (x == "metschnikowia_matae")
                               {
                                 x <- paste0(x,".max.cds")
                               }
                               str_which(list_of_option,x)
                             }
  ))
  return(list_of_option[index])
}


get.ws.local <- function(tRNA,      # tRNA gene copy number
                   s_perfect = c(0,0,0,0),  #constraints
                   s_wobble = c(0.41, 0.28, 0.9999, 0.68, 0.89),
                   sking,
                   genetic.code=1)     # super kingdom: 0-eukaryota, 1-prokaryota
{
  s <- c(s_perfect,s_wobble)
  p = 1 - s
  # initialise w vector
  W = NULL  # don't confuse w (lowercase) and W (uppercase)
  
  # obtain absolute adaptiveness values (Ws)
  for (i in seq(1, 61, by=4))
    W = c(W,
          p[1]*tRNA[i]   + p[5]*tRNA[i+1],     # INN -> NNT, NNC, NNA
          p[2]*tRNA[i+1] + p[6]*tRNA[i],       # GNN -> NNT, NNC
          p[3]*tRNA[i+2] + p[7]*tRNA[i],       # TNN -> NNA, NNG
          p[4]*tRNA[i+3] + p[8]*tRNA[i+2])     # CNN -> NNG
  
  # check methionine
  W[36] = p[4]*tRNA[36]
  
  # if bacteria, modify isoleucine ATA codon
  if(sking == 1) W[35] = p[9]
  
  # get rid of stop codons (11, 12, 15) and methionine (36)
  if(genetic.code == 12)
  {
    W = W[-c(11,12,15,20,36)] ## drop S-CTG (20) from consideration
  } else {
    W = W[-c(11,12,15,36)]
  }
  # get ws
  w = W/max(W)
  
  if(sum(w == 0) > 0) {
    ws <- w[w != 0] # zero-less ws
    gm <- exp(sum(log(ws))/length(ws)) # geometric mean
    w[w == 0] = gm # substitute 0-ws by gm
  }
  
  return(w)
}

optimizeS <- function(s_wobble,trna,count.matrix,gene.expression,genetic.code,method.corr)
{
  ws <- get.ws.local(tRNA=trna, s_perfect = c(0,0,0,0), 
                     s_wobble = s_wobble, sking=0, 
                     genetic.code = genetic.code)
  tai <- get.tai(count.matrix, ws)
  corr <- cor(gene.expression,tai,method=method.corr,use="complete.obs")
  -100 * corr ## optim minimizes, so make correlation negative
}

getPhi <- function(species,improved.fits.species)
{
  if (species %in% improved.fits.species){
    phi.df <- read_csv(file.path("/data2/Labella2019/Final_runs/Results_k_2_selectionShared/",species,"run_2/Parameter_est/gene_expression.txt"),col_types = cols()) 
  } else{
    phi.df <- read_csv(file.path("/data2/Labella2019/Final_runs/Results_k_1/",species,"run_2/Parameter_est/gene_expression.txt"),col_types = cols()) 
  }
  phi.df <- phi.df %>%
    dplyr::select(GeneID,Mean)
  return(phi.df)
}


getEmpiricalPhi <- function(target.species,empirical.data,dist.matrix,orthogroups)
{
  dist.from.target <- dist.matrix[,target.species]
  closest.species <- names(dist.from.target[which.min(dist.from.target)])
  ortho.species <- orthogroups %>% 
    dplyr::select(Orthogroups,!!as.name(target.species),!!as.name(closest.species))
  ortho.species <- ortho.species %>% 
    full_join(empirical.data[[closest.species]],
              by=setNames("GeneID",closest.species)
    )
  phi.df <- ortho.species %>% 
    dplyr::select(!!as.name(target.species),Empirical) %>% 
    dplyr::rename(GeneID = !!as.name(target.species),
                  Mean = Empirical) 
  return(phi.df)
}

find.best.weights <- function(species,
                              phi,
                              s_wobble = c(0.41, 0.28, 0.9999, 0.68, 0.89),
                              cutoff.percent = 0.0,
                              total.restarts=1,
                              method="spearman")
{
  ser.genomes <- readLines("~/Labella2019_scripts/ser_fungi.txt")
  codon.order <- read.table("~/dos_reis_tai_codon_order.txt",header=F)[,1]
  
  trna.count <- read_tsv("~/Labella2019_scripts/tRNA/cope_tgcn_per_species.tsv",col_types = cols()) %>%
    filter(Species == species)
         
  trna <- trna.count[,c(codon.order)] %>%
    unlist() %>%
    as.vector()
  genome <- initializeGenomeObject(file.path("/data2/Labella2019/Genomes/cds_cleaned/",species))
  count.matrix <- getCodonCounts(genome)
  gene.names <- rownames(count.matrix)
  count.matrix <- count.matrix[,codon.order]

  if (species %in% ser.genomes)
  {
    count.matrix <- unname(count.matrix[,-c(11,12,15,20,36)]) # drop STOP, M-ATG, S-CTG
  } else {
    count.matrix <- unname(count.matrix[,-c(11,12,15,36)]) # drop STOP, M-ATG
  }
  phi <- phi %>% 
    mutate(Mean = log10(Mean),
           Mean = ifelse(is.finite(Mean),Mean,NA))

  high.expression <- quantile(phi$Mean,probs=cutoff.percent,na.rm=T)
  phi <- phi %>%
    filter(!is.na(Mean) & is.finite(Mean) & !is.na(GeneID)) %>%
    filter(Mean > high.expression)
  count.matrix <- count.matrix[phi$GeneID,]
  gene.expression <- phi$Mean
  
  ws <- get.ws.local(tRNA=trna, s_perfect = c(0,0,0,0), 
                     s_wobble = s_wobble, sking=0, 
                     genetic.code = ifelse(species %in% ser.genomes,12,1))

  tai <- get.tai(count.matrix, ws)
  original.corr <- cor(gene.expression,tai,method=method,use="complete.obs")

  current.best <- data.frame(s5=s_wobble[1],
                             s6=s_wobble[2],
                             s7=s_wobble[3],
                             s8=s_wobble[4],
                             s9=s_wobble[5],
                             Performance.original = original.corr,
                             Performance.optimized = original.corr,
                             Num.Genes=length(gene.expression))
  for (i in 1:total.restarts)
  {
    if (i == 1)
    {
      start <- s_wobble
    } else{
      start <- runif(length(s_wobble)-1,min=0.00001,max=0.99999)
      start <- c(start,s_wobble[5])
    }
    x <- tryCatch({
      best.ws <- optim(par = start,
                  fn=optimizeS,
                  trna=trna,
                  count.matrix=count.matrix,
                  gene.expression=gene.expression,
                  genetic.code = ifelse(species %in% ser.genomes,12,1),
                  method.corr = method,
                  method="L-BFGS-B",
                  lower=rep(0.00001,5),
                  upper=rep(0.99999,5))
      data.frame(s5=best.ws$par[1],
                        s6=best.ws$par[2],
                        s7=best.ws$par[3],
                        s8=best.ws$par[4],
                        s9=best.ws$par[5],
                        Performance.original = original.corr,
                        Performance.optimized = (-1/100)*best.ws$value,
                        Num.Genes=length(gene.expression))
    },
       error=function(e){
         print(str_c(species," encountered an error."))
         data.frame(s5=NA,
                           s6=NA,
                           s7=NA,
                           s8=NA,
                           s9=NA,
                           Performance.original = original.corr,
                           Performance.optimized = NA,
                           Num.Genes=length(gene.expression))
  
       }
    )
    current <- x$Performance.optimized
    previous <- current.best$Performance.optimized
    if (!is.na(current) && current > previous && current > original.corr)
    {
      current.best <- x
    }
  }
  return(current.best)
}


find.best.weights.train.test <- function(species,
                                         phi,
                                         size.train = 0.75,
                                         s_wobble = c(0.41, 0.28, 0.9999, 0.68, 0.89),
                                         cutoff.percent = 0.0,
                                         total.restarts=1,
                                         method = "spearman")
{
  ser.genomes <- readLines("~/Labella2019_scripts/ser_fungi.txt")
  codon.order <- read.table("~/dos_reis_tai_codon_order.txt",header=F)[,1]
  
  trna.count <- read_tsv("~/Labella2019_scripts/tRNA/cope_tgcn_per_species.tsv",col_types = cols()) %>%
    filter(Species == species)
  
  trna <- trna.count[,c(codon.order)] %>%
    unlist() %>%
    as.vector()
  genome <- initializeGenomeObject(file.path("/data2/Labella2019/Genomes/cds_cleaned/",species))
  count.matrix <- getCodonCounts(genome)
  gene.names <- rownames(count.matrix)
  count.matrix <- count.matrix[,codon.order]
  
  if (species %in% ser.genomes)
  {
    count.matrix <- unname(count.matrix[,-c(11,12,15,20,36)]) # drop STOP, M-ATG, S-CTG
  } else {
    count.matrix <- unname(count.matrix[,-c(11,12,15,36)]) # drop STOP, M-ATG
  }
  phi <- phi %>% 
    mutate(Mean = log10(Mean),
           Mean = ifelse(is.finite(Mean),Mean,NA))
  high.expression <- quantile(phi$Mean,probs=cutoff.percent,na.rm=T)
  phi <- phi %>%
    filter(!is.na(Mean) & is.finite(Mean) & !is.na(GeneID)) %>%
    filter(Mean > high.expression)
 
  count.matrix <- count.matrix[phi$GeneID,]
  gene.expression <- phi$Mean
  
  train.index <- sample(1:length(gene.expression),size=floor(length(gene.expression) * size.train),replace = F)
  train.count.matrix <- count.matrix[train.index,]
  train.gene.expression <- gene.expression[train.index]
  
  test.count.matrix <- count.matrix[-train.index,]
  test.gene.expression <- gene.expression[-train.index]
  
  
  ws <- get.ws.local(tRNA=trna, s_perfect = c(0,0,0,0), 
                     s_wobble = s_wobble, sking=0, 
                     genetic.code = ifelse(species %in% ser.genomes,12,1))
  
  tai <- get.tai(train.count.matrix, ws)
  original.corr.train <- cor(train.gene.expression,tai,method=method,use="complete.obs")
  
  tai <- get.tai(test.count.matrix, ws)
  original.corr.test <- cor.test(test.gene.expression,tai,method=method,use="complete.obs")
  
  current.best <- data.frame(s5=s_wobble[1],
                             s6=s_wobble[2],
                             s7=s_wobble[3],
                             s8=s_wobble[4],
                             s9=s_wobble[5],
                             Performance.original = unname(original.corr.test$estimate),
                             Performance.train = unname(original.corr.train),
                             Performance.test = unname(original.corr.test$estimate),
                             # Performance.test.95.lower = unname(original.corr.test$conf.int[1]),
                             # Performance.test.95.upper = unname(original.corr.test$conf.int[2]),
                             Num.Genes=length(gene.expression))
  for (i in 1:total.restarts)
  {
    if (i == 1)
    {
      start <- s_wobble
    } else{
      start <- runif(length(s_wobble)-1,min=0.00001,max=0.99999)
      start <- c(start,s_wobble[5])
    }
    x <- tryCatch({
      best.ws <- optim(par = start,
                       fn=optimizeS,
                       trna=trna,
                       count.matrix=train.count.matrix,
                       gene.expression=train.gene.expression,
                       genetic.code = ifelse(species %in% ser.genomes,12,1),
                       method.corr=method,
                       method="L-BFGS-B",
                       lower=rep(0.00001,5),
                       upper=rep(0.99999,5))
      best.s.wobble <- c(best.ws$par[1],best.ws$par[2],best.ws$par[3],best.ws$par[4],best.ws$par[5])
      test.ws <- get.ws.local(tRNA=trna, s_perfect = c(0,0,0,0), 
                         s_wobble = best.s.wobble, sking=0, 
                         genetic.code = ifelse(species %in% ser.genomes,12,1))

      test.tai <- get.tai(test.count.matrix, test.ws)
      optimized.corr <- cor.test(test.gene.expression,test.tai,method=method,use="complete.obs")
      data.frame(s5=best.ws$par[1],
                 s6=best.ws$par[2],
                 s7=best.ws$par[3],
                 s8=best.ws$par[4],
                 s9=best.ws$par[5],
                 Performance.original = unname(original.corr.test$estimate),
                 Performance.train= (-1/100)*best.ws$value,
                 Performance.test = unname(optimized.corr$estimate),
                 # Performance.test.95.lower = unname(optimized.corr$conf.int[1]),
                 # Performance.test.95.upper = unname(optimized.corr$conf.int[2]),
                 Num.Genes=length(gene.expression))
    },
    error=function(e){
      print(str_c(species," encountered an error."))
      data.frame(s5=NA,
                 s6=NA,
                 s7=NA,
                 s8=NA,
                 s9=NA,
                 Performance.original = unname(original.corr.test$estimate),
                 Performance.train = NA,
                 Performance.test = NA,
                 # Performance.test.95.lower = NA,
                 # Performance.test.95.upper = NA,
                 Num.Genes=length(gene.expression))
      
    }
    )
    current <- x$Performance.test
    previous <- current.best$Performance.test
    if (!is.na(current) && current > previous && current > original.corr.test$estimate)
    {
      current.best <- x
    }
  }
  return(current.best)
}

fungi.tree <- read.tree("~/Labella2019_scripts//tree_with_cds_labels.nwk")
all.species <- readLines("~/Labella2019_scripts/all_fungi.txt")
fungi.tree <- cleanTree(fungi.tree,all.species)
dist.matrix <- cophenetic.phylo(fungi.tree)

species.w.emp <- list.files("/data2/cope/Intragenomic_variation_mutation_bias/Data/Expression/",full.names = F,pattern=".max.cds")
empirical.data.file <- list.files("/data2/cope/Intragenomic_variation_mutation_bias/Data/Expression/",full.names = T,include.dirs = F,pattern=".max.cds")
names(empirical.data.file) <- species.w.emp
dist.matrix <- dist.matrix[species.w.emp,]


improved.fits.species <- readLines("~/Labella2019_scripts/improved_species_by_gene_expression_25_percent_2023_05_03.txt")

species <- list.dirs(file.path("/data2/Labella2019/Final_runs/Results_k_2_selectionShared/"),full.names = F,recursive = F)
names(species) <- species

orthogroups <- read_tsv("~/Labella2019_scripts/shen_etal_ortholog_matrix_2021_12_17.tsv")
orthogroups <- orthogroups %>% 
  column_to_rownames("Orthogroup") %>%
  rename_with(~cleanColumnNames(.x,fungi.tree$tip.label)) %>%
  rownames_to_column("Orthogroups")

empirical.data <- purrr::map(empirical.data.file,function(x) {
  read_csv(x,col_types = cols()) %>% 
    mutate(across(contains("tpm_kallisto"),~na_if(., 0))) %>%
    mutate(across(contains("tpm_kallisto"),log)) %>%
    mutate(
      Log.Mean.Expression=rowMeans(dplyr::select(.,contains("tpm_kallisto")),na.rm=T),
      Empirical=exp(Log.Mean.Expression),
      GeneID=str_remove(GeneID,"-mRNA-1")) %>% 
    dplyr::select(-contains("tpm_kallisto"),-Log.Mean.Expression)
})

cutoff.percent <- 0.0
total.restarts <- 10
phi <- lapply(species,getPhi,improved.fits.species)

best.ws.species <- mclapply(species,
              function(x)
              {
                find.best.weights(x,phi[[x]],cutoff.percent=cutoff.percent,total.restarts=total.restarts)
              },
              mc.cores = 48
              
)
best.ws.species.df <- best.ws.species %>% bind_rows(.id="Species")
write_tsv(best.ws.species.df,"~/Labella2019_scripts/tai_optimization_w_phi_2023_05_03_10_restarts_roc_phi_no_cutoff.tsv")

