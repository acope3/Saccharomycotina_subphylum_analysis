---
title: "test"
author: "Alex Cope"
date: "8/12/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(ggside)
library(cowplot)
library(tidyverse)
#library(rtracklayer)
```

```{r fig.width=12,fig.height=12}
gffToDataFrame <- function(gff,gff.id)
{
  chr <- as.character(gff$seqid)
  tid <- gff[,gff.id]
  df <- data.frame(Chr=chr,Gene=tid,stringsAsFactors = F)
  df <- df %>% distinct()
  return(df)
}


compareGC3ToClusterFreq <- function(cluster.file,gc.file,gff.file,gff.id="transcript_id",chr.reg=NULL)
{
  gff <- rtracklayer::readGFF(gff.file)
  gff.cds <- gff[which(gff$type == "CDS"),]

 
  clust <- read_tsv(cluster.file,col_types=cols())
  gc <- read_tsv(gc.file,col_types=cols())
  
 
  gff.cds$ID <- str_remove(gff.cds[,gff.id],"-[A-Za-z0-1]+")
  gff.cds <- gff.cds[which(gff.cds$ID %in% gc$Gene),]
  #gff.cds$ID <- str_remove(gff.cds[,gff.id],"_[AB]")
  
  clust$Gene <- str_remove(clust$Gene,"(_Allele){0,1}")
  gc$Gene <- str_remove(gc$Gene,"(_Allele){0,1}")
 
  gc.clust <- gc %>% left_join(clust,by="Gene")
  max.cluster <- gc.clust  %>% group_by(Cluster) %>% 
               summarize(Mean.GC3 = mean(GC3)) %>% 
               top_n(1) %>%
               dplyr::rename(Max.Cluster=Cluster)

  gc.clust <- gc.clust %>% mutate(Max.Cluster = max.cluster$Max.Cluster)


  gc.clust <- gc.clust %>% 
              mutate(Cluster = case_when(Cluster == 1 & Max.Cluster == 1 ~ 2,
                                         Cluster == 2 & Max.Cluster == 1 ~ 1,
                                         Cluster == 1 & Max.Cluster == 2 ~ 1,
                                         Cluster == 2 & Max.Cluster == 2 ~ 2)
              )

  gene.names <- gc.clust$Gene
  
  gc.clust <- gc.clust %>% column_to_rownames("Gene")
  gc.clust <- gc.clust %>% add_column(Gene = gene.names)
  
 
  gc.clust <- gc.clust[unique(gff.cds$ID),]
  gc.clust <- gc.clust[which(!is.na(gc.clust$Gene)),]
  gff.cds <- gffToDataFrame(gff.cds,gff.id)
  if(!is.null(chr.reg))
  {
    gff.cds$Chr <- str_extract(gff.cds$Chr,pattern = chr.reg)
  }
  gc.clust <- gc.clust %>% left_join(gff.cds,by="Gene")

  gc.clust$Cluster_binary <- ifelse(gc.clust$Cluster == 1,0,1)
  gc.clust <-gc.clust %>% distinct()
  gc.clust <- gc.clust %>% group_by(Chr) %>% mutate(Gene.Number=seq(1,length(Chr)),
                                                                        `GC3%` = data.table::frollmean(GC3,n=30,fill=NA,na.rm=T,align="center"),
                                                                        `Higher GC3% Cluster Frequency`=data.table::frollmean(Cluster_binary,n=30,fill=NA,na.rm=T,align="center")) %>%
                ungroup()
gc.clust <- gc.clust %>% mutate(Cluster = case_when(Cluster == 1 ~ "Lower GC3%",
                                         Cluster == 2 ~ "Higher GC3%",
                                         is.na(Cluster) ~ "Excluded CDS")
                    )

  return(gc.clust)
}




gc.clust <- compareGC3ToClusterFreq("Clara_k_2/candida_albicans.max.cds","GC/candida_albicans.max.cds","Genomes/Genomes/gff/candida.albicans.max.cds",gff.id="ID",chr.reg="chr[0-9RM]")
gc.clust <- gc.clust %>% filter(Chr != "chrM")
gc.clust.long <- gc.clust %>%
                    pivot_longer(cols=c(`GC3%`,`Higher GC3% Cluster Frequency`),
                                 names_to="Moving.Average",
                                 values_to="Frequency")

color.scheme <- c("#000000","#E41A1C","#377EB8")
names(color.scheme) <- c("Excluded CDS","Lower GC3%","Higher GC3%")

gc3.lineplot <- ggplot() + 
           geom_xsidetile(data=gc.clust.long,aes(x=Gene.Number,y=1,height=0.5,xfill=Cluster)) +
           ggside(x.pos="bottom") +
           geom_line(data=gc.clust.long,mapping=aes(x=Gene.Number,y=Frequency,linetype=Moving.Average),show.legend = T) +
           theme_cowplot() +
           theme(aspect.ratio = 1) +
           ylim(c(0,1)) +
           scale_xfill_manual(values=color.scheme,name="Cluster") +
           scale_linetype_manual(values=c(`GC3%`=1,`Higher GC3% Cluster Frequency`=2),name = "Moving Average")+
           ggtitle("Variation in GC3% and Cluster Frequency in Chromosome 6 of C. albicans") +
           xlab("Gene Number") +
           ylab("Frequency") +
           scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
           scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
           facet_wrap(~Chr,scales="free",nrow=2,ncol=4)
gc3.lineplot

```
