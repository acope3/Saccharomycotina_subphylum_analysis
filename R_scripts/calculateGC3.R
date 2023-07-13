library(seqinr)
library(AnaCoDa)
library(parallel)

createGC3Table <- function(cds.file)
{
  species <- unlist(strsplit(cds.file,"/"))
  species <- species[length(species)]
  genome <- initializeGenomeObject(cds.file)
  size <- length(genome)
  genes <- lapply(1:size,function(x){tmp<-genome$getGeneByIndex(x,F);tmp$seq})
  genes <- lapply(genes,function(x){unlist(strsplit(x,split="",fixed=T))})
  ids <- lapply(1:size,function(x){tmp<-genome$getGeneByIndex(x,F);tmp$id})
  gc <- lapply(genes,GC)
  gc3 <- lapply(genes,GCpos,pos=3)
  gc3.df <- data.frame(Gene=unlist(ids),GC=unlist(gc),GC3=unlist(gc3),stringsAsFactors = F)
  write.table(gc3.df,file.path("GC",species),sep="\t",col.names=T,row.names=F,quote=F)

}

cds.files <- list.files("Genomes/Genomes/cds/",pattern = ".max.cds",full.names = T)
mclapply(cds.files,createGC3Table,mc.cores=12)


