library(stringr)

convertoCdsIDs <- function(prot.ids)
{
  if (!is.na(prot.ids) && !is.null(prot.ids))
  {
    present <- gregexpr("gene=\\S+",prot.ids)
    if (present[[1]][1] == -1)
    {
      prot.ids.split <- unlist(strsplit(prot.ids,","))
      y <- unlist(lapply(strsplit(unlist(lapply(prot.ids.split,trimws))," "),function(i){i[1]}))
    }else{
      x<-unlist(strsplit(unlist(regmatches(prot.ids,present)),"gene=",fixed=T))
      y <- x[which(x != "")]
    }
    if (any(is.na(y))) {print(prot.ids);print(y);stop()}
    if (any(y=="NA")) {print(prot.ids);print(y);stop()}
    return(paste(y,collapse=","))
  }else{
    return("")
  }
}

geom.mean <- function(x)
{
  num.na <- length(which(is.na(x)))
  l <- length(x) - num.na
  return(prod(x)^(1/l))
}

getPhi <- function(spec.prot,phi)
{
  phi.reorder <- numeric(length(spec.prot))
  for (i in 1:length(spec.prot))
  {

    if (is.na(spec.prot[i]) || spec.prot[i]=="")
    {
      phi.reorder[i] <- NA
    } else{
      proteins <- unlist(strsplit(spec.prot[i],","))
      
      phi.prots <- phi[proteins,"Mean",drop=F]
      if (nrow(phi.prots) > 1)
      {
        phi.reorder[i] <- geom.mean(phi.prots$Mean)
        #if (is.na(phi.reorder[i])){print(proteins);print(phi.prots);print(proteins[which(!proteins %in% rownames(phi.prots))]);stop()}
      } else{
        phi.reorder[i] <- phi.prots$Mean
      }
    }
  }
  return(phi.reorder)
}


updateOrthoFinderWCDS <- function()
{
  ortho <- read.table("/data/cope/Labella2019/Genomes/pep/OrthoFinder/Results_Jul29_2/Orthogroups/Orthogroups.tsv",sep="\t",header=T,stringsAsFactors = F,row.names = 1,na.strings=c("NA",""))
  leu.genomes <- read.table("all_fungi.txt",header=F,stringsAsFactors = F)
  leu.genomes <- unlist(strsplit(leu.genomes[,1],".cds"))
  ortho <- ortho[,leu.genomes]
  index <- rownames(ortho)
  ortho.prot <- matrix(rep(NA,nrow(ortho)*ncol(ortho)),nrow = nrow(ortho),ncol=ncol(ortho))
  ortho.prot <- as.data.frame(ortho.prot,stringsAsFactors=F)
  rownames(ortho.prot) <- rownames(ortho)
  colnames(ortho.prot) <- colnames(ortho)
  for (i in leu.genomes)
  {
    print(i)
    tmp <- ortho[,i,drop=F]
    z <- apply(tmp,1,convertoCdsIDs)
    ortho.prot[,i] <- z[index]
  }
  write.table(ortho.prot,"ortholog_matrix_all_fungi_fixed.tsv",sep="\t",row.names=T,col.names=T,quote=F)
  return(ortho.prot)
}


#ortho.prot <- updateOrthoFinderWCDS()
prot <- read.table("ortholog_matrix_all_fungi.tsv",sep="\t",header=T,row.names=1,stringsAsFactors = F)
all.genomes <- read.table("all_fungi.txt",header=F,stringsAsFactors = F)
bad.fits <- read.table("improved_species_2021_11_24.txt",header=F,stringsAsFactors=F)
all.genomes.no.cds <- unlist(strsplit(all.genomes[,1],".cds"))
bad.fits.no.cds <- unlist(strsplit(bad.fits[,1],".cds"))
prot <- prot[,all.genomes.no.cds]

# prot[] <- lapply(prot,
#                  function(x){
#                    paralogs <- str_detect(x,",");
#                    x[which(paralogs)] <- NA
#                    return(x)
#                  })
head.directory <- "/data2/Labella2019/Results/"

prot.phi <- prot

for (i in all.genomes.no.cds)
{
  if (i %in% bad.fits.no.cds)
  {
    phi.file <- file.path("/data2/Labella2019/Final_runs/Results_k_2_selectionShared/",paste0(i,".cds"),"run_2","Parameter_est","gene_expression.txt")
  } else
  {
    phi.file <- file.path("/data2/Labella2019/Final_runs/Results_k_1/",paste0(i,".cds"),"run_2","Parameter_est","gene_expression.txt")
  } 
  print(i)
  phi <- read.table(phi.file,sep=",",header=T,stringsAsFactors = F,row.names=1)
  rownames(phi) <- unlist(strsplit(rownames(phi),"-mRNA-1(_1)?"))
  #prot.phi[,i] <- phi[prot[,i],1]
  prot.phi[,i] <- getPhi(prot[,i],phi)
}
colnames(prot.phi) <- paste0(colnames(prot.phi),".cds")

for (i in colnames(prot.phi))
{
  not.na <- which(!is.na(prot.phi[,i]))
  if (length(not.na) == 0)
  {
    print(i)
  }
}

write.table(prot.phi,"phi_matrix_all_fungi_geom_mean_for_paralogs_2021_11_24.tsv",sep="\t",row.names=T,col.names=T,quote=F)
