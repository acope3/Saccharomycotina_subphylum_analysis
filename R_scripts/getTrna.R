canonical.wobble <- 0.61 ## purine-purine or pyrmidine-pyrimidine
non.canonical.wobble <- 0.64 ## purine-pyrimidine (GT/AC)

reverse.complement <- function(seq)
{
  seq.split <- unlist(strsplit(seq,split=""))
  reverse <- rev(seq.split)
  trna <- unlist(lapply(reverse,complement.rules))
  trna <- paste(trna,collapse='')
  return(trna)
  
}


complement.rules <- function(nucleotide)
{
  if (nucleotide == "A")
  {
    return("T")
  } else if (nucleotide == "C"){
    return("G")
  } else if(nucleotide == "G"){
    return("C")
  } else if (nucleotide == "T"){
    return("A")
  }
}

wobble.rules <- function(nucleotide)
{
  if (nucleotide == "A")
  {
    return(c("T","C","G","A"))
  } else if (nucleotide == "C"){
    return(c("C"))
  } else if(nucleotide == "G"){
    return(c("G","A"))
  } else if (nucleotide == "T"){
    return(c("A","G","C","T"))
  }
}

## assumes seq is the trna sequences ordered from 5'-end to 3'-end (will match codon from 3'-end to 5'-end)
wobble.neighbors <- function(seq)
{
  seq.split <- unlist(strsplit(seq,split=""))
  wobbles <- wobble.rules(seq.split[1])
  second <- rep(seq.split[2],length(wobbles))
  third <- rep(seq.split[3],length(wobbles))
  wobble.trna <- paste0(wobbles,second,third)
  return(wobble.trna)
}

get.neighbors <- function(codon)
{
  neighbors <- character(length=12)
  nuc <- unlist(strsplit(codon,split='',fixed=T))
  neighbors[1:4] <- paste0(nuc[1],c("A","G","T","C"),nuc[3])
  neighbors[5:8] <- paste0(c("A","G","T","C"),nuc[2],nuc[3])
  neighbors[9:12] <- paste0(nuc[1],nuc[2],c("A","G","T","C"))
  return(neighbors)
}

checkWobble <- function(codon,trna)
{
  reverse <- reverse.complement(codon)
  trna.split <- unlist(strsplit(trna,""))
  codon.split <- unlist(strsplit(codon,""))
  reverse.split <- unlist(strsplit(reverse,""))
  if(reverse.complement(codon) == trna)
  {
    return(1)
  } else{
    index <- which(reverse.split != trna.split)
    if(index == 1)
    {
      t.index <- 1
      c.index <- 3
    } else if (index == 2){
      t.index <- 2
      c.index <- 2
    } else if (index == 3){
      t.index <- 3
      c.index <- 1
    } else{
      print("Error")
    }
    if (trna.split[t.index] == "C" && codon.split[c.index] == "T")
    {
      return(canonical.wobble)
    } else if (trna.split[t.index] == "T" && codon.split[c.index] == "C"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "G" && codon.split[c.index] == "A"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "A" && codon.split[c.index] == "G"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "G" && codon.split[c.index] == "T"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "T" && codon.split[c.index] == "G"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "A" && codon.split[c.index] == "C"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "C" && codon.split[c.index] == "A"){
      return(non.canonical.wobble)
    } else if (trna.split[t.index] == "A" && codon.split[c.index] == "A"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "C" && codon.split[c.index] == "C"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "G" && codon.split[c.index] == "G"){
      return(canonical.wobble)
    } else if (trna.split[t.index] == "T" && codon.split[c.index] == "T"){
      return(canonical.wobble)
    }
  }
}
## Given list of tRNA, identify cognate, pseudo-cognate, and near-cognate tRNA for all codons
getCognate <- function(target.file,output)
{
  trna <- read.table(target.file,sep="\t",header=T,stringsAsFactors = F)
  trna[,"AntiCodon"] <- unlist(lapply(trna$Codon,reverse.complement))
  trna <- trna[which(is.na(trna$tRNA) == F),]
  trna.present <- trna[which(trna$tRNA > 0),]
  anticodons <- trna.present$AntiCodon
  df <- data.frame(Codon = trna$Codon,AA=trna$AA,Cognates=character(length(nrow(trna))),Pseudo=character(length(nrow(trna))),Near=character(length(nrow(trna))),stringsAsFactors=F)
  for (i in 1:length(anticodons))
  {
    aa <- trna.present[i,"AA"]
    wobbles <- unlist(lapply(wobble.neighbors(anticodons[i]),reverse.complement))
    for (j in wobbles)
    {
      if(j %in% c("TAA","TAG","TGA")) next
      if (trna[which(trna$Codon== j),"AA"] == aa)
      {
        df[which(df$Codon == j),"Cognates"] <- paste(unlist(df[which(df$Codon == j),"Cognates"]),anticodons[i],sep=",",collapse = "")
        
      } else{
        df[which(df$Codon == j),"Near"] <- paste(unlist(df[which(df$Codon == j),"Near"]),anticodons[i],sep=",",collapse = "")
      }
    }
    
  }
  
  codons <- trna$Codon
  for (i in 1:length(codons))
  {
    aa <- trna[i,"AA"]
    trna.neighbors <- get.neighbors(codons[i])
    
    for (j in trna.neighbors)
    {
      if(j %in% c("TAA","TAG","TGA")) next
      cognates <- unlist(strsplit(df[i,"Cognates"],","))
      reverse <- reverse.complement(j)
      if (j %in% trna.present$Codon && (!reverse %in% cognates))
      {
        
        if (trna[which(trna$Codon== j),"AA"] == aa)
        {
          tmp <- df[i,"Pseudo"]
          tmp.split <- unlist(strsplit(tmp,","))
          
          if (!reverse %in% tmp.split)
          {
            df[i,"Pseudo"] <- paste(df[i,"Pseudo"],reverse,sep=",",collapse = "")
          }
          
        } else {
          tmp <- df[i,"Near"]
          tmp.split <- unlist(strsplit(tmp,","))
          if (!reverse %in% tmp.split)
          {
            df[i,"Near"] <- paste(df[i,"Near"],reverse,sep=",",collapse = "")
          }
        }
      }
    }
  }
  for (i in 1:nrow(df))
  {
    
    tmp <- unlist(strsplit(df[i,"Cognates"],split = ",",fixed=T))
    if (length(tmp) >= 2)
    {
      df[i,"Cognates"] <- paste(tmp[2:length(tmp)],sep=",",collapse=",")
    }
    tmp <- unlist(strsplit(df[i,"Pseudo"],split = ",",fixed=T))
    
    if (length(tmp) >= 2)
    {
      df[i,"Pseudo"] <- paste(tmp[2:length(tmp)],sep=",",collapse=",")
    }
    tmp <- unlist(strsplit(df[i,"Near"],split = ",",fixed=T))
    
    if (length(tmp) >= 2)
    {
      df[i,"Near"] <- paste(tmp[2:length(tmp)],sep=",",collapse=",")
    }
  }
  write.table(df,output,sep="\t",row.names=F,col.names=T,quote=F)
  
}

updatetRNATable <- function(trna.file,cognate.file,include.pseudo=F)
{
  trna <- read.table(trna.file,sep="\t",header=T,stringsAsFactors = F)
  cognate <- read.table(cognate.file,sep="\t",header=T,stringsAsFactors = F)
  ctg <- trna[trna$Codon == "CTG","tRNA"]
  
  trna <- trna[which(is.na(trna$tRNA) == F),]
  rownames(trna) <- trna[,2]
  for (i in 1:nrow(cognate))
  {
    row <- cognate[i,]
    codon <- row$Codon
    cognates <- unlist(strsplit(row$Cognate,","))
    pseudo <- unlist(strsplit(row$Pseudo,","))
    near <- unlist(strsplit(row$Near,","))
    r_c <- 0
    for (j in cognates)
    {
      r_c <- r_c + (trna[reverse.complement(j),"tRNA"] *checkWobble(codon,j))
    }
    if (include.pseudo)
    {
      for (j in pseudo)
      {
        r_c <- r_c + (trna[reverse.complement(j),"tRNA"] * checkWobble(codon,j))
      }
    }
    cognate[i,"w.tRNA"] <- r_c
  }
  write.table(cognate,cognate.file,sep="\t",row.names=F,col.names = T,quote = F)
  
}


createtRNAFiles <- function()
{
  files <- list.files(path="tRNA/tGCN/",pattern=".tsv",full.names = F,recursive = F)
  for (i in 1:length(files))
  {
    f <- files[i]
    species <- unlist(strsplit(f,".tsv"))[1]
    getCognate(file.path("tRNA","tGCN",f),output=file.path("tRNA","Cognate",f))
    updatetRNATable(file.path("tRNA","tGCN",f),file.path("tRNA","Cognate",f))
  }
}



createtRNAMAtrix <- function(results.dir)
{
  runs <- list.files(results.dir,full.names = F,recursive = F)
  runs <- gsub("\\.tsv","",runs)
  csp.matrix <- matrix(rep(0,40*length(runs)),nrow=40,ncol = length(runs))
  
  csp.matrix <- as.data.frame(csp.matrix)
  colnames(csp.matrix) <- runs
  aa <- AnaCoDa::aminoAcids()
  codons <- c()
  for (a in aa)
  {
    aa.codons <- AnaCoDa::AAToCodon(a,T)
    codons <- c(codons,aa.codons)
  }
  rownames(csp.matrix) <- codons
  for (i in runs)
  {
    trna.file <- file.path(results.dir,paste0(i,".tsv"))
    trna <- read.table(trna.file,sep="\t",header=T,stringsAsFactors = F,row.names = 1)
    csp.matrix[,i] <-trna[rownames(csp.matrix),"relative.w.tRNA"]
    
  }
  return(csp.matrix)
}



trna.rel.rates <- createtRNAMAtrix("tRNA/Cognate_w_trna/")
write.table(trna.rel.rates,"trna_relative_rates.tsv",col.names=T,row.names=T,quote=F)
