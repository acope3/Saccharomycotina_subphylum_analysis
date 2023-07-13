library(AnaCoDa)


complement <- function(nuc)
{
  if (nuc == "A")
  {
    return("T")
  } else if (nuc == "C")
  {
    return("G")
  } else if (nuc == "G" )
  {
    return("C")
  } else if (nuc == "T")
  {
    return("A")
  }
}


reverseComplement <- function(anticodon)
{
  nucleo <- unlist(strsplit(anticodon,""))
  nucleo <- rev(nucleo)
  codon <- unlist(lapply(nucleo,complement))
  codon <- paste0(codon,collapse="")
  return(codon)
}



three.to.one.code <- function(three.key,codon)
{
  if(three.key == "Ala")
  {
    return("A")
  } else if (three.key == "Cys") {
    return("C")
  } else if (three.key == "Asp"){
    return("D")
  } else if(three.key == "Glu"){
    return("E")
  } else if (three.key == "Phe"){
    return("F")
  } else if(three.key == "Gly"){
    return("G")
  } else if(three.key == "His"){
    return("H")
  } else if(three.key == "Ile"){
    return("I")
  } else if(three.key == "Lys"){
    return("K")
  } else if(three.key == "Leu"){
    return("L")
  } else if(three.key == "Met"){
    return("M")
  } else if(three.key == "Asn"){
    return("N")
  } else if(three.key == "Pro"){
    return("P")
  } else if(three.key == "Gln"){
    return("Q")
  } else if(three.key == "Arg"){
    return("R")
  } else if(three.key == "Ser"){
    
    if(codon == "AGT" || codon == "AGC" || codon=="CTG")
    {
 
      return("Z")
    }else{
     return("S")
    }
  } else if(three.key == "Thr"){
    return("T")
  }else if(three.key == "Val"){
    return("V")
  }else if(three.key == "Trp"){
    return("W")
  }else if(three.key == "Tyr"){
    return("Y")
  }else if(three.key == "Stop"){
    return("X")
  }
}
converttRNAColnames <- function(trna)
{
  anti <- colnames(trna)
  y<-strsplit(x=anti,split=".",fixed=T)
  codon.aa <- lapply(y,function(x){paste(reverseComplement(x[1]),x[2],sep=".")})
  colnames(trna) <- codon.aa
  return(trna)
}

to.remove <- c("CAG.Ala","TTA.Stop","TCA.Stop","CTA.Stop")
to.include <- read.table("all_fungi.txt",sep="",header=F,stringsAsFactors=F)
#fungi.tree <- read.tree("tree_with_cds_labels.nwk")
trna <- read.table("Labella2019_tRNA_per_genome.tsv",sep="\t",header=T,stringsAsFactors = F,row.names=1)
trna <- trna[to.include[,1],]
trna <- trna[,which(!colnames(trna) %in% to.remove)]
trna <- trna[,4:ncol(trna)]
trna <- converttRNAColnames(trna)
header <- colnames(trna)

trna.aa <- strsplit(header,split = ".",fixed=T)

df <- data.frame("AA"=character(length=61),"Codon"=character(length=61),"tRNA"=numeric(length=61),stringsAsFactors = F)

for (i in 1:length(trna.aa))
{
  df[i,"AA"] <- three.to.one.code(trna.aa[[i]][2],trna.aa[[i]][1])
  df[i,"Codon"] <- trna.aa[[i]][1]
}


for (i in 1:nrow(trna))
{
  row <- trna[i,]
  filename <- rownames(trna)[i]
  df[,"tRNA"] <- unname(unlist(row))
  write.table(df,paste0("tRNA/tGCN/",filename,".tsv"),col.names = T,row.names = F,quote=F,sep="\t")
}