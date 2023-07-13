library(AnaCoDa)
library(dplyr)
library(ggplot2)
library(cowplot)

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

# wobble.rules <- function(nucleotide)
# {
#   if (nucleotide == "A")
#   {
#     return(c("T","C","G","A"))
#   } else if (nucleotide == "C"){
#     return(c("C"))
#   } else if(nucleotide == "G"){
#     return(c("G","A"))
#   } else if (nucleotide == "T"){
#     return(c("A","G","C","T"))
#   }
# }

wobble.rules <- function(nucleotide)
{
  if (nucleotide == "A")
  {
    return(c("T","G","A"))
  } else if (nucleotide == "C"){
    return(c("C"))
  } else if(nucleotide == "G"){
    return(c("G","A"))
  } else if (nucleotide == "T"){
    return(c("C","T"))
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

calculateCorrelations <- function(aa.include=NULL)
{
  
  deta <- read.table("Post_analysis/deta_matrix_2022.tsv",sep="\t",header=T,stringsAsFactors = F,row.names=1)
  files <- list.files(path="tRNA/tGCN/",pattern=".tsv",full.names = F,recursive = F)
  aa <- aminoAcids()
  
  cor.df <- data.frame(Species=colnames(deta),Correlation=numeric(length(colnames(deta))),Pvalue=numeric(length(colnames(deta))),stringsAsFactors = F)
 
  ctg.ser <- F
  for (i in 1:length(files))
  {
    f <- files[i]
    species <- unlist(strsplit(f,".tsv"))[1]
    cognate <- read.table(file.path("tRNA","Cognate",f),sep="\t",header=T,stringsAsFactors = F,row.names = 1)
    deta.spec <- deta[,species,drop=F]
    
    deta.spec <- deta.spec[which(!is.na(deta.spec)),,drop=F]
    if (!"CTG" %in% rownames(deta.spec))
    {
      ctg.ser <- T
      cognate["CTG","AA"] <- "J"
    } else{
      ctg.ser <- F
    }
    deta.trna <- merge(deta.spec,cognate,by=0,all.y = T)
    deta.trna[is.na(deta.trna)] <- 0
    deta.trna <- deta.trna[which(!deta.trna$AA %in% c("M","X","W","J")),]
    for (a in aa)
    {
      if (a == "M" || a == "X" || a == "W") next
     if(a == "L" && ctg.ser)
     {
        codons.w.ref <- c("CTA","CTC","CTT","TTA","TTG")
        codons.wo.ref <- c("CTA","CTC","CTT","TTA")
     } else {
        codons.w.ref <- AAToCodon(a,F)
        codons.wo.ref <- AAToCodon(a,T)
      }
      ref.codon <- codons.w.ref[which(!codons.w.ref %in% codons.wo.ref)] 
      ref.row <- which(deta.trna$Row.names == ref.codon)
      codon.rows <- which(deta.trna$Row.names %in% codons.w.ref)
      deta.trna[codon.rows,"1/w.tRNA"] <- 1/deta.trna[codon.rows,"w.tRNA"] - 1/deta.trna[ref.row,"w.tRNA"]
      write.table(deta.trna,file.path("tRNA","Cognate_w_trna",f),col.names=T,row.names=F,quote=F,sep="\t")
    }
    if (!is.null(aa.include))
    {
      keep <- which(deta.trna$AA %in% aa.include)
    } else{
      keep <- which(deta.trna$AA %in% deta.trna$AA)
    }
    x <- cor.test(deta.trna[keep,"1/w.tRNA"],deta.trna[keep,2],use = "pairwise.complete.obs",method="spearman")
    cor.df[i,"Species"] <- species
    cor.df[i,"Correlation"] <- unname(x$estimate)
    cor.df[i,"Pvalue"] <- unname(x$p.value)
    if (x$estimate < 0)
    {
      print(f)
    }
    
  }
  return(cor.df)
}

plotCorrelationDist <- function(cor.df,title="Selection coefficient vs. Relative Waiting Time Estimate\nbased on tRNA Abundances")
{
  cor.mean <- aggregate(cor.df$Correlation,list(cor.df$CTG),median)
  p <- (ggplot(cor.df) + geom_histogram(binwidth=0.1,alpha=0.5,aes(x=Correlation,fill=CTG)) 
        + theme_cowplot()
        + xlab("Spearman Correlation")
        + ylab("Count")
        + geom_vline(aes(xintercept=cor.mean[which(cor.mean$Group.1 == "Leucine"),"x"],color="Leucine"), linetype = "longdash")
        + geom_vline(aes(xintercept=cor.mean[which(cor.mean$Group.1 == "Serine"),"x"],color="Serine"), linetype = "longdash")
        + scale_color_manual(name="Median Spearman Correlation",values=c(Leucine="red",Serine="blue"))
        + ggtitle(title)
  )
  return(p)
}

createtRNAFiles()
aa.groups <- list(c("C","D","E","F","H","I","K","N","Q","Y","Z","A","G","P","S","T","V","R","L")
                  ,c("C","D","E","F","H","I","K","N","Q","Y","Z"),
                                     c("A","G","P","S","T","V"),
                                     c("R","L"))
ser <- read.table("ser_fungi.txt",sep="",header=F,stringsAsFactors = F)

base.title <- "Selection coefficient vs. Relative Waiting Time Estimate\nbased on tRNA Abundances"

hist.plots <- vector(mode = "list",length=length(aa.groups))

for (i in 1:length(aa.groups))
{
  if (i == 1)
  {
    suffix.title <- "All Amino Acids"
  } else if (i == 2) {
    suffix.title <- "2/3 Codon Amino Acids"
  } else if (i == 3) {
    suffix.title <- "4 Codon Amino Acids"
  } else if (i == 4) {
    suffix.title <- "5/6 Codon Amino Acids"
  }
  cor.df <- calculateCorrelations(aa.groups[[i]])
  cor.df[,"CTG"] <- ifelse(cor.df$Species %in% ser[,1],"Serine","Leucine")
  p <- plotCorrelationDist(cor.df,suffix.title)
  hist.plots[[i]] <- p
}

title <- title <- ggdraw() + draw_label(
  base.title,
  fontface = 'bold',
  x = 0,
  hjust =0,size=24
) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 2)
  )
top <- plot_grid(hist.plots[[1]], hist.plots[[2]], labels = c('A', 'B'), label_size = 16)
bottom <- plot_grid(hist.plots[[3]], hist.plots[[4]], labels = c('C', 'D'), label_size = 16)

total.plot <- plot_grid(title,top,bottom,ncol=1,rel_heights=c(0.1,1,1))

ggsave2(filename="correlation_selection_coefficients_wait_times_by_trna_ser_reruns.png",plot=total.plot,dpi=600,width=16,height=16)

# files <- list.files(path="tRNA/Cognate_w_trna/",pattern=".tsv",full.names = F,recursive = F)
# num_files <- length(files)
# diff.species.corr <- c()
# 
# bad.ser <- read.table("bad_ctg_ser.txt",sep="",header=F,stringsAsFactors = F)
# 
# for (i in 1:(num_files-1))
# {
#   for (j in (i+1):num_files)
#   {
#     col.name.1 <- unlist(strsplit(files[i],split = ".tsv",fixed = T))[1]
#     col.name.2 <- unlist(strsplit(files[j],split = ".tsv",fixed = T))[1]
#     #if (!col.name.1 %in% ser[,1] || !col.name.2 %in% ser[,1] || !col.name.1 %in% bad.ser[,1] || !col.name.2 %in% bad.ser[,1]) next
#     df.1 <- read.table(file.path("tRNA","Cognate_w_trna",files[i]),sep="\t",header=T,stringsAsFactors = F)
#     df.2 <- read.table(file.path("tRNA","Cognate_w_trna",files[j]),sep="\t",header=T,stringsAsFactors = F)
#     df.comb <- merge(df.1,df.2,by="Row.names",suffixes=c("_X","_Y"))
# 
#     deta.diff <- df.comb[,col.name.1] - df.comb[,col.name.2]
#     trna.diff <- df.comb[,"X1.w.tRNA_X"] - df.comb[,"X1.w.tRNA_Y"]
#     y <- cor(trna.diff,deta.diff,use = "pairwise.complete.obs",method="spearman")
#     diff.species.corr <- c(diff.species.corr,y)
#   }
# }
# 
# pdf("correlation_between_changes_in_trna_and_deta_between_species_ser_reruns.pdf")
# hist(diff.species.corr,nclass=40,main="Correlation between changes in selection coefficients and\n tRNA gene copy number between species",xlab="Spearman Rank Correlation")
# dev.off()
