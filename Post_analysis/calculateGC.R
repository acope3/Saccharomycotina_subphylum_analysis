library(seqinr)
library(stringi)
library(ape)

plotComp <- function(x,y,...)
{
  xlim <- range(x)
  ylim <- range(y)
  corr <- round(cor(x,y),3)
  height <- ylim[2] - ylim[1]
  width <- xlim[2] - xlim[1]
  plot(x,y,...)
  if (corr < 0)
  {
    text(x=xlim[1]+0.1*width,y=ylim[1]+0.1*height,labels=bquote(rho~"="~.(corr)))
    print(xlim[1]+0.1*width)
    print(ylim[1]+0.1*height)
  } else{
    text(x=xlim[2]-0.1*width,y=ylim[1]+0.1*height,labels=bquote(rho~"="~.(corr)))
  }
}



all.genomes <- read.table("../all_fungi.txt",header=F,stringsAsFactors = F)

cds.dir <-"/data/cope/Labella2019/Genomes/cds_cleaned/"

cds <- list.files(cds.dir,pattern = "*.cds",full.names = T)

cds.genome <- list.files(cds.dir,pattern = "*.cds",full.names = F)

seq.fasta <- lapply(cds, read.fasta,set.attributes=F,seqonly=T,as.string=F,forceDNAtolower = TRUE)

combined.cds <- lapply(seq.fasta,paste0,collapse="")

combined.cds.split <- lapply(combined.cds,strsplit,split="")


gc.content <- lapply(combined.cds.split,function(x){GC(unlist(x))})

gc.content.unlist <- unlist(gc.content)

df <- data.frame(Species=cds.genome,GC=gc.content.unlist,stringsAsFactors = F)

pep.dir <-"/data/cope/Labella2019/Genomes/pep/"

pep <- list.files(pep.dir,pattern = "*.pep",full.names = T)

pep.genome <- list.files(pep.dir,pattern = "*.pep",full.names = F)

seq.fasta <- lapply(pep, read.fasta,set.attributes=F,seqtype = "AA",seqonly=T,as.string=F)

combined.pep <- lapply(seq.fasta,paste0,collapse="")

lys.count <- stri_count_regex(combined.pep,"K")
arg.count <- stri_count_regex(combined.pep,"R")

total.pos <- lys.count + arg.count
freq.lys <- lys.count/arg.count

df[,"Lysine:Arginine"] <- freq.lys

rownames(df) <- df$Species

tree <- read.tree("../tree_with_cds_labels.nwk")
df <- df[tree$tip.label,]

gc.pic <- pic(df$GC,tree)
lys.pic <- pic(df$`Lysine:Arginine`,tree)

pdf("comparing_gc_content_to_lys_arg_ratio.pdf")
plotComp(df$GC,df$`Lysine:Arginine`,xlab="CDS GC%",ylab="Lysine:Arginine",main="Comparing GC% to K:R")
plotComp(gc.pic,lys.pic,xlab="CDS GC% (PIC)",ylab="Lysine:Arginine (PIC)",main="Comparing GC% to K:R\nPhylogenetic Independent Contrasts")

deta <- read.table("leu_dEta.tsv",sep="\t",header=T,stringsAsFactors = F,row.names=1)
deta <- deta[c("AAA","AGA","AGG","CGA","CGC","CGG"),]
deta.t <- t(deta)

df.deta <- merge(df,deta.t,by=0)
rownames(df.deta) <- df.deta$Row.names
leu.genomes <- read.table("../leu_fungi.txt",header=F,stringsAsFactors = F)
leu.tree <- drop.tip(tree,tree$tip.label[which(!tree$tip.label %in% leu.genomes[,1])])
df.deta<-df.deta[leu.tree$tip.label,]
for (i in 5:ncol(df.deta))
{
  gc.pic <- pic(df.deta$GC,leu.tree)
  sel.pic <- pic(df.deta[,i],leu.tree)
  plotComp(df.deta$GC,df.deta[,i],xlab="GC%",ylab="Selection Coefficient",main=colnames(df.deta)[i])
  plotComp(gc.pic,sel.pic,xlab="GC% (PIC)",ylab="Selection Coefficient (PIC)",main=paste0(colnames(df.deta)[i],"\nPhylogenetic Independent Contrasts"))
}

dev.off()
