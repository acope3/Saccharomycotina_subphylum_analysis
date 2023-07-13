library(ape)
library(phylogram)
library(cluster)
library(tidyverse)
library(factoextra)
library(dendextend)
library(gplots)
library(RColorBrewer)
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

converttRNAColnames <- function(trna)
{
  anti <- colnames(trna)[5:70]
  y<-strsplit(x=anti,split=".",fixed=T)
  codon.aa <- lapply(y,function(x){paste(reverseComplement(x[1]),x[2],sep=".")})
  colnames(trna)[5:70] <- codon.aa
  return(trna)
}


cleanTree <- function(tree,target.species)
{
  tips.to.drop <- tree$tip.label[which(!(tree$tip.label %in% target.species))]
  cleaned <- drop.tip(tree,tips.to.drop)
  return(cleaned)
}

createDataMatrix <- function(results.dir,target.runs,param="_Selection.csv" ,aa=NULL)
{
  runs <- list.dirs(results.dir,full.names = F,recursive = F)
  leu.runs <- runs[which(runs %in% target.runs)]
  if (is.null(aa))
  {
    csp.matrix <- matrix(rep(0,40*length(leu.runs)),nrow=40,ncol = length(leu.runs))
  } else{
    numCodons <- 0
    for (a in aa)
    {
      codons <- AAToCodon(a,T)
      numCodons <- numCodons + length(codons)
    }
    csp.matrix <- matrix(rep(0,numCodons*length(leu.runs)),nrow=numCodons,ncol = length(leu.runs))
  }
  csp.matrix <- as.data.frame(csp.matrix)
  colnames(csp.matrix) <- leu.runs
  
  for (i in leu.runs)
  {
    filename <- file.path(results.dir,i,"restart_4","Parameter_est",paste0(i,param))
    sel <- read.table(filename,header=T,stringsAsFactors = F,sep=",")

    sel <- sel[which(sel$Mean != 0),]
    if (!is.null(aa))
    {
      sel <- sel[which(sel$AA %in% aa),]
    }
    csp.matrix[,i] <- sel$Mean
    rownames(csp.matrix) <- sel$Codon
  }
  return(csp.matrix)
}

createDataMatrixAll <- function(results.dir,ser.runs,param="_Selection.csv",column = "Mean")
{
  runs <- list.dirs(results.dir,full.names = F,recursive = F)
  csp.matrix <- matrix(rep(0,40*length(runs)),nrow=40,ncol = length(runs))

  csp.matrix <- as.data.frame(csp.matrix)
  colnames(csp.matrix) <- runs
  aa <- aminoAcids()
  codons <- c()
  for (a in aa)
  {
    aa.codons <- AAToCodon(a,T)
    codons <- c(codons,aa.codons)
  }
  rownames(csp.matrix) <- codons
  for (i in runs)
  {
    if (!i %in% ser.runs)
    {
      restart.list <- sort(list.files(file.path(results.dir,i)))
      filename <- file.path(results.dir,i,restart.list[length(restart.list)],"Parameter_est",paste0(i,param))
      sel <- read.table(filename,header=T,stringsAsFactors = F,sep=",",row.names=2)
      sel <- sel[which(sel[,column] != 0),]
      
    } else{
      restart.list <- sort(list.files(file.path("/data/cope/Labella2019/Results_ser_reruns",i)))
      filename <- file.path("/data/cope/Labella2019/Results_ser_reruns",i,restart.list[length(restart.list)],"Parameter_est",paste0(i,param))
      sel <- read.table(filename,header=T,stringsAsFactors = F,sep=",",row.names=2)
      sel <- sel[which(sel[,column] != 0),]
      sel[40,"AA"] <- "J"
      sel[40,2:5] <- c(NA,NA,NA,NA)
    }
    sel <- sel[which(sel[,column] != 0),]
    csp.matrix[rownames(sel),i] <- sel[rownames(sel),column]
  }
  csp.matrix[csp.matrix == 0] <- NA
  return(csp.matrix)
}


createBKPlot <- function(dend.1,dend.2,...)
{
  Bk_plot(dend.1,dend.2,...)
}


compareHClustToPhylo <- function(phylo.dendro,csp.matrix,hclust.method,title="",...)
{
  cormat <- cor(t(csp.matrix),use="pairwise.complete.obs",method="spearman")
  hc <- hclust(as.dist(1-cormat),method = hclust.method)
  hc.dendro <- as.dendrogram(hc)
  dend.list <- dendlist(phylo.dendro,hc.dendro)
  corr <- cor.dendlist(dend.list,method="cophenetic")
  
  tanglegram(dend.list,highlight_distinct_edges = FALSE, # Turn-off dashed lines
             common_subtrees_color_lines = FALSE, # Turn-off line colors
             common_subtrees_color_branches = TRUE, # Color common branches
             main_left = paste(title,"\nCorrlelation =", round(corr[1,2], 2))
  )
  createBKPlot(phylo.dendro,hc.dendro,...)
}


# aa <- aminoAcids()

# codon.list <- vector(mode="list",length=3)

# for (a in aa)
# {
#   codons <- AAToCodon(a)
#   numCodons <- length(codons)
#   if (numCodons == 2){
#     codon.list[[1]] <- c(codon.list[[1]],a) 
#   } else if (numCodons == 4){
#     codon.list[[2]] <- c(codon.list[[2]],a) 
#   } else if (numCodons == 6){
#     codon.list[[3]] <- c(codon.list[[3]],a) 
#   }
# }

all.genomes <- read.table("../../all_fungi.txt",header=F,stringsAsFactors = F)
ser.genomes <- read.table("../../ser_fungi.txt",header=F,stringsAsFactors=F)
# deta <- createDataMatrixAll(results.dir = "/data/cope/Labella2019/Results/",ser.runs = ser.genomes[,1],param="_Selection.csv",column = "Mean")
# dm <- createDataMatrixAll(results.dir = "/data/cope/Labella2019/Results/",ser.runs = ser.genomes[,1],param="_Mutation.csv",column = "Mean")

# deta.sd <- createDataMatrixAll(results.dir = "/data/cope/Labella2019/Results/",ser.runs = ser.genomes[,1],param="_Selection.csv",column = "Std.Dev")
# dm.sd <- createDataMatrixAll(results.dir = "/data/cope/Labella2019/Results/",ser.runs = ser.genomes[,1],param="_Mutation.csv",column = "Std.Dev")


# write.table(deta,"../all_dEta_rerun_ser.tsv",sep="\t",col.names=T,row.names=T,quote=F)
# write.table(dm,"../all_dM_rerun_ser.tsv",sep="\t",col.names=T,row.names=T,quote=F)

# write.table(deta.sd,"../all_dEta_sd_rerun_ser.tsv",sep="\t",col.names=T,row.names=T,quote=F)
# write.table(dm.sd,"../all_dM_sd_rerun_ser.tsv",sep="\t",col.names=T,row.names=T,quote=F)


results.dir <- "/data/cope/Labella2019/Results/"
# 
fungi.tree <- read.tree("../../tree_with_cds_labels.nwk")
fungi.tree <- cleanTree(fungi.tree,all.genomes[,1])
fungi.dend <- as.dendrogram(fungi.tree)
# for (i in 1:3)
# {
#   if (i == 1)
#   {
#     suffix <- "2_codon_aa"
#   } else if (i == 2){
#     suffix <- "4_codon_aa"
#   } else if (i == 3){
#     suffix <- "6_codon_aa"
#   }
#   csp.matrix <- createDataMatrix(results.dir,leu.genomes[,1],param="_Selection.csv",aa=codon.list[[i]])
#   csp.matrix <- t(csp.matrix)
#   pdf(paste0("tanglegram_selection_",suffix,".pdf"),width = 12,height=20)
#   compareHClustToPhylo(fungi.dend.leu,csp.matrix,"single",title="HClust Algo: Single")
#   compareHClustToPhylo(fungi.dend.leu,csp.matrix,"complete",title="HClust Algo: Complete")
#   compareHClustToPhylo(fungi.dend.leu,csp.matrix,"average",title="HClust Algo: Average")
#   compareHClustToPhylo(fungi.dend.leu,csp.matrix,"ward.D",title="HClust Algo: Ward.D")
#   compareHClustToPhylo(fungi.dend.leu,csp.matrix,"ward.D2",title="HClust Algo: Ward.D2")
#   dev.off()
#   
#   csp.matrix <- createDataMatrix(results.dir,leu.genomes[,1],param="_Mutation.csv",aa=codon.list[[i]])
#   csp.matrix <- t(csp.matrix)
#   pdf(paste0("tanglegram_mutation_",suffix,".pdf"),width = 12,height=20)
#   compareHClustToPhylo(fungi.dend.leu,csp.matrix,"single",title="HClust Algo: Single")
#   compareHClustToPhylo(fungi.dend.leu,csp.matrix,"complete",title="HClust Algo: Complete")
#   compareHClustToPhylo(fungi.dend.leu,csp.matrix,"average",title="HClust Algo: Average")
#   compareHClustToPhylo(fungi.dend.leu,csp.matrix,"ward.D",title="HClust Algo: Ward.D")
#   compareHClustToPhylo(fungi.dend.leu,csp.matrix,"ward.D2",title="HClust Algo: Ward.D2")
#   dev.off()
#   
# }
# 
csp.matrix <- read.table("../all_dEta_rerun_ser.tsv",sep="\t",header=T,stringsAsFactors=F,row.names=1)
csp.matrix <- t(csp.matrix)
csp.matrix <- csp.matrix[fungi.tree$tip.label,]
pdf(paste0("Images/tanglegram_selection_all_fungi.pdf"),width = 12,height=20)
compareHClustToPhylo(fungi.dend,csp.matrix,"single",title="HClust Algo: Single",main="Bk plot:\nSingle")
compareHClustToPhylo(fungi.dend,csp.matrix,"complete",title="HClust Algo: Complete",main="Bk plot:\nComplete")
compareHClustToPhylo(fungi.dend,csp.matrix,"average",title="HClust Algo: Average",main="Bk plot:\nAverage")
compareHClustToPhylo(fungi.dend,csp.matrix,"ward.D",title="HClust Algo: Ward.D",main="Bk plot:\nWard.D")
compareHClustToPhylo(fungi.dend,csp.matrix,"ward.D2",title="HClust Algo: Ward.D2",main="Bk plot:\nWard.D2")
dev.off()


csp.matrix <- read.table("../all_dM_rerun_ser.tsv",sep="\t",header=T,stringsAsFactors=F,row.names=1)
csp.matrix <- t(csp.matrix)
csp.matrix <- csp.matrix[fungi.tree.$tip.label,]
pdf(paste0("Images/tanglegram_mutation_all_fungi.pdf"),width = 12,height=20)
compareHClustToPhylo(fungi.dend,csp.matrix,"single",title="HClust Algo: Single",main="Bk plot:\nSingle")
compareHClustToPhylo(fungi.dend,csp.matrix,"complete",title="HClust Algo: Complete",main="Bk plot:\nComplete")
compareHClustToPhylo(fungi.dend,csp.matrix,"average",title="HClust Algo: Average",main="Bk plot:\nAverage")
compareHClustToPhylo(fungi.dend,csp.matrix,"ward.D",title="HClust Algo: Ward.D",main="Bk plot:\nWard.D")
compareHClustToPhylo(fungi.dend,csp.matrix,"ward.D2",title="HClust Algo: Ward.D2",main="Bk plot:\nWard.D2")
dev.off()

# # 
# # 
to.remove <- c("CAG.Ala","CAG.Ser","TTA.Stop","TCA.Stop","CTA.Stop")
trna <- read.table("../Labella2019_tRNA_per_genome.tsv",sep="\t",header=T,stringsAsFactors = F,row.names = 1)
trna.filt <- trna[,which(!colnames(trna) %in% to.remove)]
trna.filt <- trna.filt[fungi.tree.leu$tip.label,]
dimm <- dim(trna.filt)
tgcn <- trna.filt[,4:dimm[2]]

tgcn.hc<- hclust(dist(tgcn),method = "complete")
tgcn.dendro <- as.dendrogram(tgcn.hc)


# # 
# # pdf("tanglegram_trna.pdf",width=12,height = 20)
# # compareHClustToPhylo(fungi.dend.leu,tgcn,"single",title="HClust Algo: Single")
# # compareHClustToPhylo(fungi.dend.leu,tgcn,"complete",title="HClust Algo: Complete")
# # compareHClustToPhylo(fungi.dend.leu,tgcn,"average",title="HClust Algo: Average")
# # compareHClustToPhylo(fungi.dend.leu,tgcn,"ward.D",title="HClust Algo: Ward.D")
# # compareHClustToPhylo(fungi.dend.leu,tgcn,"ward.D2",title="HClust Algo: Ward.D2")
# # dev.off()
# # 
# pdf("heatmap_trna.pdf",width=20,height=20)
# createCorrelationHeatMap(tgcn,fungi.tree.leu,main="Comparing tRNA Gene Copy Numbers",cexRow = 0.3,cexCol=0.3,margins=c(10,10))
# dev.off()


# phi <- read.table("../phi_matrix_all_fungi.tsv",sep="\t",header=T,row.names = 1,stringsAsFactors = F)
# phi <- phi[,fungi.tree$tip.label]
# phi <- t(phi)


# pdf(paste0("Images/tanglegram_phi_all_fungi.pdf"),width = 12,height=20)
# compareHClustToPhylo(fungi.dend,phi,"single",title="HClust Algo: Single",main="Bk plot:\nSingle")
# compareHClustToPhylo(fungi.dend,phi,"complete",title="HClust Algo: Complete",main="Bk plot:\nComplete")
# compareHClustToPhylo(fungi.dend,phi,"average",title="HClust Algo: Average",main="Bk plot:\nAverage")
# compareHClustToPhylo(fungi.dend,phi,"ward.D",title="HClust Algo: Ward.D",main="Bk plot:\nWard.D")
# compareHClustToPhylo(fungi.dend,phi,"ward.D2",title="HClust Algo: Ward.D2",main="Bk plot:\nWard.D2")
# dev.off()
# pdf("Images/heatmap_phi_all_fungi.pdf",width=20,height=20)
# createCorrelationHeatMap(phi,fungi.tree,main="Comparing Evolutionary Average Protein Production Rates\nPairwise Single-copy Orthologs",cexRow = 0.3,cexCol=0.3,margins=c(10,10))
# dev.off()


# csp.matrix <- read.table("leu_dEta.tsv",sep="\t",header=T,stringsAsFactors=F,row.names=1)
# csp.matrix <- t(csp.matrix)
# csp.matrix <- csp.matrix[fungi.tree.leu$tip.label,]
# 
# deta.hc<- hclust(dist(csp.matrix),method = "complete")
# deta.dendro <- as.dendrogram(deta.hc)
# 
# dend.list <- dendlist(deta.dendro,phi.dendro)
# corr <- cor.dendlist(dend.list,method="cophenetic")
# pdf("tanglegram_phi_vs_deta.pdf",width = 15,height=20)
# tanglegram(dend.list,highlight_distinct_edges = FALSE, # Turn-off dashed lines
#            common_subtrees_color_lines = FALSE, # Turn-off line colors
#            common_subtrees_color_branches = TRUE, # Color common branches
#            main_left = paste("Comparing Selection Coefficients and Phi\nCorrlelation =", round(corr[1,2], 2)))
# dev.off()
# 
# 
# dend.list <- dendlist(deta.dendro,tgcn.dendro)
# corr <- cor.dendlist(dend.list,method="cophenetic")
# pdf("tanglegram_tgcn_vs_deta.pdf",width = 15,height=20)
# tanglegram(dend.list,highlight_distinct_edges = FALSE, # Turn-off dashed lines
#            common_subtrees_color_lines = FALSE, # Turn-off line colors
#            common_subtrees_color_branches = TRUE, # Color common branches
#            main_left = paste("Comparing Selection Coefficients and tGCN\nCorrlelation =", round(corr[1,2], 2)))
# dev.off()
# 

# deta.matrix <- read.table("leu_dEta.tsv",sep="\t",header=T,stringsAsFactors=F,row.names=1)
# deta.matrix <- t(deta.matrix)
# deta.matrix <- deta.matrix[fungi.tree.leu$tip.label,]
# 
# dm.matrix <- read.table("leu_dM.tsv",sep="\t",header=T,stringsAsFactors=F,row.names=1)
# dm.matrix <- t(dm.matrix)
# dm.matrix <- dm.matrix[fungi.tree.leu$tip.label,]
# 
# csp.all <- cbind(as.data.frame(deta.matrix),as.data.frame(dm.matrix))
# 
# pdf("tanglegram_csp_all.pdf",width = 12,height=20)
# compareHClustToPhylo(fungi.dend.leu,csp.all,"single",title="HClust Algo: Single")
# compareHClustToPhylo(fungi.dend.leu,csp.all,"complete",title="HClust Algo: Complete")
# compareHClustToPhylo(fungi.dend.leu,csp.all,"average",title="HClust Algo: Average")
# compareHClustToPhylo(fungi.dend.leu,csp.all,"ward.D",title="HClust Algo: Ward.D")
# compareHClustToPhylo(fungi.dend.leu,csp.all,"ward.D2",title="HClust Algo: Ward.D2")
# dev.off()
# 
# 