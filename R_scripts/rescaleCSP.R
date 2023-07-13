library(AnaCoDa,lib.loc="~/R_dev")
library(tidyverse)


changeCodon <- function(codon)
{
  nucleotides <- unlist(strsplit(codon,""))
  if (nucleotides[3] == "G")
  {
    codon <- paste0(nucleotides[1],nucleotides[2],nucleotides[3],"|A",collapse="")
  } else if (nucleotides[3] == "C"){
    codon <- paste0(nucleotides[1],nucleotides[2],nucleotides[3],"|T",collapse="")
  }
  return(codon)
}


rescaleAGCT <- function(sel.file,parameter.file,genome.file,codon.table=1,param.type=0,mixture=1,percent.to.keep=0.75)
{
  data <- read_csv(sel.file,sep=",",header=T,stringsAsFactors=F,row.names = 2)
  new.data <- data.frame("AA"=c(),"Mean"=c(),"Std.Dev"=c(),"X0.025."=c(),"X0.975."=c())
  parameter <- loadParameterObject(parameter.file)
  trace <- parameter$getTraceObject()
  aa <- aminoAcids()
  genome <- initializeGenomeObject(genome.file,codon_table=codon.table)
  for (a in aa)
  {
    if (a == "X" || a=="M" || a == "W") next
    data.tmp <- data[which(data$AA == a),]
    if (a == "R")
    {
      trace.aga <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,"AGA",param.type,T)
      trace.agg <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,"AGG",param.type,T)
      trace.cga <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,"CGA",param.type,T)
      trace.cgc <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,"CGC",param.type,T)
      trace.cgg <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,"CGG",param.type,T)
      samples <- length(trace.cgc)
      trace.agg.a <- trace.agg[(samples * percent.to.keep):samples] - trace.aga[(samples * percent.to.keep):samples]
      trace.cgg.a <- trace.cgg[(samples * percent.to.keep):samples] - trace.cga[(samples * percent.to.keep):samples]
      trace.cgc.t <- trace.cgc.t 
      data.tmp["AGG",c("AA","Mean","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.agg.a),sd(trace.agg.a),quantile(trace.agg.a,c(0.025,0.975)))
      data.tmp["CGG",c("AA","Mean","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.cgg.a),sd(trace.cgg.a),quantile(trace.cgg.a,c(0.025,0.975)))
      data.tmp["CGC",c("AA","Mean","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.cgc.t),sd(trace.cgc.t),quantile(trace.cgc.t,c(0.025,0.975)))
      data.tmp["CGT",] <- NA
      data.tmp["CGA",] <- NA
      data.tmp["AGA",] <- NA
    }
    if (a == "L")
    {
      trace.ctc <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,"CTC",param.type,T)
      trace.ctt <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,"CTT",param.type,T)
      trace.tta <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,"TTA",param.type,T)
      trace.ctc.t <- trace.ctc[(samples * percent.to.keep):samples] - trace.ctt[(samples * percent.to.keep):samples]
      trace.ttg.a <- -1*trace.tta[(samples * percent.to.keep):samples]
      data.tmp["CTC",c("AA","Mean","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.ctc.t),sd(trace.ctc.t),quantile(trace.ctc.t,c(0.025,0.975)))
      data.tmp["TTG",c("AA","Mean","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.ttg.a),sd(trace.ttg.a),quantile(trace.ttg.a,c(0.025,0.975)))
      data.tmp["CTA",] <- NA
      data.tmp["CTT",] <- NA
      data.tmp["TTA",] <- NA
      if (codon.table != 12) 
      {
        trace.cta <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,"CTA",param.type,T)
        trace.ctg <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,"CTG",param.type,T)
        trace.cta.g <- trace.ctg[(samples * percent.to.keep):samples] - trace.cta[(samples * percent.to.keep):samples]
        data.tmp["CTG",c("AA","Mean","Std.Dev","X0.025.","X0.975.")] <-  c(a,mean(trace.ctg.a),sd(trace.ctg.a),quantile(trace.ctg.a,c(0.025,0.975)))
      }
      
    }
    if (a == "I")
    {
      data.tmp["ATC",] <- data.tmp["ATC",]
      data.tmp["ATA", ] <- NA
      data.tmp["ATT", ] <- NA
    }
    if (nrow (data.tmp) == 4)
    {
      codon.a <- row.names(data.tmp)[1]
      codon.c <- row.names(data.tmp)[2]
      codon.g <- row.names(data.tmp)[3]
      codon.t <- row.names(data.tmp)[4]
      trace.a <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,codon.a,param.type,T)
      trace.g <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,codon.g,param.type,T)
      trace.g.a<- trace.g[(samples * percent.to.keep):samples] - trace.a[(samples * percent.to.keep):samples]
      data.tmp[codon.g,c("AA","Mean","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.g.a),sd(trace.g.a),quantile(trace.g.a,c(0.025,0.975)))
      data.tmp[codon.c,] <- data.tmp[codon.c,]
      data.tmp[codon.a,] <- NA
      data.tmp[codon.t,] <- NA
    }
    if (nrow (data.tmp) == 2)
    {
      codon.1 <- row.names(data.tmp)[1]
      codon.2 <- row.names(data.tmp)[2]
      if (substring(codon.1,first = 3,last = 3) == "A")
      {
        trace.a <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture,codon.a,param.type,T)
        trace.g.a<- -1 * trace.a[(samples * percent.to.keep):samples]
        data.tmp[codon.2,c("AA","Mean","Std.Dev","X0.025.","X0.975.")] <- c(a,mean(trace.g.a),sd(trace.g.a),quantile(trace.g.a,c(0.025,0.975)))
        data.tmp[codon.1,] <- NA
      } else{
        data.tmp[codon.1,] <- data.tmp[codon.1,]
        data.tmp[codon.2,] <- NA
      }
    }
    data.tmp <- data.tmp[which(!is.na(data.tmp$Mean)),]
    new.data <- rbind(new.data,data.tmp)
  }
  
  new.data[,"Codon"] <- row.names(new.data)
  new.data <- new.data[,c("AA","Codon","Mean","Std.Dev","X0.025.","X0.975.")]
  return(new.data)
}


head.directory <- "../Final_runs/Results_k_1/"
species <- list.files()
targets <- c("Secondary_structure_paired_est_dM")

for (f in targets)
{
  structure.loc <- file.path(head.directory,f)
  structures <- list.dirs(structure.loc,recursive=F)
  for (struct in structures)
  {
    sel.file.loc <- file.path(struct,"restart_5","Parameter_est/")
    sel.file <- list.files(sel.file.loc,pattern="*_Selection.csv",full.names=T)
    parameter.file <- file.path(struct,"restart_5","R_objects","parameter.Rda")
    print(sel.file.loc)
    sel.rescaled <- rescaleAGCT(sel.file,genome.file=genome.file)
    write.table(sel.rescaled,file.path(sel.file.loc,"selection_rescaled_to_genome_optimal.csv"),sep=",",col.names=T,row.names=F,quote=F)
  }	
}