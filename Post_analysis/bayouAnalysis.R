library(bayou)
library(parallel)

cleanTree <- function(tree,target.species)
{
  tips.to.drop <- tree$tip.label[which(!(tree$tip.label %in% target.species))]
  cleaned <- drop.tip(tree,tips.to.drop)
  return(cleaned)
}

ou.analysis <- function(trait,sd.trait,tree,codon)
{
  priorOU <- make.prior(tree, 
                        dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                   dk="cdpois", dtheta="dnorm"),
                        param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                   dk=list(lambda=10, kmax=50), dsb=list(bmax=1, prob=1), 
                                   dtheta=list(mean=mean(trait), sd=1.5*sd(trait))))
  startpars <- priorSim(priorOU, tree, plot=TRUE)$pars[[1]]
  priorOU(startpars)
  dir.create(file.path("/data2","Labella2019","Phylogenetic_analysis",paste0(codon,"mutation_MCMC_objects")))
  mcmcOU <- bayou.makeMCMC(tree, trait, prior=priorOU,SE=sd.trait, 
                           new.dir=file.path("/data2","Labella2019","Phylogenetic_analysis",paste0(codon,"mutation_MCMC_objects")), outname=paste0("r0001_mcmc_deta_",codon), plot.freq=NULL) # Set up the MCMC
  mcmcOU$run(4000000) # Run the MCMC
  
  pdf(paste0(codon,".pdf"))
  chainOU <- mcmcOU$load()
  
  chainOU <- set.burnin(chainOU, 0.001)
  summary(chainOU)
  plot(chainOU, auto.layout=FALSE)
  plotSimmap.mcmc(chainOU, burnin = 0.001, pp.cutoff = 0.001,cex=0.02)
  plotBranchHeatMap(tree, chainOU, "theta", burnin = 0.001, pal = cm.colors,cex=0.02)
  phenogram.density(tree, trait, burnin = 0.001, chainOU, pp.cutoff = 0.3)
  dev.off()
  saveRDS(mcmcOU,file=paste0("mutation_mcmc_",codon,".rds"))
}


all.genomes <- read.table("../all_fungi.txt",header=F,stringsAsFactors = F)

fungi.tree <- read.tree("../tree_with_cds_labels.nwk")
fungi.tree <- cleanTree(fungi.tree,all.genomes[,1])

deta <- read.table("all_dM_rerun_ser.tsv",sep="\t",header=T,stringsAsFactors = F,row.names = 1)
deta.t <- t(deta)

deta.sd <- read.table("all_dM_sd_rerun_ser.tsv",sep="\t",header=T,stringsAsFactors = F,row.names = 1)
deta.sd.t <- t(deta.sd)

deta.t <- deta.t[fungi.tree$tip.label,]
deta.sd.t <- deta.sd.t[fungi.tree$tip.label,]

ctg <- which(colnames(deta.t)=="CTG")
deta.t <- deta.t[,-ctg]
deta.sd.t <- deta.sd.t[,-ctg]


mclapply(colnames(deta.t),
         function(x){
           ou.analysis(trait=deta.t[,x],sd=deta.sd.t[,x],tree=fungi.tree,codon=x)
           }
         ,mc.cores=39)