treeTraitPlot <- function(tree.plot,trait.df,trait.name,panel.name,significance = NULL,xlim=c(-1,1))
{
  
  if (!is.null(significance))
  {
    trait.df$Significance <- ifelse(trait.df[,significance] < 0.05,"Significant","Not Significant")
    trait.plot <- facet_plot(tree.plot,panel=panel.name,data=trait.df,geom=geom_segment,aes(x=0,xend=!!as.name(trait.name),y=y,yend=y,color=Clade,linetype=Significance))
  } else{
    trait.plot <- facet_plot(tree.plot,panel=panel.name,data=trait.df,geom=geom_segment,aes(x=0,xend=!!as.name(trait.name),y=y,yend=y,color=Clade))
  }
  trait.plot <- trait.plot +
    xlim_expand(xlim, panel.name) +
    scale_linetype_manual(values=c("Significant"=1,"Not Significant"=2))
  return(trait.plot)
}

cleanTree <- function(tree,target.species)
{
  tips.to.drop <- tree$tip.label[which(!(tree$tip.label %in% target.species))]
  cleaned <- drop.tip(tree,tips.to.drop)
  return(cleaned)
}