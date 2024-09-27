
## A wrapper around cowplot ggsave2 to make handling unicode characters for PDF output simpler
savePlot <- function(filename,plot,device=cairo_pdf,family="Arial Unicode MS",...)
{
  ggsave2(filename,plot = plot,device=device,family=family,...)
}


cleanColumnNames <- function(col_name,list_of_option)
{
  ## find potential column name from list of options
  index <- unlist(purrr::map(col_name, 
                             function(x){
                               ## this species will match "metschnikowia_matae_maris.max.cds" and "metschnikowia_matae.max.cds", so force it to match the latter
                               if (x == "metschnikowia_matae")
                               {
                                 x <- paste0(x,".max.cds")
                               }
                               str_which(list_of_option,x)
                             }
  ))
  return(list_of_option[index])
}




cleanSpeciesLabels <- function(x)
{
  y <- x %>% str_remove(".max.cds") %>% 
    str_extract("[a-z]+_[a-z]+(_maris)*") %>%
    str_split("_") %>% 
    lapply(function(x){
      if (length(x) == 2)
      {
        paste0(str_to_sentence(x[1]),"_",x[2])
      } else {
        paste0(str_to_sentence(x[1]),"_",x[2],"_",x[3])
      }
    }
    )
  
  unlist(y)
}

cleanSpeciesLabelsForPlots <- function(x)
{
  y <- x %>% str_remove(".max.cds") %>% 
    str_extract("[a-z]+_[a-z]+") %>%
    str_split("_") %>% 
    lapply(function(x){paste0(toupper(substring(x[1],1,1)),". ",x[2])})
  #lapply(function(x){paste(str_to_sentence(x[1]),x[2])})
  unlist(y)
}


calculateAcrossSpeciesCorrelations <- function(df,tree,x.var,y.var,...) 
{
    x.var.sym <- as.name(x.var)
    y.var.sym <- as.name(y.var)
    tmp <- df %>% 
      column_to_rownames("Species") %>% 
      dplyr::select(!!x.var.sym,!!y.var.sym) %>% 
      dplyr::filter(!is.na(!!x.var.sym) & !is.na(!!y.var.sym)) %>%
      as.matrix()
    print(tmp)
    cleaned <- cleanData(tree,tmp)
    x <- cleaned$data
    tmp.tree <- cleaned$phy
    x.values <- x[,x.var]
    y.values <- x[,y.var]
    x.pic <- calculatePIC(x.values[tmp.tree$tip.label],tree = tmp.tree)
    y.pic <- calculatePIC(y.values[tmp.tree$tip.label],tree=tmp.tree)
    cor.test(x.pic[,1],y.pic[,1],...) %>% broom::tidy()
}

calculatePIC <- function(data,tree)
{
  
  pic.value <- pic(data,tree,var.contrasts = T)
  return(pic.value)
}

cleanData <- function(phy, data) {
  cleaned <- geiger::treedata(phy,data,warnings=T)
  return(cleaned)
}

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