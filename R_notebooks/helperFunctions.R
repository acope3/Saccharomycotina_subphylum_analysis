
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




calculateAcrossSpeciesCorrelations <- function(df,tree,x.var,y.var,approach="standard",...) 
{
    x.var.sym <- as.name(x.var)
    y.var.sym <- as.name(y.var)
    tmp <- df %>% 
      column_to_rownames("Species") %>% 
      dplyr::select(!!x.var.sym,!!y.var.sym) %>% 
      dplyr::filter(!is.na(!!x.var.sym) & !is.na(!!y.var.sym)) %>%
      as.matrix()
    if (nrow(tmp) == 0)
    {
      return(NA)
    } else {
    tryCatch({
        cleaned <- cleanData(tree,tmp)
        x <- cleaned$data
        tmp.tree <- cleaned$phy
        x.values <- x[,x.var]
        y.values <- x[,y.var]
        x.values <- x.values[tmp.tree$tip.label]
        y.values <- y.values[tmp.tree$tip.label]
        if (approach == "standard")
        {
          x.pic <- calculatePIC(x.values,tree = tmp.tree)
          y.pic <- calculatePIC(y.values,tree=tmp.tree)
          df <- data.frame(X = x.pic[,1],Y = y.pic[,1])
          return(cor.test(df$X,df$Y,...) %>% broom::tidy())
        } else {
          result <- picRegression(tmp.tree,x.values,y.values,method="rank", sigTest="permutation",nperm=9999)
          df <- data.frame(estimate=result$testStat,p.value=result$pVal)
          return(df)
        }
    }, error=function(e){return(NA)})
    }
}

picRegression <- function (tree, x, y, method="standard", sigTest="permutation",
                           ...) 
{
  ## original version was missing nperm argument so I added it (LR)
  if(hasArg(nperm)) nperm<-list(...)$nperm
  else nperm<-100
  if (!inherits(tree, "phylo")) 
    stop("tree should be object of class \"phylo\".")
  
  if (method %in% c("standard", "sign", "rank") == FALSE) {
    cat("  Invalid model. Setting model=\"standard\"\n\n")
    model <- "standard"
  }
  
  if (sigTest %in% c("analytic", "permutation") == FALSE) {
    cat("  Invalid model. Setting model=\"permutation\"\n\n")
    sigTest <- "permutation"
  }
  
  if((method=="sign" | method=="rank")&sigTest=="analytic") {
    cat("  No analytic p-value method exists for method ", method, ". Setting model=\"permutation\"\n\n")
    sigTest <- "permutation"
  }
  
  if(method=="standard") 
  {
    icx<-ape::pic(x, tree)
    icy<-ape::pic(y, tree)
    
    if(sigTest=="analytic") 
    {
      res<-summary(lm(icy~icx+0))
      testStat<-res$coefficients[1,3]
      pVal<-res$coefficients[1,4]
    }
    if(sigTest=="permutation") 
    {
      dicx<-c(icx, -icx)	
      dicy<-c(icy, -icy)
      res<-summary(lm(dicy~dicx))
      testStat<-res$coefficients[2,3]
      nullDist<-numeric(nperm)
      for(i in 1:nperm) 
      {
        picx<-sample(icx)
        dpicx<-c(picx, -picx)
        res<-summary(lm(dicy~dpicx))
        nullDist[i]<-res$coefficients[2,3]
      }
      
      pValueHigh<-2*(sum(nullDist >= testStat)+1)/(nperm+1)
      pValueLow<-2*(sum(nullDist <= testStat)+1)/(nperm+1)
      
      pVal<-min(pValueHigh, pValueLow)
    }
  }
  
  if(method=="sign") 
  {
    icx<-ape::pic(x, tree)
    icy<-ape::pic(y, tree)
    xPos<-icx>0
    yPos<-icy>0
    testStat<-sum(xPos==yPos)
    pValueLow<-2*pbinom(testStat, length(icx), prob=0.5)
    pValueHigh<-2*pbinom(length(icx)-testStat, length(icx), prob=0.5)
    pVal<-min(pValueHigh, pValueLow)
  }
  
  if(method=="rank") 
  {
    icx<-ape::pic(x, tree)
    icy<-ape::pic(y, tree)
    dicx<-c(icx, -icx)	
    dicy<-c(icy, -icy)
    
    rx<-rank(dicx, ties.method="average")
    ry<-rank(dicy, ties.method="average")
    res<-summary(lm(ry~rx))
    testStat<-res$coefficients[2,1]
    testStat <- cor(rx,ry)
    nullDist<-numeric(nperm)
    for(i in 1:nperm) 
    {
      picx<-sample(icx)
      dpicx<-c(picx, -picx)
      prx<-rank(dpicx, ties.method="average")
      res<-summary(lm(ry~prx))
      nullDist[i]<-res$coefficients[2,1]
      #nullDist[i]<- cor(prx,ry)
    }
    
    pValueHigh<-2*(sum(nullDist >= testStat)+1)/(nperm+1)
    pValueLow<-2*(sum(nullDist <= testStat)+1)/(nperm+1)
    
    pVal<-min(pValueHigh, pValueLow)
    
  }
  obj <- list(testStat=testStat, pVal=pVal)
  class(obj) <- "picRegression"
  obj
}


calculatePIC <- function(data,tree)
{
  
  pic.value <- pic(data,tree,var.contrasts = T)
  return(pic.value)
}

cleanData <- function(phy, data) {
  cleaned <- geiger::treedata(phy,data,warnings=F)
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