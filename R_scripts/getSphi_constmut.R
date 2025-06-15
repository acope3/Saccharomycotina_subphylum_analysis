library(AnaCoDa)

samples <- 20000
all.genomes <- readLines("../Data/all_fungi.txt")
improved.fits.species <- readLines("../Data/improved_species_by_gene_expression_25_percent_2023_05_03.txt")

for (species in improved.fits.species)
{
  if (species %in% improved.fits.species)
  {
    results.dir <- "Results_k_1"
  } else {
    results.dir <- "Results_k_1"
  }
  # param <- loadParameterObject(file.path("..","Final_runs",results.dir,species,"run_2","R_objects","parameter.Rda"))
  # mean.sphi <- param$getStdDevSynthesisRatePosteriorMean(samples,0)
  # write(mean.sphi,file = file.path("..","Final_runs",results.dir,species,"run_2","Parameter_est","sphi.txt"),ncolumns = 1)
  # 
  tryCatch({
    param <- loadParameterObject(file.path("..","Final_runs",results.dir,species,"run_2","R_objects","parameter.Rda"))
    mean.sphi <- param$getStdDevSynthesisRatePosteriorMean(samples,0)
    write(mean.sphi,file = file.path("..","Final_runs",results.dir,species,"run_2","Parameter_est","sphi.txt"),ncolumns = 1)
  }, error=function(e){print(species)},
  warning = function(w){print(species)}
  )
}
