
getPCADataFrom_ropls <- function(obj.rops){

  x <- data.frame(obj.rops@scoreMN)

  return(x)
}



plot_cpd_pathview <- function(cpd,
                              pathway.id,
                              dir.to.save ,
                              kegg.dir = tempdir(),
                              ...){

  library(pathview)
  wd <- getwd()
  setwd(dir.to.save)
  pathview::pathview(
    cpd.data = cpd,pathway.id = pathway.id,
    kegg.dir = kegg.dir,
    high = list(gene = "red", cpd =
                  "red"),
           ...)
  setwd(wd)

}

