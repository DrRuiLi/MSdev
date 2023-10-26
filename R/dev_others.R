
getPCADataFrom_ropls <- function(obj.rops){

  x <- data.frame(obj.rops@scoreMN)

  return(x)
}



plot_cpd_pathview <- function(cpd,
                              pathway.id,
                              dir.to.save ){

  wd <- getwd()
  setwd(dir.to.save)
  pathview(cpd.data = cpd,pathway.id = pathway.id,
           kegg.dir = "D:/pathview/")
  setwd(wd)

}

