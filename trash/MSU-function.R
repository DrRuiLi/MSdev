saveMSU <- function(object){
  MSU.obj <- object
  save.dir <- dirname(object@projectInfo$MSUFile)
  if (!dir.exists(save.dir)) {
    dir.create(save.dir,recursive = T)
  }
  save(MSU.obj, file =  object@projectInfo$MSUFile)
  invisible(MSU.obj)
}
