
saveMSdev <- function(object){
  MSdev <- object
  save(MSdev, file =  object@projectInfo$MSdevFile)
}
