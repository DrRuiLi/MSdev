
setOldClass("object_size")
setClass("onDiskData",
         slots = list("path" = "character",
                      "size" = "object_size"))

onDiskData <- function(data,path = tempfile()){


  path <- normalizePath(path,winslash = "/",mustWork = F)
  data.size <- object.size(data)
  #data.size <- as.numeric(data.size)/1024^2
  saveRDS(data,file = path)
  new("onDiskData",
      path=path,
      size = data.size)
}
setMethod("show","onDiskData",
          definition = function(object){

            message("onDiskData with size of ",format(object@size,units = "MB"))
            message("Store in ",format(object@path))
          })


onDiskData_retrieve <- function(object){

  if (class(object)=="onDiskData") {
    data <- readRDS(object@path)
    return(data)
  }else{
    return(object)
  }
}
