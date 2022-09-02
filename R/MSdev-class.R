setClass("MSdev",
         slots = list(
           projectInfo = "data.frame",
           processingInfo = "list",
           sampleInfo = "data.frame",
           experimentInfo ="MS_Exp",
           xcmsData = "list",
           spectra = "list",
           annotation = "list",
           statData = "list"
           ))



MSdev <- function()
  new("MSdev")


setMethod("initialize" , "MSdev",
          function(.Object){
            .Object@projectInfo = data.frame()
            .Object

          })



