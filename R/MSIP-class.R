

setClass(Class = "MSIPIsotopologueData",
         contains = "list")


MSIPIsotopologueData <- function(...){


  x <- new("MSIPIsotopologueData")

  input <- list(...)
  for (i in seq_along(input)) {
    x[names(input)] <- input[[i]]
  }
  return(x)

}


setMethod(f = "show",signature = "MSIPIsotopologueData",
          definition = function(object){
            xn <- names(object)
            xn <- xn[grepl("M",xn)]
            l <- paste0(format_isotopologue(xn,"+"),collapse = ";")
            cat(length(xn),"isotopologues: ",l )

          })


setClass(Class = "MSIPMetaboliteData",
         slots = list(
           "CompoundInfo" = "list",
           "Spectra" = "list",
           "MSIPIsotopologueDatas" = "MSIPIsotopologueData"
         ))

MSIPMetaboliteData <- function(CompoundInfo = list(),
                               Spectra = list(),
                               MSIPIsotopologueDatas = MSIPIsotopologueData()
){

  new("MSIPMetaboliteData",
      CompoundInfo = CompoundInfo,
      Spectra = Spectra,
      MSIPIsotopologueDatas = MSIPIsotopologueDatas
  )


}

setMethod(f = "show",signature = "MSIPMetaboliteData",definition = function(object){

  cat("MSIPMetaboliteData")
})
