setClass(Class = "MSIPFragmentMap",
         slots = list(
           "FG.atom.matrix" = "matrix",
           "FG.ratio.matrix" = "matrix",
           "FG.data" = "data.frame"
         ))
setClass(Class = "MSIPCoreData",
         slots = list(
           "Spectra_data" = "ANY",
           "FG_map" = "MSIPFragmentMap",
           "solve" = "list"
         ))

setMethod(f = "show",signature = "MSIPCoreData",
          definition = function(object){

            show(object@FG_map)
})


setMethod(f = "isEmpty",signature = "MSIPCoreData",
          definition = function(x){
            isEmpty(x@FG_map)
          })

#### MSIPFragmentMap



setMethod(f = "show",
          "MSIPFragmentMap",definition = function(object){

            x <- paste0("Map of ",
                        crayon::yellow(nrow(object@FG.atom.matrix)),
                        " fragment group, ",
                        crayon::yellow( ncol(object@FG.atom.matrix)),
                        " atom, ",
                        crayon::yellow(ncol(object@FG.ratio.matrix)),
                        " labeled")
            message(x)
          })

setMethod(f = "isEmpty",signature = "MSIPFragmentMap",
          definition = function(x){
  nrow(x@FG.atom.matrix)==0
})

#### MSIPIsotopomerMap
setClass(Class = "MSIPIsotopomerMap",
         slots = list(
           "isotopomer.defination" = "list",
           "isotopomer.map" = "matrix",
           "Labeled.FG.data" = "data.frame",
           "isotopomer.probability" = "numeric",
           "solve" = "list"
         ))


setMethod(f = "show",
          "MSIPIsotopomerMap",definition = function(object){

            x <- paste0("Map of ",
                        crayon::yellow(length(object@isotopomer.defination)),
                        " isotopomers, ",
                        crayon::yellow(nrow(object@isotopomer.map)),
                        " labeled fragment")
            message(x)

          })


setMethod(f = "isEmpty",signature = "MSIPIsotopomerMap",
          definition = function(x){
            length(x@isotopomer.defination)==0
          })

