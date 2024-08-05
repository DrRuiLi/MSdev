setClass(Class = "MSIPCoreData",
         slots = list(
           "Spectra_data" = "data.frame",
           "FG_map" = "MSIPFragmentMap",
           "solve" = "list"
         ))

setMethod(f = "show",signature = "MSIPCoreData",
          definition = function(object){

            show(object@FG_map)
})


setMethod(f = "isEmpty",signature = "MSIPCoreData",
          definition = function(object){
            isEmpty(object@FG_map)
          })

#### MSIPFragmentMap
setClass(Class = "MSIPFragmentMap",
         slots = list(
           "fragment.atom.matrix" = "matrix",
           "fragment.ratio.matrix" = "matrix",
           "fragment.intensity" = "numeric",
           "fragment.include" = "logical"
         ))


setMethod(f = "show",
          "MSIPFragmentMap",definition = function(object){

            x <- paste0("Map of ",
                        yellow(nrow(object@fragment.atom.matrix)),
                        " fragment group, ",
                        yellow( ncol(object@fragment.atom.matrix)),
                        " atom, ",
                        yellow(ncol(object@fragment.ratio.matrix)),
                        " labeled")
            message(x)
          })


setMethod(f = "isEmpty",signature = "MSIPFragmentMap",
          definition = function(object){
  nrow(object@fragment.atom.matrix)==0
})

#### MSIPIsotopomerMap
setClass(Class = "MSIPIsotopomerMap",
         slots = list(
           "isotopomer.defination" = "list",
           "isotopomer.map" = "matrix",
           "isotopomer.ratio" = "numeric",
           "isotopomer.probability" = "numeric",
           "isotopomer.intensity" = "numeric",
           "solve" = "list"
         ))


setMethod(f = "show",
          "MSIPIsotopomerMap",definition = function(object){

            x <- paste0("Map of ",
                        yellow(length(object@isotopomer.defination)),
                        " iso-forms, ",
                        yellow(nrow(object@isotopomer.map)),
                        " labeled fragment")
            message(x)

          })


setMethod(f = "isEmpty",signature = "MSIPIsotopomerMap",
          definition = function(object){
            length(object@isotopomer.defination)==0
          })

