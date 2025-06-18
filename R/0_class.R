setGeneric(
  "vdata",
  def = function(object) {
    igraph::as_data_frame(object, "vertices")
  }
)

setGeneric(
  "edata",
  def = function(object) {
    igraph::as_data_frame(object, "edges")
  }
)



setGeneric(
  "vdata<-",
  def = function(object, value) {
    igraph::vertex.attributes(object) <- as.list(value)
    object
  }
)
setGeneric(
  "edata<-",
  def = function(object, value) {
    value <- value[, !grepl("^from|to$", colnames(value))]
    igraph::edge.attributes(object) <- as.list(value)
    object
  }
)

#' atom
#'
#' @export
#'
setGeneric("atom",
           def = function(object,
                          element = "ANY"){

             if (element=="ANY") {
               data(element_table)
               element <- element_table$element
             }
             rownames(atomblock(object))[bonds(object)$atom%in%element ]
           }
)

setGeneric("element",
           def = function(object,...){
             bonds(object)$atom
           }
)


#' Example dataset: element_table
#'
#' This is a description of your dataset.
#'
#' @format A data frame with N rows and M columns.
#' @source Simulated data.
"element_table"
