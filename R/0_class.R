
#if (!isGeneric("vdata"))
  setGeneric(
  "vdata",
  def = function(object) {
    igraph::as_data_frame(object, "vertices")
  }
)

#if (!isGeneric("edata"))
setGeneric(
  "edata",
  def = function(object) {
    igraph::as_data_frame(object, "edges")
  }
)



#if (!isGeneric("vdata<-"))
setGeneric(
  "vdata<-",
  def = function(object, value) {
    igraph::vertex.attributes(object) <- as.list(value)
    object
  }
)
#if (!isGeneric("edata<-"))
setGeneric(
  "edata<-",
  def = function(object, value) {
    value <- value[, !grepl("^from$|^to$", colnames(value))]
    igraph::edge.attributes(object) <- as.list(value)
    object
  }
)

#' atom
#'
#' @export
#'
#if (!isGeneric("atom"))
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

#if (!isGeneric("element"))



  setGeneric("get_element",
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
