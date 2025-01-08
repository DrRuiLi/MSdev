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

setGeneric("atom",
           def = function(object,
                          element = element_table$element){
             rownames(atomblock(object))[bonds(object)$atom%in%element ]
           }
)

setGeneric("element",
           def = function(object,...){
             bonds(object)$atom
           }
)
