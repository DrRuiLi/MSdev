#' Import graph generics from MSCC after chemistry migration
#'
#' @importFrom MSCC atom
#' @importFrom MSCC edata
#' @importFrom MSCC get_element
#' @importFrom MSCC vdata
#' @rawNamespace importFrom(MSCC,`edata<-`)
#' @rawNamespace importFrom(MSCC,`vdata<-`)
#' @name mscc-graph-imports
#' @keywords internal
NULL

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  methods::setMethod(
    "vdata",
    signature = "Metabolic_flux_network",
    definition = function(object) {
      vdata(object@metabolic_network)
    },
    where = ns
  )
  methods::setMethod(
    "vdata<-",
    "Metabolic_flux_network",
    definition = function(object, value) {
      vdata(object@metabolic_network) <- value
      object
    },
    where = ns
  )
  methods::setMethod(
    "edata",
    signature = "Metabolic_flux_network",
    definition = function(object) {
      edata(object@metabolic_network)
    },
    where = ns
  )
  methods::setMethod(
    "edata<-",
    "Metabolic_flux_network",
    definition = function(object, value) {
      edata(object@metabolic_network) <- value
      object
    },
    where = ns
  )
}
