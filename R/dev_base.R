t.test_dev <- function(x,y){

  try.catch <-try(p.value <- t.test(x,y)$p.value)
  if (grepl("Error " , try.catch)) {
    return(1)
  }
  return(p.value)
}



#' @title date_suffix
#' @description generate a character of current date
#' @return
#' @export
#'
#' @examples
date_suffix <- function(){

  paste0("_",gsub(Sys.Date(),pattern = "-",replacement = "_"))

}



#' @title matrixSub
#' @description expand 2 vector to matrix and "`-`"
#' @param v1 vector
#' @param v2 vector
#'
#' @return matrix
#'
#' @export
#'
#' @examples
#' a <- 3:8
#' b <- 1:2
#' matrixSub(a,b)
matrixSub <- function(v1,v2){

  m1 <- matrix(rep(v1,length(v2)),nrow = length(v1),byrow = F)
  m1
  m2 <- matrix(rep(v2,length(v1)),ncol = length(v2),byrow = T)
  m2
  x <- m1-m2
  rownames(x) <- names(v1)
  colnames(x) <- names(v2)
  return(x)
}



#' @title load_as_var
#' @description load Rdata file and return, only one variable in `file_to_load`
#' @param file_to_load
#'
#' @return
#' @export
#'
#' @examples
load_as_var <- function(file_to_load){

  var <- load(file_to_load)
  if (length(var)!=1) {
    stop("Too many variabls in ",file_to_load)
  }
  eval(str2expression(var))

}


