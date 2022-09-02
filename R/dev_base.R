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
