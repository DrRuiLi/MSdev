#' @title str_short
#' @description cut string with length >n and replace char > n with "..."
#' @param x
#' @param n
#'
#' @return
#' @export
#'
#' @examples
str_short <- function(x , n = 10){

  x.length <- nchar(x)
  x.new <- stringr::str_sub(x , 0 , n)%>%
    paste0(ifelse(x.length > n,"...",""))

}



#' @title groupStringFactor
#' @description convert a vector to a factor,
#' levels will be ordered by 1: string, "con" or "wt" will be placed in first;
#' 2: others will be ordered according to0 number
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
groupStringFactor <- function(x){

  x.con <- x[grepl(pattern = "con|wt",x,ignore.case = T)]%>%
    sort()
  x.other <- x[!x%in% x.con]

  ### number
  x.other.num <- str_extract(x.other,"[0-9]+")%>%as.numeric()
  x.new <- c(x.con,x.other[order(x.other.num)])


  ### letter
 # x.new <- c(x.con,sort(x.other))
  factor(x,levels = unique(x.new))
}


vector2str <- function(v){

  paste0("c(\"",paste0(na.omit(v),collapse = "\",\""),"\")")%>%
    cat()
}
