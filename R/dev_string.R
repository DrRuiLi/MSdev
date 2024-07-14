#' @title str_short
#' @description cut string with length >n and replace char > n with "..."
#' @param x str
#' @param n char number
#'
#' @return str
#' @export
#'

str_short <- function(x , n = 10){

  x.length <- nchar(x)
  x.new <- stringr::str_sub(x , 0 , n)%>%
    paste0(ifelse(x.length > n,"...",""))
  return(x.new)

}



#' @title groupStringFactor
#' @description convert a vector to a factor,
#' levels will be ordered by 1: string, "con" or "wt" will be placed in first;
#' 2: others will be ordered according to0 number
#'
#' @param x str
#'
#' @return str with factor
#' @export
#'

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


#' @title vector2str
#' @param v vector
#'
#' @return str
#' @export
#'

vector2str <- function(v,verbose = F){

  x <- paste0("c(\"",paste0(na.omit(v),collapse = "\",\""),"\")")
  if (verbose)
   cat(x)
    return(x)
}



num2str <- function(x,n.digit = NA){

  if (is.na(n.digit)) n.digit <- max(nchar(x))
  sp.exp <- paste0("%0",n.digit,"d")
  sprintf(sp.exp,x)

}

str_extract_num <- function(x){

  stringr::str_extract(x,"[:digit:]+")%>%
    as.numeric()
}



str_add <- function(x , n = 1, add.type = "numeric"){

  if (add.type == "numeric") {

      x.char <- str_extract(x,pattern = "[^0-9]+")
      x.num <- str_extract(x , "[0-9]+")
      x <- paste0(x.char,num2str(as.numeric(x.num)+n,
                                 n.digit = max(nchar(x.num))))

    }
return(x)

}


str_formate <- function(str){

  str
  gsub(x = str , pattern =  "[^A-z0-9]" ,
       replacement = "_",perl = T)

}


str_format_len <- function(str,to = "max"){

  nc <- nchar(str)
  if (to=="max") {
    nc.max <- max(nc)
    to.add <- sapply(nc.max-nc,
                     function(x){
                       paste0(rep("  ",x),collapse = "")
                     })
    str1 <- paste0(str,to.add)
  }

  return(str1)
}



str_search_files <- function(str,file_path = paste0(getwd(),"/R")){

  if (any(file.info(file_path)$isdir) ) {
    file_path <- dir(file_path,full.names = T,recursive = T)
  }

  for (i in seq_along(file_path)) {

    txt <- readr::read_lines(file_path[i])
    txt.exist <-  grepl(str,txt)
    if (any(txt.exist)) {

      #message(i)
      message(crayon::red(file_path[i]))
      message(crayon::red("Line ", which(txt.exist), ":",
                          yellow(txt[txt.exist]),
                          collapse = "\n") )

    }
  }


}
