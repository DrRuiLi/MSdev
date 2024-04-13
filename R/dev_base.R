
#' @title date_suffix
#' @description generate a character of current date
#' @return str
#' @export
#'

date_suffix <- function(){

  paste0("_",gsub(Sys.Date(),pattern = "-",replacement = "_"))

}






#' @title load_as_var
#' @description load Rdata file and return, only one variable in `file_to_load`
#' @param file_to_load file path
#'
#' @return data
#' @export
#'

load_as_var <- function(file_to_load){

  var <- load(file_to_load)
  if (length(var)!=1) {
    stop("Too many variabls in ",file_to_load)
  }
  eval(str2expression(var))

}


#' @title colorMix
#' @description mix color by RGB and weighted by alpha
#' @param ... colors
#'
#' @return color
#' @export
#'

colorMix <- function(...){
  col.list <- list(...)
  col.df <- lapply(col.list,function(x){
    if (is.na(x)) {
      x <- "#FFFFFF00"

    }
    data.frame(t(col2rgb(x,alpha = T)))

  })%>%data.table::rbindlist()%>%
    dplyr::mutate(r = red*alpha/255/255,
                  g = green *alpha/255/255,
                  b = blue*alpha/255/255,
                  a = alpha/255)%>%
    summarise_all(sum)
  rgb(red = col.df$r,
      green = col.df$g,
      blue = col.df$b,
      maxColorValue = col.df$a)

}




#' @title open_dir
#' open dir by Windows
#'
#' @param dir a path of dir or file
#'
#' @return null
#' @export
#'

open_dir <- function(dir = getwd()){

  if (!file.info(dir)$isdir) {
    dir <- dirname(dir)

  }
  system(sprintf("open %s", shQuote(dir)))
  return(dir)

}

open_dir_ActivedFilePath <- function(){
  path <- dirname(rstudioapi::getSourceEditorContext()$path)
  open_dir(path)
}

open_R_libPath <- function(){

  open_dir(.libPaths()[1])


}



open_file <- function(dir){
  system(sprintf("open %s", shQuote(dir)))
  return(dir)

}

open_script <- function(){

  if (!file.exists("./Script/Test_Script.R")) {
    file.create("./Script/Test_Script.R")
  }
  rstudioapi::documentOpen("./Script/Test_Script.R")
}



#' split_df
#' random split data.frame
#'
#' @param df data.frame
#' @param n number
#'
#' @return list
#' @export
#'

split_df <- function(df,n = 2){


  split(df,f = sample(1:n,nrow(df),replace = T))


}



seq_unique <- function(x){

  seq_along(na.omit(unique(x)))

}



get_file_formate <- function(file){


  str_extract(pattern  = "\\.[^\\.]*$",string = basename(file))%>%
    sub("\\.",x = .,"")

}




