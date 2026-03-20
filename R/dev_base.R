
#' @title date_suffix
#' @description generate a character of current date
#' @return str
#' @export
#'

date_suffix <- function(){

  paste0("_",gsub(Sys.Date(),pattern = "-",replacement = "_"))

}

between.range <- function(x, r){

  rm <- matrix(r,ncol =2,byrow= F)
  if (nrow(rm)==1|nrow(rm)==length(x)) {
    between(x,rm[,1],rm[,2])
  }else{
    stop( "size of x and r not match")
  }
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


#' @export
message_with_time <- function(...){

  message(format(Sys.time(),digits = 0)," ",...)
}


#' @title open_dir
#' open dir by Windows
#'
#' @param dir a path of dir or file
#'
#' @return null
#' @export
#'
setGeneric("open_dir",
           def = function(x = getwd()){

             if (!file.info(x)$isdir) {
               x <- dirname(x)

             }
             shell (sprintf("open %s", shQuote(x)))
             return(x)

           })


setMethod("open_dir",signature = "MSdev",
          definition = function(x){
            open_dir(x@projectInfo$MSdevFile)
          })

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

open_MSIP <- function(){
  msip.file <- "c:/Users/91879/OneDrive/Code/R/Projecct/2024.01.11.MSIP/Workflow.2025.07.09.R"
  rstudioapi::documentOpen(msip.file)
  return(msip.file)
}

open_PAVE <- function(){

  pave.file <- get_dir_expand_from_onedrive("Documents/YLF_Lab/Project/2025.10.10.PAVE/code/PAVE_data_Analysis.R")
  rstudioapi::documentOpen(pave.file)
  setRStudioDir()
  return(pave.file)
}

open_Project_dir <- function(){

  proj.file <- get_dir_expand_from_onedrive("Documents/YLF_Lab/Project/")
  setRStudioDir(proj.file)
  return(proj.file)
}

open_MSIP_dir <- function(){

  open_dir("C:/Users/91879/OneDrive/Code/R/data/MSIP_data")

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



object.size.mb <- function(x){

  format(object.size(x),"MB")
}

size_of <- function(x){
  object.size.mb(x)
}

unlist_to_df <- function(x,name_to = "name",value_to = "value"){
  x.df <- data.frame(
    rep(names(x),times =
          lengths(x)),
    unlist(x)
  )
  colnames(x.df) <- c(name_to,value_to)
  x.df
}


make_vector <- function(x = NA ,name = NULL){

  if (length(x)==length(name)) {
    names(x)<-name
  }
  if (length(x)==1&length(name)) {
    x <- rep(x,length(name))
    names(x)<-name
  }

  return(x)

}


match_path <- function(p1,p2){

  p1 <- normalizePath(p1,mustWork = F)
  p2 <- normalizePath(p2,mustWork = F)
  match(p1,p2)


}


which.na <- function(x){
  which(is.na(x))
}


get_matrix_value_fill_with_NA <- function(mat,
                                          rownames_vec = rownames(mat),
                                          colnames_vec = colnames(mat),
                                          drop = T) {
  # Initialize a matrix to store the results
  result_matrix <- matrix(NA, nrow = length(rownames_vec), ncol = length(colnames_vec))

  # Set row and column names of the result matrix
  rownames(result_matrix) <- rownames_vec
  colnames(result_matrix) <- colnames_vec

  # Loop through each combination of row and column names
  for (i in seq_along(rownames_vec)) {
    for (j in seq_along(colnames_vec)) {
      rowname <- rownames_vec[i]
      colname <- colnames_vec[j]

      # Check if both rowname and colname exist in the matrix
      if (rowname %in% rownames(mat) && colname %in% colnames(mat)) {
        result_matrix[i, j] <- mat[rowname, colname]
      } else {
        result_matrix[i, j] <- NA
      }
    }
  }

  if (drop&length(result_matrix)==1) {
    result_matrix <- as.vector(result_matrix)
  }

  return(result_matrix)
}


mean_matrix <- function(mat1,mat2){

  rn <- unique(c(rownames(mat1),rownames(mat2)))%>%
    groupStringFactor()%>%levels()
  cn <- unique(c(colnames(mat1),colnames(mat2)))%>%
    groupStringFactor()%>%levels()
  m1 <- get_matrix_value_fill_with_NA(mat1,rn,cn)
  m2 <- get_matrix_value_fill_with_NA(mat2,rn,cn)
  combined_array <- array(c(m1, m2), dim = c(nrow(m1), ncol(m2), 2))
  mean_matrix <- apply(combined_array, c(1, 2),
                       function(x) mean(x, na.rm = TRUE))
  dimnames(mean_matrix) <- dimnames(m1)
  return(mean_matrix)
}

sum_matrix <- function(mat1,mat2){

  rn <- unique(c(rownames(mat1),rownames(mat2)))%>%
    groupStringFactor()%>%levels()
  cn <- unique(c(colnames(mat1),colnames(mat2)))%>%
    groupStringFactor()%>%levels()
  m1 <- get_matrix_value_fill_with_NA(mat1,rn,cn)
  m2 <- get_matrix_value_fill_with_NA(mat2,rn,cn)
  combined_array <- array(c(m1, m2), dim = c(nrow(m1), ncol(m2), 2))
  sum_matrix <- apply(combined_array, c(1, 2),
                       function(x) sum(x, na.rm = TRUE))
  dimnames(sum_matrix) <- dimnames(m1)
  return(sum_matrix)
}


setMethod("isEmpty","NULL",definition = function(x){
  return(T)
})



try_until_success <- function(func, ...) {
  success <- FALSE
  result <- NULL

  while (!success) {
    try({
      result <- func(...)
      success <- TRUE
    }, silent = TRUE)
  }

  return(result)
}



rm_except <- function(..., env = parent.frame()) {
  vars <- sapply(substitute(list(...))[-1], deparse)
  all_vars <- ls(envir = env)
  rm_vars <- setdiff(all_vars, vars)
  if (length(rm_vars) > 0) {
    rm(list = rm_vars, envir = env)
  }
}


show_col <- scales::show_col


paste0_without_na <- function(..., sep = "") {
  args <- lapply(list(...), function(x) ifelse(is.na(x), "", x))
  do.call(paste, c(args, sep = sep))
}


get_progress_bar <- function(total_iterations = 100){

  pbar <- progress::progress_bar$new(
    format = "[:bar] :percent in :elapsed, current index: :current",
    total = total_iterations,
    clear = FALSE
  )

  return(pbar)


}


get_dir_expand_from_onedrive <- function(d = "."){

  user.dir <- normalizePath(Sys.getenv("USERPROFILE"),winslash = "/")
  normalizePath(paste0(user.dir,"/Onedrive/",d),winslash = "/")

}


calculate_pulse <- function(x = 1000, mt = F, qr = T, guru = 0.028,rc = 1.02*0.994){

  r.base <- x * 0.024
  r.qr <- ifelse(qr,x * 0.02,0)
  r.food <- ifelse(mt, x * 0.05,0)
  r.food <- ifelse(r.food > 100, 100, r.food)
  r.guru <- x * guru
  r  <- r.base + r.qr + r.food + r.guru

  x  - r * rc

}
