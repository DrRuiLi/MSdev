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





groupMz <- function(x,ppm.thresh = 10){

  x.table <- data.frame( mz = x)%>%
    dplyr::mutate(raw.order = 1:length(x))%>%
    dplyr::arrange(mz)%>%
    dplyr::mutate(mz.diff = c(diff(mz),0),
                  mz.ppm = mz.diff/mz*1e6,
                  mz.group = "")
  i <- 1
  i.group <- 1
  this.group.idx <- c()
  while(i <= nrow(x.table)){
    x.table$mz.ppm[i]
    if (x.table$mz.ppm[i] <ppm.thresh) {
      this.group.idx <- c(this.group.idx,i)
      x.table$mz.group[this.group.idx] <- paste0("ion_group",sprintf("%06d",i.group))
    }else{
      this.group.idx <-  c(this.group.idx,i)
      x.table$mz.group[this.group.idx] <- paste0("ion_group",sprintf("%06d",i.group))
      i.group <- i.group+1
      i <- i+1
      this.group.idx <- c()
      next
    }

    i <- i+1
    next

  }

  x.table <-x.table %>%
    dplyr::group_by(mz.group)%>%
    dplyr::mutate(mz.width = max(mz)-min(mz),
                  mz.width.ppm = mz.width/mz)%>%
    dplyr::ungroup()%>%
    dplyr::arrange(raw.order)
  return(x.table)


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


#' @title colorMix
#' @description mix color by RGB and weighted by alpha
#' @param ... colors
#'
#' @return
#' @export
#'
#' @examples
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
#' @return
#' @export
#'
#' @examples
open_dir <- function(dir = getwd()){

  system(sprintf("open %s", shQuote(dir)))

}




cor.mtest <- function(x){

 index.map <- expand.grid(1:ncol(x),1:ncol(x))

 cor.test.x <- function(y,x){

   p <- cor.test(x[,index.map[y,1]],x[,index.map[y,2]])$p.value
  r <-cor(x[,index.map[y,1]],x[,index.map[y,2]])
  return(list(p=p,r=r))

 }

 cor.df<- BiocParallel::bplapply(1:nrow(index.map),cor.test.x,x)%>%
   data.table::rbindlist()

  cor.r.m <- cor.p.m <-
     matrix(ncol = ncol(x),nrow = ncol(x))%>%
    `rownames<-`(colnames(x))%>%
    `colnames<-`(colnames(x))

  cor.r.m[1:25] <- cor.df$r
  cor.p.m[1:25] <- cor.df$p
  return(list(p = cor.p.m,
              r = cor.r.m))



}




#' @title median_part
#' @description similar to head() and tail(), return the median of a vector
#'
#' @param x vector
#' @param n number of element
#'
#' @return
#' @export
#'
#' @examples
median_part <- function(x, n = 10){

  x.indice <- seq_along(x)
  x.median <-x[floor((median(x.indice)-n/2) : (median(x.indice)+n/2))]
  if (length(x.median)>n) {
    x.median <- x.median[1:n]
  }
  x.median
}

