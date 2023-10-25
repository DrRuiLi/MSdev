t.test_dev <- function(x,y){

  try.catch <-try(p.value <- t.test(x,y)$p.value)
  if (grepl("Error " , try.catch)) {
    return(1)
  }
  return(p.value)
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
    dplyr::mutate(mz.center = median(mz),
                  mz.diff = abs(mz.center - mz),
                  mz.ppm = mz.diff/mz*1e6,
                  mz.width = max(mz)-min(mz),
                  mz.width.ppm = mz.width/mz)%>%
    dplyr::ungroup()%>%
    dplyr::arrange(raw.order)%>%
    dplyr::select(-raw.order)
  return(x.table)


}


cluster_ion <- function(mz,rt , ppm.thresh =10, rt.tol = 15){

  .cluster_rt <- function(rt,rt.tol){
    if (length(rt)==1) {
      return(1)
    }
    dist(rt)%>%
      hclust()%>%
      cutree(h = rt.tol)

  }

  ion.df <- data.frame(mz,rt)%>%
    dplyr::mutate(groupMz(mz,ppm.thresh = ppm.thresh))%>%
    dplyr::group_by(mz.group)%>%
    dplyr::mutate(rt.cluster = .cluster_rt(rt,rt.tol))%>%
    dplyr::group_by(mz.group,rt.cluster)%>%
    dplyr::mutate(ion.cluster = paste0("Ion_cluster",sprintf("%08d",cur_group_id())))%>%
    dplyr::ungroup()%>%
    dplyr::pull(ion.cluster)
  return(ion.df)




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



#' @title normalize_max_min
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
normalize_max_min<- function(x){

  f <- function(z){
    (z-min(z,na.rm = T))/(max(z,na.rm = T)-min(z,na.rm = T))
  }
  if (is.vector(x)) {
    return(f(x))
  }
  if (is.matrix(x)|is.data.frame(x)) {
    apply(x, 1, f)%>%t%>%return()
  }

}

mz.range.ppm <- function(mz = 200, ppm = 5){

  mz.range <- c(mz - mz * ppm /1e6,mz+mz*ppm/1e6)
  return(mz.range)
}


expand_range <- function(x= c(5,10),add = 0,multi = 0){

  x <- c(x[1]-diff(x)*multi-add,
         x[2]+diff(x)*multi+add)
  return(x)
}




gaussian_functioin <- function(x , a =1,b = 0,c = 0.5){

  a * exp(-(x-b)^2/2/c^2)
}

