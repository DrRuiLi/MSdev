max_min_normalize<- function(x){

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
