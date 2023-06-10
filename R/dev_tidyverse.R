add_multi_column <- function(x, column_name){

  if (length(column_name)==0) {
    message("no column to add")
    return(x)
  }
  df.to.add <- matrix(nrow = nrow(x),
                      ncol = length(column_name))
  colnames(df.to.add) <- column_name
  cbind(x , df.to.add)%>%
    as_tibble()

}





#' @title list2df
#'
#' @description convert list to data.frame. In each list,
#' sub-list will be replaced with "large list" and missing value will be fill by NA.
#' See `data.table::rbindlist`
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
list2df <- function(x){

  list.vars <- sapply(x, names)%>%
    unlist()%>%
    unique()

  lapply(x , function(z){
    z <- z[list.vars]
    z.is.list <- sapply(z, is.list)
    z[z.is.list] <- "large list"
    z.is.null <- sapply(z , is.null)
    z[z.is.null] <- NA
    z
  })%>%
    data.table::rbindlist(fill = T)%>%
    as.data.frame()->x.df
  ### empty column with name NA will lead to error, unknown reason
  x.df <- x.df[,which(!is.na(colnames(x.df)))]
  x.df

}







dplyr_copy_row <- function(df,nrow.idx,n){

  df.to.add <- df[nrow.idx,] %>%
    add_multi_column(1:n)%>%
    tidyr::pivot_longer(as.character(1:n))%>%
    dplyr::select(-name,-value)
  rbind(df,df.to.add)


}
