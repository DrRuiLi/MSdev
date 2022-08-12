add_multi_column <- function(x, column_name){

  if (length(column_name)==0) {
    message("no column to add")
    return(x)
  }
  df.to.add <- matrix(nrow = nrow(x),
                      ncol = length(column_name))
  colnames(df.to.add) <- column_name
  cbind(x , df.to.add)%>%
    as.tibble()

}
