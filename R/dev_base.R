t.test_dev <- function(x,y){

  try.catch <-try(p.value <- t.test(x,y)$p.value)
  if (grepl("Error " , try.catch)) {
    return(1)
  }
  return(p.value)
}



date_suffix <- function(){

  gsub(Sys.Date(),pattern = "-",replacement = "_")

}
