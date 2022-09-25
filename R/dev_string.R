fixStringLength <- function(x , n = 10){

  x.length <- nchar(x)
  x.new <- str_sub(x , 0 , n)%>%
    paste0(ifelse(x.length > n,"...",""))

}



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
