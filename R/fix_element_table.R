fix_element_table <- function(){

  data("isotopes",package = "enviPat")
  isotopes_new <- isotopes%>%
    mutate(a =  str_replace(isotope,"[0-9]+",paste0("[",str_match(isotope,"[0-9]+"),"]")))%>%
    filter(!element == a)%>%
    select(element = a , isotope,mass,abundance,ratioC)
  isotopes <- rbind(isotopes,isotopes_new)
  use_data(isotopes,overwrite = T)

  data("elem_table",package = "lc8")
  elem_table_new <- elem_table%>%
    mutate(a =  str_replace(isotope,"[0-9]+",paste0("[",str_match(isotope,"[0-9]+"),"]")))%>%
    filter(!element == a)%>%
    select(element = a , isotope,mass,abundance,ratioC)
  elem_table <- rbind(elem_table,elem_table_new)
  use_data(elem_table,overwrite = T)

}
