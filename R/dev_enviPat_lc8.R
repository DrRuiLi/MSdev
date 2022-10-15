formula_calculate_lc8 <- function(Formula1 = "C2H4O1S2P1",Formula2 = "N1H1O-1",sign = 1,Valid_formula = FALSE ){

  to.return <- lc8::my_calculate_formula(Formula1 ,Formula2 ,sign ,Valid_formula )
  ### when formula calculate result to 0, such as "CH - CH = NULL", it will return "NANA"
  to.return[to.return == "NANA"] <- "C0"
  return(to.return)


}




#' @title isotopes_pattern_enviPat
#' @description update of enviPat::isopattern, which do not output the formula of isotopologues
#' this function return isotopologues with abundance more than 1%
#' @param formula
#'
#' @return
#' @export
#'
#' @examples
isotopes_pattern_enviPat <- function(formula) {
  #formula <- "C80[13]C3H33[2]H12"
  data("isotopes",package = "enviPat")
  isopat <-
    enviPat::isopattern(isotopes = isotopes,
                        chemforms = formula,
                        threshold = 1)[[1]]
  iso.matrix <- isopat[, 3:ncol(isopat)]

  ele <- colnames(iso.matrix)
  #ele <- str_replace(ele,"[0-9]+",paste0("[",str_match(ele,"[0-9]+"),"]"))

  formulat_list <- c()
  for (i in 1:nrow(iso.matrix)) {
    formulat_list[i] <-
      data.frame(element = ele ,
                 n = iso.matrix[i, ]) %>%
      mutate(element = elem_table$element[match(element , elem_table$isotope)]) %>%
      group_by(element) %>%
      summarise(sum(n)) %>%
      filter(`sum(n)`!=0)%>%
      rowwise() %>%
      mutate(f = paste0(element, `sum(n)`)) %>%
      pull(f) %>%
      paste0(collapse = "")

  }
  formula_raw <- formula
  select_elemet <- function(x ){
    ele_table <- lc8::my_break_formula(x)%>%as.data.frame()%>%
      filter(count >0)%>%
      mutate(f = paste0(elem,count))%>%
      pull(f)%>%
      paste0(collapse = "")
    ele_table
  }
  isopat <- data.frame(formula = formulat_list ,
                       isopat[, 1:2])%>%
    rowwise()%>%
    mutate(isotope_element = formula_calculate_lc8(formula,formula_raw , -1),
           isotope_element = select_elemet(isotope_element))


  return(isopat)

}
