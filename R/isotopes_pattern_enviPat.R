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
      rowwise() %>%
      mutate(f = paste0(element, `sum(n)`)) %>%
      pull(f) %>%
      paste0(collapse = "")

  }
  isopat <- data.frame(formula = formulat_list ,
                       isopat[, 1:2])
  return(isopat)

}
