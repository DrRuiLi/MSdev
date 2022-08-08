#' @title  Expand compound to adduct
#' @description Using formula, exact.mass and ion_mode, expand adduct by enviPat::adduct
#'
#' @param compound.record
#' @param ion_mode
#'
#' @return
#' @export
#'
#' @examples
expand_adduct_from_compounds <- function(compound.record, ion_mode) {
  data("adduct_MSdev")
  adduct.rule <- adduct_MSdev %>%
    dplyr::filter(Ion_mode == ion_mode)

  fff <- function(x) {
    #x <- compound.record[1,]
    to.return <- list()
    to.return[["compound"]] <- x
    x.formula <- x["formula"]
    x.mass <- check_chemform(isotopes,x.formula)$monoisotopic_mass
    x.isopat <- isotopes_pattern_enviPat(x.formula)%>%
      mutate(form = paste0(isotope_element , "M"))

    x_adduct <-data.frame()
    for (i in 1:nrow(x.isopat)) {

      m_formula <- x.isopat$formula[i]
      m_mass <- x.isopat$m.z[i]
      m_form <- x.isopat$form[i]
      m_adduct <- adduct.rule %>%
        rowwise() %>%
        mutate(
          formula =  multiform(formula1 = m_formula , fact = Multi),
          formula = formula_calculate_lc8(formula , Formula_diff),
          adduct = Name,
          adduct = sub(pattern = "M",replacement = m_form , x = adduct),
          charge = Charge,
          multi = Multi,
          ion_mode = Ion_mode,
          form = "adduct",
          exact.mz = m_mass * Multi / abs(Charge) + Mass
        ) %>%
        select("formula"  ,
               "adduct",
               "charge",
               "multi",
               "ion_mode",
               "form",
               "exact.mz")
      m_adduct
      x_adduct <- rbind(x_adduct , m_adduct)

    }



    to.return[["adduct.candidate"]] <- x_adduct
    return(to.return)

  }
  MS.network <- apply(compound.record, 1,fff  )

return(MS.network)

}
