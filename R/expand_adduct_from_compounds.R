#' @title  Expand compound to adduct
#' @description Using formula, exact.mass and polarity, expand adduct by enviPat::adduct
#'
#' @param compound.record
#' @param polarity
#'
#' @return
#' @export
#'
#' @examples
expand_adduct_from_compounds <- function(compound.record, polarity) {
  adduct.rule <- adducts %>%
    dplyr::filter(Ion_mode == polarity) %>%
    dplyr::mutate(
      Formula_add = case_when(Formula_add == "FALSE" ~  "C0" , T ~ Formula_add),
      Formula_ded = case_when(Formula_ded == "FALSE" ~  "C0" , T ~ Formula_ded)
    )
  fff <- function(x) {
    to.return <- list()
    to.return[["compound"]] <- x
    x.formula <- x["formula"]
    x.mass <- check_chemform(isotopes,x.formula)$monoisotopic_mass
    to.return[["adduct.candidate"]] <- adduct.rule %>%
      rowwise() %>%
      mutate(
        formula =  multiform(formula1 = x.formula , fact = Multi),
        formula = mergeform(formula, Formula_add),
        formula = subform(formula, Formula_ded),
        adduct = Name,
        charge = Charge,
        multi = Multi,
        ion_mode = Ion_mode,
        form = "adduct",
        exact.mz = x.mass * Multi / abs(Charge) + Mass
      ) %>%
      select("formula"  ,
             "adduct",
             "charge",
             "multi",
             "ion_mode",
             "form",
             "exact.mz")

    return(to.return)

  }
  MS.network <- apply(compound.record, 1,fff  )

return(MS.network)

}
