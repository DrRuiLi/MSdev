

make_MS_network_rule <- function() {
  data("isotopes", package = "enviPat")
  data("adducts", package = "enviPat")
  data("empirical_rules_from_NetID")

  adduct_MSdev <- adducts %>%
    rowwise()%>%
    mutate(Formula_add = case_when(Formula_add == "FALSE"~ "C0",
                                   T ~ Formula_add),
           Formula_ded = case_when(Formula_ded == "FALSE"~ "C0",
                                   T ~ Formula_ded),
           Formula_diff = formula_calculate_lc8(Formula_add,Formula_ded, -1))%>%
    select(-Formula_add,-Formula_ded)



}
