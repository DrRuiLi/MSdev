formula_calculate_lc8 <- function(Formula1 = "C2H4O1S2P1",Formula2 = "N1H1O-1",sign = 1,Valid_formula = FALSE ){

  to.return <- lc8::my_calculate_formula(Formula1 ,Formula2 ,sign ,Valid_formula )
  to.return[to.return == "NANA"] <- "C0"
  return(to.return)


}
