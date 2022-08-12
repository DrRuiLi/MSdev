

#' @title show_ms_condition
#' @description show indicated ms_condition record,
#' or return all record
#' @param MSC_id *character* or *integer*, indicate MSC_id or its index
#'
#' @return
#' @export
#'
#' @examples
show_ms_condition <- function(MSC_id = "all") {
  data("ms_condition")
  if (MSC_id == "all") {
    return(as_tibble(ms_condition))
  }
  if (is.numeric(MSC_id)) {
    MSC_id <- ms_condition$MSC_id[MSC_id]
  }
  if (!MSC_id %in% ms_condition$MSC_id) {
    return()
  }

  to_show  <- ms_condition[which(ms_condition$MSC_id == MSC_id), ]
  column <- to_show$column[[1]]
  phaseA <- to_show$phaseA[[1]]
  phaseB <- to_show$phaseB[[1]]
  gradient <- to_show$gradient[[1]]


  ggplot(gradient) +
    geom_line(aes(x = time , y = Contentration_B)) +
    labs(
      title = to_show$name,
      subtitle = paste0(
        "Column : ",
        column$values[1],
        " , ",
        column$values[4],
        " x ",
        column$values[5],
        " mm , ",
        column$values[3],
        " μm\n",
        "Phase A : ",
        paste0(
          str_c(phaseA$Concentration , phaseA$Compound, sep = " "),
          collapse = " + "
        ),
        "\nPhase B : ",
        paste0(
          str_c(phaseB$Concentration , phaseB$Compound, sep = " "),
          collapse = " + "
        )
      )
    ) +
    ylim(c(0, 100)) +
    theme_bw()

}
