#' @title create_ms_condition
#'
#' @description create a set of ms param, call an excel windows to edit and then record in inner data
#' @return
#' @export
#'
#' @examples
create_ms_condition <- function() {
  ### template
  {
    ms_condition_table <- data.frame(
      params = c(
        "MSC_id",
        "name",
        "creat_time",
        "instrument",
        "data_aquisition"
      ),
      values = c(
        "MSC0001",
        "Metabolomics",
        paste0(Sys.time()),
        "SCIEX TripleTOF 6600",
        "DDA TOP10"
      )
    )
    ms_column <-
      data.frame(
        param = c(
          "Column",
          "Manufacturer",
          "Paricle_size",
          "Length",
          "Diameter"
        ),
        values = c("Kinetex C18", "Phenomenex", "1.7", "50", "2.1")
      )


    ms_phase_A <- data.frame(
      Compound = c("H2O", "Formic acid"),
      Type = c("solvent", "solute"),
      Concentration = c(100, 0.001)
    )

    ms_phase_B <- data.frame(
      Compound = c("ACN", "Formic acid"),
      Type = c("solvent", "solute"),
      Concentration = c(100, 0.001)
    )

    ms_gradient <- data.frame(
      time = c(0, 60, 100, 500),
      Contentraion_B = c(0, 50, 100, 100)
    )
  }
  ### manipulate in Excel
  {
    temp.xlsx <- paste0(tempdir(), "/temp.xlsx")
    openxlsx::write.xlsx(
      list(
        "Exp_Condition" = ms_condition_table,
        "Column" = ms_column,
        "phase A" = ms_phase_A,
        "phase B" = ms_phase_B,
        "Gradient" = ms_gradient
      ),
      file = temp.xlsx
    )
    openxlsx::openXL(temp.xlsx)
    readline("Edit MS condition")
    ms_workbook <- openxlsx::loadWorkbook(file = temp.xlsx)
  }

  ### load new data
  {
    ms_condition_table <-
      openxlsx::read.xlsx(ms_workbook, sheet = "Exp_Condition")
    ms_column <- openxlsx::read.xlsx(ms_workbook, sheet = "Column")
    ms_phase_A <- openxlsx::read.xlsx(ms_workbook, sheet = "phase A")
    ms_phase_B <- openxlsx::read.xlsx(ms_workbook, sheet = "phase B")
    ms_gradient <-
      openxlsx::read.xlsx(ms_workbook, sheet = "Gradient")


    ms_condition <- DataFrame(MSC_id = NA)
    for (i in 1:nrow(ms_condition_table)) {
      exprss <-
        paste0(
          "ms_condition$",
          ms_condition_table$params[i],
          " <- \"",
          ms_condition_table$values[i],
          "\""
        )
      eval(parse(text = exprss))
    }
    ms_condition$column <- ms_column %>% list()
    ms_condition$phaseA <- ms_phase_A %>% list()
    ms_condition$phaseB <- ms_phase_B %>% list()
    ms_condition$gradient <- ms_gradient %>% list()
  }

  ### add to database
  {
    to.add <- ms_condition
    ms_condition_file <- system.file("data",
      "ms_condition.rda",
      package = "MSdev"
    )
    load(ms_condition_file)
    to.add$MSC_id <- paste0(
      "MSC",
      sprintf(
        "%04d",
        1 + max(
          grep(ms_condition$MSC_id, pattern = "[0-9]")
        )
      )
    )
    ms_condition <- rbind(
      ms_condition,
      to.add
    )
    save(ms_condition, file = ms_condition_file)
  }
}
