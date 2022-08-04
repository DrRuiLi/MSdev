edit_ms_condition <- function(MSC_id = "MSC0001"){

  ms_condition_file <- system.file("data",
                                   "ms_condition.rda",
                                   package = "MSdev"
  )
  load(ms_condition_file)
  if (is.numeric(MSC_id)) {
    MSC_id <- ms_condition$MSC_id[MSC_id]
  }
  ms_condition_to_edit <- ms_condition[which(ms_condition$MSC_id== MSC_id),]


  ms_condition_to_edit$column[[1]] -> ms_column
  ms_condition_to_edit$phaseA[[1]] -> ms_phase_A
  ms_condition_to_edit$phaseB [[1]]-> ms_phase_B
  ms_condition_to_edit$gradient [[1]]-> ms_gradient

  ms_condition_table <- ms_condition_to_edit[-which(colnames(ms_condition_to_edit) %in% c("column","phaseA","phaseB","gradient"))]

  ms_condition_table <- data.frame(ms_condition_table,row.names = "values")%>%
      t%>%as.data.frame()%>%rownames_to_column("params")

    ### edit in excel
    temp.xlsx <- paste0(tempdir(), "/temp_",paste0(sample(letters,5),collapse = ""),".xlsx")
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




    ### load new data
    {
      ms_condition_table <-
        openxlsx::read.xlsx(ms_workbook, sheet = "Exp_Condition")
      ms_column <- openxlsx::read.xlsx(ms_workbook, sheet = "Column")
      ms_phase_A <- openxlsx::read.xlsx(ms_workbook, sheet = "phase A")
      ms_phase_B <- openxlsx::read.xlsx(ms_workbook, sheet = "phase B")
      ms_gradient <-
        openxlsx::read.xlsx(ms_workbook, sheet = "Gradient")


      ms_condition_to_edit <- DataFrame(MSC_id = NA)
      for (i in 1:nrow(ms_condition_table)) {
        exprss <-
          paste0(
            "ms_condition_to_edit$",
            ms_condition_table$params[i],
            " <- \"",
            ms_condition_table$values[i],
            "\""
          )
        eval(parse(text = exprss))
      }
      ms_condition_to_edit$column <- ms_column %>% list()
      ms_condition_to_edit$phaseA <- ms_phase_A %>% list()
      ms_condition_to_edit$phaseB <- ms_phase_B %>% list()
      ms_condition_to_edit$gradient <- ms_gradient %>% list()
      }


if (ms_condition_to_edit$MSC_id == "NA") {
  ms_condition <- ms_condition[-which(ms_condition$MSC_id== MSC_id),]
}else{

  ms_condition_to_edit -> ms_condition[which(ms_condition$MSC_id== MSC_id),]
}



  save(ms_condition, file = ms_condition_file)

}
