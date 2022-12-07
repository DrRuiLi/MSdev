#' @title df_to_wb
#'
#' @param df
#' @param sheet_name
#'
#' @return
#' @export
#'
#' @examples
df_to_wb <- function(df, sheet_name = "data"){

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb,sheetName = sheet_name)
  openxlsx::writeData(wb ,sheet = 1, x = df)
  wb
}

wb_to_df <- function(wb,sheet_name = 1){
  openxlsx::read.xlsx(wb,sheet = sheet_name)
}

#' @title edit_df_in_excel
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
edit_df_in_excel <- function(df){
  wb <- df_to_wb(df)
  temp.xlsx <- paste0(tempdir(), "/temp_",paste0(sample(letters,5),collapse = ""),".xlsx")
  openxlsx::saveWorkbook(wb , file = temp.xlsx)
  openxlsx::openXL(temp.xlsx)
  readline("press any key to continue")
  wb <- openxlsx::loadWorkbook(file = temp.xlsx)
  df <- wb_to_df(wb)
df
}

#' write.xlsx
#'
#' @param data
#' @param file.dir
#'
#' @return
#' @export
#'
#' @examples
write.xlsx <- function(data,file.dir ){

  dir.create(dirname(file.dir),recursive = T)
  openxlsx::write.xlsx(data , file.dir)

}
