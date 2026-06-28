#' @description Df to wb.
#' @title df_to_wb
#'
#' @param df data.frame
#' @param sheet_name str
#'
#' @return wb
#' @export
#'

df_to_wb <- function(df, sheet_name = "data",rowname=T){

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb,sheetName = sheet_name)
  openxlsx::writeData(wb ,sheet = 1, x = df,rowNames = rowname,na.string  = NA)
  wb
}

wb_to_df <- function(wb,sheet_name = 1){
  openxlsx::read.xlsx(wb,sheet = sheet_name)
}

#' @description Edit df in excel.
#' @title edit_df_in_excel
#'
#' @param df data.frame
#' @param rowname Logical; write row names to the worksheet.
#' @param read_line Logical; if \code{TRUE} (default), wait for user input and
#'   return the edited data frame. If \code{FALSE}, open Excel and return
#'   \code{NULL} without reading the file back.
#'
#' @return Edited \code{data.frame}, or \code{NULL} when \code{read_line = FALSE}.
#' @export
#'

edit_df_in_excel <- function(df = data.frame(), rowname = T, read_line = T){
  wb <- df_to_wb(df,rowname=rowname)
  temp.xlsx <- paste0(tempdir(), "/temp_",paste0(sample(letters,5),collapse = ""),".xlsx")
  openxlsx::saveWorkbook(wb , file = temp.xlsx)
  openxlsx::openXL(temp.xlsx)
  if (!read_line) {
    return(NULL)
  }
  readline("press any key to continue")
  wb <- openxlsx::loadWorkbook(file = temp.xlsx)
  df <- wb_to_df(wb)
  df
}

#' Write data to Excel file
#'
#' Writes a data frame or other object to an Excel file, creating the
#' directory if it does not exist.
#'
#' @param data Data frame or object to write to Excel.
#' @param file.dir File path for the output Excel file.
#'
#' @return Invisible NULL.
#' @export
#'

write.xlsx <- function(data,file.dir ){

  dir.create(dirname(file.dir),recursive = T,showWarnings = F)
  openxlsx::write.xlsx(data , file.dir)

}



#' Write list of data frames to Excel
#'
#' Writes a list of data frames (or other objects) to a single Excel
#' workbook, each element as a separate sheet. Sheet names are taken from
#' list names, with invalid characters replaced by underscores.
#'
#' @param df.list List of data frames or objects to write.
#' @param file File path for the output Excel workbook.
#'
#' @return Invisible NULL.
#' @export
#'

xlsx.write.list <- function(df.list ,file){

  wb <- openxlsx::createWorkbook()
  for (i in 1:length(df.list)) {

    sheet.name <- names(df.list)[i]
    sheet.name <- gsub(x = sheet.name,pattern = "[:;/]",replacement = "_")
    openxlsx::addWorksheet(wb,sheetName = sheet.name)
    openxlsx::writeData(wb ,sheet = i, x = df.list[[i]],na.string =NA)


  }
  openxlsx::saveWorkbook(wb,file = file,overwrite = T)

  return(invisible())

}








