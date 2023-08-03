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
#' @description write xlsx, dir will be created if not existed
#' @param data
#' @param file.dir
#'
#' @return
#' @export
#'
#' @examples
write.xlsx <- function(data,file.dir ){

  dir.create(dirname(file.dir),recursive = T,showWarnings = F)
  openxlsx::write.xlsx(data , file.dir)

}



#' xlsx.write.list
#'
#' Write a list of data.frame (or other object could be write into wroksheet), every term will be writen to a single sheet
#'
#' @param df.list
#' @param file
#'
#' @return
#' @export
#'
#' @examples
xlsx.write.list <- function(df.list ,file){

  wb <- openxlsx::createWorkbook()
  for (i in 1:length(df.list)) {

    sheet.name <- names(df.list)[i]
    sheet.name <- gsub(x = sheet.name,pattern = "[:;/]",replacement = "_")
    openxlsx::addWorksheet(wb,sheetName = sheet.name)
    openxlsx::writeData(wb ,sheet = i, x = df.list[[i]])


  }
  openxlsx::saveWorkbook(wb,file = file,overwrite = T)

  return(invisible())

}








