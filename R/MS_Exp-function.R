


.remove_MS_Exp <- function(x , i){
  new.object <- MS_Exp()
  new.object@General <- x@General[-i,]
  new.object@Pre_process <- x@Pre_process[-i]
  new.object@Moblie_phase_A <- x@Moblie_phase_A[-i]
  new.object@Moblie_phase_B <- x@Moblie_phase_B[-i]
  new.object@Chroma_column <- x@Chroma_column[-i]
  new.object@Chroma_gradient <- x@Chroma_gradient[-i]
  new.object@Mass_Spectrum <- x@Mass_Spectrum[-i]
  new.object@Internal_Standard <- x@Internal_Standard[-i]
  new.object

}
.plot_MS_Exp <- function(x){

  to_show <- x@General
  column <- x@Chroma_column%>%as.data.frame()
  phaseA <- x@Moblie_phase_A%>%as.data.frame()%>%
    arrange(Type)
  phaseB <-  x@Moblie_phase_B%>%as.data.frame()%>%
    arrange(Type)
  gradient <- x@Chroma_gradient%>%as.data.frame()

  ggplot(gradient) +
    geom_line(aes(x = time , y = Contentration_B)) +
    labs(
      title = to_show$Name,
      subtitle = paste0(
        "Column : ",
        column$Column_name,
        " , ",
        column$Length,
        " x ",
        column$Diameter,
        " , ",
        column$Paricle_size,
        "\n",
        "Phase A : ",
        paste0(
          stringr::str_c(phaseA$Concentration , phaseA$Compound, sep = " "),
          collapse = " + "
        ),
        "\nPhase B : ",
        paste0(
          stringr::str_c(phaseB$Concentration , phaseB$Compound, sep = " "),
          collapse = " + "
        )
      )
    ) +
    ylim(c(0, 100)) +
    theme_bw()
}
.MS_Exp_to_workbook <- function(x){

  if (length(x) >1) {
    stop("length of MS_Exp > 1")
  }

  MS_workbook <- openxlsx::createWorkbook()
  sapply(slotNames(x), function(y){
    openxlsx::addWorksheet( wb = MS_workbook,sheetName = y)
  } )
  x@General %>%t %>%as.data.frame()%>%
    select(value = 1)%>%
    openxlsx::writeData(wb = MS_workbook,sheet = 1,rowNames = T)
  x@Pre_process %>%as.data.frame()%>%
    openxlsx::writeData(wb = MS_workbook,sheet = 2,rowNames = F)
  x@Moblie_phase_A %>%as.data.frame()%>%
    openxlsx::writeData(wb = MS_workbook,sheet = 3,rowNames = F)
  x@Moblie_phase_B %>%as.data.frame()%>%
    openxlsx::writeData(wb = MS_workbook,sheet = 4,rowNames = F)
  x@Chroma_column%>%unlist%>%as.data.frame()%>%
    select(value = 1)%>%
    openxlsx::writeData(wb = MS_workbook,sheet = 5,rowNames = T)
  x@Chroma_gradient %>%as.data.frame()%>%
    openxlsx::writeData(wb = MS_workbook,sheet = 6,rowNames = F)
  x@Mass_Spectrum%>%unlist%>%as.data.frame()%>%
    select(value = 1)%>%
    openxlsx::writeData(wb = MS_workbook,sheet = 7,rowNames = T)
  x@Internal_Standard%>%as.data.frame()%>%
    openxlsx::writeData(wb = MS_workbook,sheet = 8,rowNames = F)

  MS_workbook

}
.workbook_to_MS_Exp<-function(wb){
  x <- MS_Exp()
  x@General <- openxlsx::read.xlsx(wb , sheet = 1,rowNames = T)%>%t%>%
    as.tibble()
  x@Pre_process <- openxlsx::read.xlsx(wb , sheet = 2,rowNames = F)%>%as.tibble()%>%list()
  x@Moblie_phase_A <- openxlsx::read.xlsx(wb , sheet = 3,rowNames = F)%>%as.tibble()%>%list()
  x@Moblie_phase_B <- openxlsx::read.xlsx(wb , sheet = 4,rowNames = F)%>%as.tibble()%>%list()
  x@Chroma_column <- openxlsx::read.xlsx(wb , sheet = 5,rowNames = T)%>%t%>%
    as.tibble()%>%list()
  x@Chroma_gradient <- openxlsx::read.xlsx(wb , sheet = 6,rowNames = F)%>%as.tibble()%>%list()
  x@Mass_Spectrum <-openxlsx::read.xlsx(wb , sheet = 7,rowNames = T)%>%t%>%
    as.tibble()%>%list()
  x@Internal_Standard <-openxlsx::read.xlsx(wb , sheet = 8,rowNames = F)%>%
    as.data.frame()%>%list()
  x

}
create_MS_Exp_record <-function( copy_from = 1,edit = F){

  MS_Experiment.file <- system.file("data",
                                   "MS_Experiment.rda",
                                   package = "MSdev"
  )
  load(MS_Experiment.file)
  x <- MS_Experiment[copy_from]
  if (copy_from==0) {
    x <- MS_Exp()
  }
  x@General$Creat_time <- as.character(Sys.time())
  if(edit){
    wb <- .MS_Exp_to_workbook(x)
    temp.xlsx <- paste0(tempdir(), "/temp_",paste0(sample(letters,5),collapse = ""),".xlsx")
    openxlsx::saveWorkbook(wb , file = temp.xlsx)
    openxlsx::openXL(temp.xlsx)
    readline("press any key to continue")
    wb <- openxlsx::loadWorkbook(file = temp.xlsx)
    x <- .workbook_to_MS_Exp(wb)
  }
  if (x@General$MSE_id %in% MS_Experiment@General$MSE_id) {
    x@General$MSE_id <-MS_Experiment@General$MSE_id %>%
      str_extract(pattern = "[0-9]+")%>%
      as.numeric()%>%
      setdiff(1:1e5,.)%>%
      min()%>%
      sprintf("%05d",.)%>%
      paste0("MSE",.)
  }
  MS_Experiment <- c(MS_Experiment,x)
  save(MS_Experiment ,file =  MS_Experiment.file)
  return(MS_Experiment)
  }
remove_MS_Exp_record <- function(i){
  MS_Experiment.file <- system.file("data",
                                    "MS_Experiment.rda",
                                    package = "MSdev"
  )
  load(MS_Experiment.file)
  MS_Experiment <- .remove_MS_Exp(MS_Experiment,i)
  save(MS_Experiment ,file =  MS_Experiment.file)
  return(MS_Experiment)

}
edit_MS_Exp_record <- function(i ){
  MS_Experiment.file <- system.file("data",
                                    "MS_Experiment.rda",
                                    package = "MSdev"
  )
  load(MS_Experiment.file)
  x <- MS_Experiment[i]
  wb <- .MS_Exp_to_workbook(x)
  temp.xlsx <- paste0(tempdir(), "/temp_",paste0(sample(letters,5),collapse = ""),".xlsx")
  openxlsx::saveWorkbook(wb , file = temp.xlsx)
  openxlsx::openXL(temp.xlsx)
  readline("press any key to continue")
  wb <- openxlsx::loadWorkbook(file = temp.xlsx)
  x <- .workbook_to_MS_Exp(wb)
  seq <- 1:length(MS_Experiment)
  MS_Experiment <- MS_Experiment[seq[seq <i]]%>%
    c(x)%>%
    c( MS_Experiment[seq[seq >i]])
  save(MS_Experiment ,file =  MS_Experiment.file)
  return(MS_Experiment)
}
show_MS_Exp_record <- function( i = "all"){
  MS_Experiment.file <- system.file("data",
                                    "MS_Experiment.rda",
                                    package = "MSdev"
  )
  load(MS_Experiment.file)

  if (i == "all") {
    return(MS_Experiment)
  }
  if ( is.character(i)) {
    i <- which(MS_Experiment@General$MSE_id == i)
  }
  if (rlang::is_empty(i) | i > length(MS_Experiment)) {
    stop("Record not exist")
  }
  .plot_MS_Exp(MS_Experiment[i])

}
update_MS_Exp_record <- function(){

  MS_Experiment.file <- system.file("data",
                                    "MS_Experiment.rda",
                                    package = "MSdev"
  )
  load(MS_Experiment.file)
  x <- MS_Exp()
  colnames(MS_Experiment@General)
  colnames(x@General)
  to.add <-  setdiff(colnames(x@General),
                     colnames(MS_Experiment@General))
  to.remove <-setdiff(colnames(MS_Experiment@General),
                      colnames(x@General))
  MS_Experiment@General <- MS_Experiment@General %>%
    add_multi_column(to.add)%>%
    select(colnames(x@General))%>%
    as.tibble()

  save(MS_Experiment ,file =  MS_Experiment.file)
  MS_Experiment
}
