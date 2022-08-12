setClass("MS_Exp",
         slots =list(
           "General" = "tbl",
           "Pre_process" = "tbl",
           "Moblie_phase_A" = "tbl",
           "Moblie_phase_b" = "tbl",
           "Chroma_gradient" = "tbl",
           "Mass_Spectrum" = "tbl"
         ),
         prototype = list(
           General = tibble(
             "MSC_id" =  "MSC0001",
             "name" =  "Metabolomics",
             "creat_time" = paste0(Sys.time()),
             "instrument" = "SCIEX TripleTOF 6600",
             "data_aquisition" = "DDA TOP10",
             "link" = "",
             "note" = ""
           ),
           Pre_process = tibble(

           )




         ))

MS_Exp <- function() {
  new("MS_Exp")
  }


setMethod("show" ,"MS_Exp", definition = function(object){

  general_data <- object@General
  show(general_data)

})


object <- MS_Exp()



