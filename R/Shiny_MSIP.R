#' Title
#'
#' @param object MSIP obj
#'
#' @return NULL
#' @export

MSIP_shiny_start <- function(object){


  ### Start Shiny APP
  {
    shinyApp(ui = MSIP_shiny_ui(),
             server = MSIP_shiny_server(object),
             options = list(host = "0.0.0.0",
                            #port = 6548,
                            launch.browser = T))
    }

}


#' MSIP_shiny_Acq
#'
#' @param object MSIP obj
#'
#' @return NULL
#' @export
MSIP_shiny_Acq <- function(object){

  ### Start Shiny APP
  {
    shinyApp(ui = MSIP_shiny_Acq_ui(),
             server = MSIP_shiny_Acq_server(object),
             options = list(host = "0.0.0.0",
                            launch.browser = T,
                            port = 6547))%>%
      runApp()->object


  }




}



test_fun <- function(input){

  print(names(input))
  print("\n")
  return(0)
}

MFN_manul_Shiny <- function(object){

  ### Start Shiny APP
  {
    shinyApp(ui = MFN_manul_Shiny_ui(),
             server = MFN_manul_Shiny_server(object),
             options = list(host = "0.0.0.0",
                            launch.browser = T,
                            port = 6548))%>%
      runApp()->object


  }
}
