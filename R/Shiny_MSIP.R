#' Title
#'
#' @title Msip Shiny Start
#' @description MSIP shiny start.
#' @param object MSIP obj
#'
#' @return NULL
#' @export

MSIP_shiny_start <- function(object,port = NULL){


  ### Start Shiny APP
  {
    shinyApp(ui = MSIP_shiny_ui(),
             server = MSIP_shiny_server(object),
             options = list(host = "0.0.0.0",
                            port = port,
                            launch.browser = T))
    }

}


#' MSIP_shiny_Acq
#'
#' @title Msip Shiny Acq
#' @description MSIP shiny Acq.
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



shiny_test_fun <- function(input){

  print(names(input))
  print("\n")
  return(0)
}

MFN_manul_Shiny <- function(object,port = NULL){

  ### Start Shiny APP
  {
    shinyApp(ui = MFN_manul_Shiny_ui(),
             server = MFN_manul_Shiny_server(object),
             options = list(host = "0.0.0.0",
                            launch.browser = T,
                            port = port))%>%
      runApp()->object


  }
}
