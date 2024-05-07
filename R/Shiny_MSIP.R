#' Title
#'
#' @param iso.msip.list from msip
#'
#' @return NULL
#' @export
#' @import  shiny
MSIP_shiny_start <- function(object){

  ### load temp
  {
    iso.msip.list <- object@statData$iso.msip.list
  }


  ### Start Shiny APP
  {
    shinyApp(ui = MSIP_shiny_ui(),
             server = MSIP_shiny_server(iso.msip.list),
             options = list(host = "0.0.0.0",
                            launch.browser = F,
                            port = 6548))


  }

}


MSIP_shiny_Acq <- function(object){



  ### Start Shiny APP
  {
    shinyApp(ui = MSIP_shiny_Acq_ui(),
             server = MSIP_shiny_Acq_server(object),
             options = list(host = "0.0.0.0",
                            launch.browser = T,
                            port = 6548))%>%
      runApp()


  }




}

