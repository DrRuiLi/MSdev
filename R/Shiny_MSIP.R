#' Title
#'
#' @param iso.msip.list from msip
#'
#' @return NULL
#' @export
#' @import  shiny
MSIP_shiny_start <- function(iso.msip.list){

  ### load temp
  {
    load("temp.rda")
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


