#' MSIP_start_server
#'
#' @title Msip Start Server
#' @description MSIP start server.
#' @returns NULL
#' @export
#'
#' @examples MSIP_start_server()
MSIP_start_server <- function(){


  message_with_time("Start FRP ")
  system2("cmd", args = "/k start C:/Users/91879/OneDrive/Software/EnvironmentPath/FRP/frp_0.61.2_windows_amd64/frpc.exe -c C:/Users/91879/OneDrive/Software/EnvironmentPath/FRP/frp_0.61.2_windows_amd64/frpc.toml")


  message_with_time("Start MSIP ")
  system2("cmd", args = c("/k", "start", "Rscript", "C:/Users/91879/OneDrive/Software/EnvironmentPath/FRP/frp_0.61.2_windows_amd64/MSIP.R"))


  message_with_time("Start MFNA ")
  system2("cmd", args = c("/k", "start", "Rscript", "C:/Users/91879/OneDrive/Software/EnvironmentPath/FRP/frp_0.61.2_windows_amd64/MFNA.R"))


}
