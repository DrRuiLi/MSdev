#' Generate UI for shiny
#'
#' @return null
#' @export
#' @import visNetwork MSCC DT

get_cfm_shiny_ui <- function(){

  ui <- fluidPage(
    #theme = bslib::bs_theme(bootswatch = "darkly"),
    theme =shinythemes::themeSelector() ,
    ### title--------
    titlePanel(title = "MSIP demo"),


    ### session1 show var info for dev----
    wellPanel(fluidRow(
      mainPanel(h1("Show info for debug"),
                actionButton(inputId = "test_button",
                             label = "Test"),
                verbatimTextOutput(outputId = "show_sp_clicked"),
      )
    )),

    ### session2 ----

      fluidRow(
        column(3,
               h1("Metabolites"),
               DTOutput(outputId = "metabolite_table")
        ),
        column(5,
               h1("Multi-CE Spectra"),
               plotlyOutput(outputId = "sp_plot")
        ),
        column(3,
                 h1("Molecular Structure"),
               textOutput("mol_formula"),
                 wellPanel(visNetworkOutput(outputId = "mol_graph",
                                            height = 300))
          )

      ),
    wellPanel(
      #h1("---"),
      br()),

    ### session3----
    fluidRow(
      ### s3c1 ----
      column(3,
             mainPanel(h1("Atom Tracing"),
                       DTOutput("tracing_tablle"))
      ),
      ### s3c2 ----
      column(5,
             h1("Fragment Transition"),
             wellPanel(
               visNetworkOutput(outputId = "Fragment_transition")
             )
      ),
      ### s3c3 ----
      column(3,
             h1("Fragment Structure"),
             textOutput("frag_formula"),
             #plotOutput(outputId = "Fragment_sdf")
             wellPanel(

               visNetworkOutput(outputId = "Fragment_sdf")
             )
      )

    )



  )

  return(ui)

}

