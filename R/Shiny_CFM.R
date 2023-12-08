#' Generate UI for shiny
#'
#' @return
#' @export
#' @import shiny plotly
#' @import shinythemes visNetwork
#' @examples
get_ui <- function(){

  ui <- fluidPage(
    theme =themeSelector() ,
    ### title--------
    titlePanel(title = "VisNet"),
    ### session1 ----
    fluidRow(
      mainPanel(h1("Inout var name"),
                textOutput(outputId = "show_input_name"),
                h1("Inout var value"),
                textOutput(outputId = "show_input_val")
                )
    ),

    ### session2 ----
    fluidRow(
      ### s2c1 ----
      column(3,
             visNetworkOutput(outputId = "Fragment_transition")
      ),
      ### s2c2 ----
      column(3,
             plotOutput(outputId = "Fragment_sdf")
      ),
      ### s2c3 ----
      column(3,
             mainPanel(textOutput("message_to_show"))
      )
    )


  )

  return(ui)
    ### session3----
    fluidRow()

}

#' Generate Server for shiny
#'
#' @return
#' @export
#' @import shiny visNetwork plotly
#' @examples
get_server <- function(shiny.server.data){

  server <- function(input, output) {


    output$Fragment_transition <- renderVisNetwork({

      visNetwork(shiny.server.data$fragment_transition$nodes[1:20,],
                 shiny.server.data$fragment_transition$edges)%>%
        visOptions(
          nodesIdSelection = list(
          enabled=T,
          selected = "Fragment001",
          useLabels = F
        ),
        manipulation = TRUE)
    })

    output$Fragment_sdf <- renderPlot({

      node.selected <- input$Fragment_transition_selected
      if (length(node.selected)==0) {
        node.selected <- "Fragment001"
      }
      smile.sdf <- shiny.server.data$fragment_transition$nodes.sdf[[node.selected]]
      plot(smile.sdf)

    })

    output$show_input_name <- renderPrint({
      names(input)
    })

    output$show_input_val <- renderPrint({

      #input$Fragment_transition_selected
      #input$Fragment_transition_selected  %in% shiny.server.data$fragment_transition$nodes$fragment_id
      length(input$Fragment_transition_selected)

    })


  }

  return(server)
}


get_shiny <- function(cfm.data = read_CFM_annotate_result()){


  ### data preprocess
  nodes <- cfm.data$fragment_define2%>%
    dplyr::mutate(id = fragment_id,
                  checked = check_smile(smile),
                  formula = get_smile_formula(smile),
                  label = formula
    )
  nodes.sdf <- smiles2sdf(nodes$smile%>%
                            `names<-`(nodes$fragment_id))
  edges <- cfm.data$fragment_transition%>%
    dplyr::mutate(arrows = "to",
                  formula = get_smile_formula(smile))
  shiny.server.data <- list()
  shiny.server.data$fragment_transition$nodes <- nodes
  shiny.server.data$fragment_transition$nodes.sdf <- nodes.sdf
  shiny.server.data$fragment_transition$edges <- edges


  ### shiny
  ui <- get_ui()
  sever <- get_server(shiny.server.data )
  shinyApp(ui,sever)
}

