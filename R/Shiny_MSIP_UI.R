MSIP_shiny_ui <- function() {
  tags <- shiny::tags
  fluidPage(
    title = "MSIP",
    tags$style(
      HTML("
      .visnet-container {
        border: 1px solid grey;   /* Border color and thickness */
        border-radius: 10px;       /* Rounded corners */
        padding: 10px;             /* Space inside the border */

      }


    ")
    ),
    column(
      width = 5,
      fluidRow(
        column(width = 10,
               h1("MSIP"),
               h3("Mass Spectrometry based Isotopomers Parsing")
               ),
        column(width = 2,
               div(
                 shiny::actionLink("help_act", "Help", icon = icon("external-link"), style = "font-size: 20px; color: #000000;"),
                 style = "font-size: 20px; color: black; text-align: center;
                          justify-content: center; align-items: center;display: flex;height: 80px;"
               )
               )
      ),
      wellPanel(DT::DTOutput(outputId = "metabolite_table",height = "200%"))
    ),
    br(),
    column(
      width = 7,
      fluidRow(
        fluidRow(column(
          3,
          selectInput(
            inputId = "select_sample",
            label = "Sample",
            choices = "Glu"
          )
        ),
        column(
          3,
          selectInput(
            inputId = "select_iso_count",
            label = "Isotopologue",
            choices = paste0("M1")
          )
        ),
        column(
          3,
          checkboxInput(
            inputId = "spectra_show_rawdata",
            label = "show raw data"
          )
        ),
        column(
          3,
          shiny::verbatimTextOutput(
            outputId = "compound_info"
          )
        )),
        fluidRow(
          column(width = 10,
            plotly::plotlyOutput(outputId = "plotly_ms2_sp")
          ),
          column(width = 2,
            plotly::plotlyOutput(outputId = "plotly_natural_ratio",height  = "200px"),
            plotly::plotlyOutput(outputId = "plotly_fragment_iso_distribution",height  = "200px"),
            style = "padding:0px; margin:0px; background-color: #FFFFFF"
          )
        )
      ),
      navbarPage(
        title = NULL,
        {
          shiny::tabPanel(
            "Atom map",

            fluidRow(
              column(width = 6,
                     h3("Molecule structure and atom map")
              ),
              column(width = 6,
                     column(width = 6,
                            fluidRow(align = "left",
                                     selectInput(
                                       inputId = "select_fragment_id",
                                       label = "Possible fragment structure",
                                       choices = NULL
                                     )

                            )),
                     column(width = 6,
                            div(
                              shiny::textOutput(outputId = "frag_formula")),
                            style = "font-size: 20px; color: black; text-align: center;
                          justify-content: center; align-items: center;display: flex;height: 60px;"
                     )
              )
            ),
            fluidRow(
              column(
                width = 6,
                div(class = "visnet-container",
                    visNetwork::visNetworkOutput(outputId = "mol_graph_atom_map")
                )
              ),
              column(
                width = 6,
                div(class = "visnet-container",
                    visNetwork::visNetworkOutput(outputId = "frag_graph")
                )

              ))
          )
        },### ATOM MAP
        {
          shiny::tabPanel(
            "FAR map",
            fluidRow(
              fluidRow(
                column(width = 5,
                       shiny::sliderInput(inputId = "int_thresh",
                                          label = "Intensity threshold (log10)",
                                          ticks  = F,step = 0.01,
                                          min = 0,max = 10,value = 3 )
                ),
                column(width = 5,
                       shiny::sliderInput(inputId = "certainty_thresh",
                                          label = "Certainty threshold",
                                          ticks  = F,step = 0.01,
                                          min = 0,max = 1,value = 0.8)),
                column(width = 2,
                       shiny::actionButton(
                         inputId = "Re_calc_button",
                         label = "Re calculate")
                )
              ),
              shiny::plotOutput(
                outputId = "heatmap_fg_map",
                height = "800px",
                width = "800px"
              )
            )
          )
        },### FAR MAP
        {
          shiny::tabPanel(
            title = "Isotopomers",
            fluidPage(
              column(
                width = 7,
                #DT::DTOutput(outputId = "FSIS_table")
                DT::DTOutput(outputId = "isotopomer_table")
              ),
              column(
                width = 5,
                div(class = "visnet-container",
                    #visNetwork::visNetworkOutput(outputId = "Vis_isotopomer_set")
                    visNetwork::visNetworkOutput(outputId = "Vis_isotopomer")
                )

              )
            )
          )
        },### isotopomer
        {
          shiny::tabPanel("Atom label probability", fluidRow(#h2("Molecule Structrue"),
            fluidPage(
              column(
                width = 3,
                br(),
                br(),
                br(),
                br(),
                wellPanel(
                  checkboxInput(
                    inputId = "show_atom_id",
                    label = "Show Atom ID"
                  ),
                  br(),
                  shiny::tableOutput(outputId = "atom_prob_table")
                )
              ),
              column(
                width = 8,
                #h5("Atom labeled probability"),
                fluidRow(
                  align = "right",
                  shiny::plotOutput(
                    outputId = "atom_prob_legend",
                    inline = F,
                    width = "200px",
                    height  = "70px"
                  )
                ),
                fluidRow(
                  div(class = "visnet-container",
                      visNetwork::visNetworkOutput(outputId = "mol_graph_atom_prob")
                  )
                )
              ),

            )))
        },###
        {
          shiny::tabPanel(
            title = "Natural distribution",
            fluidRow(
              plotly::plotlyOutput(
                outputId = "pred_nat_prob",
                width = "600px",
                height = "600px",
              )
            )


          )
        }






      )
    ),
    tags$script(HTML("
    Shiny.addCustomMessageHandler('openPopup', function(url) {
      window.open(url, '_blank', 'width=1000,height=600,scrollbars=yes,resizable=yes');
    });
  "))





  )



}




MSIP_shiny_Acq_ui <- function() {

  fluidPage(
    title = "MSIP Acq vis",
    shiny::verbatimTextOutput(outputId = "test_info"),
    column(
      width = 4,
      selectInput(
        inputId = "select_polarity",
        label = "Polarity",
        choices = c("Positive", "Negative")
      ),
      DT::DTOutput(outputId = "feature_tab", height  = "700px"),
      fluidRow(
        column(
          width = 2,
          offset = 8,
          shiny::actionButton(inputId = "save_button", label = "Save")
        ),
        column(
          width = 2,
          offset = 0,
          shiny::actionButton(inputId = "quit_button", label = "Quit")
        )
      )
    ),
    #br(),
    column(
      width = 8,
      h1("chromatograms"),
      plotly::plotlyOutput(outputId = "feature_chrom", height  = "400px"),
      fluidRow(
        column(width = 3, (
          plotly::plotlyOutput(outputId = "p1", height  = "300px")
        )),
        column(width = 3, (
          plotly::plotlyOutput(outputId = "p2", height  = "300px")
        )),
        column(width = 3, (
          plotly::plotlyOutput(outputId = "p3", height  = "300px")
        )),
        column(width = 3, (
          plotly::plotlyOutput(outputId = "p4", height  = "300px")
        ))
      )
    )


  )

}



MFN_manul_Shiny_ui <- function(){


  fluidPage(
    title = "Mannual Atom transfer",
    #fluidRow(verbatimTextOutput(outputId = "test_output")),
    fluidRow(

      column(6,
             fluidRow(
               column(6,
                      textInput(inputId = "save_name",label = NULL,value = "Metabolic_flux_network",width = "50%")),
               column(6,
                      actionButton(inputId = "reverse_edge",label = "Reverse",width = "30%"))

             ),
             visNetwork::visNetworkOutput(
               outputId = "Metabolic_flux_network_vis",
               height = "800px", width = "100%"
             ),
             tags$style(
               HTML("
          #Metabolic_flux_network_vis {
            border: 3px solid #000; /* Black border */
            border-radius: 10px;    /* Rounded corners */
            padding: 5px;          /* Optional padding */
          }
        ")
             )
      ),
      column(6,
             fluidRow(
               column(6,
                 selectInput(inputId = "Atom_transfer_id",
                             label = "Select a transfer",
                             choices = "Atom_map1"),
               ),
               column(3,
                      actionButton(inputId = "save_buttion",label = "Export",
                                   width = "100%")
               )
             ),
             fluidRow(
               column(8,
                      visNetwork::visNetworkOutput(
                        outputId = "Atom_transfer_vis",
                        height = "400px", width = "100%"
                      )),
               column(4,
                      visNetwork::visNetworkOutput(
                        outputId = "Metabolite_isotopomer_statu_vis",
                        height = "400px", width = "100%"
                      )
                      )

             ),
             fluidRow(
               column(9,
                      verbatimTextOutput("Reaction_info")

                      )
             ),
             fluidRow(
                      DT::DTOutput(
                        "Metabolite_isotopomer_statu_table",
                        height = "400px"
                      )
             ),
             tags$style(
               HTML("
          #Atom_transfer_vis {
            border: 3px solid #000; /* Black border */
            border-radius: 10px;    /* Rounded corners */
            padding: 5px;          /* Optional padding */
          }
          #Metabolite_isotopomer_statu_vis {
            border: 3px solid #000; /* Black border */
            border-radius: 10px;    /* Rounded corners */
            padding: 5px;          /* Optional padding */
          }
        ")
             )
      )
    )

  )

}
