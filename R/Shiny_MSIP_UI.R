MSIP_shiny_ui <- function() {
  fluidPage(
    title = "MSIP V1",
    column(
      width = 5,
      #shinythemes::themeSelector(),
      h1("Metabolite table"),
      verbatimTextOutput(outputId = "test_info"),
      wellPanel(DTOutput(outputId = "metabolite_table")),
      wellPanel(
        h3("Help"),
        h4("1. Click Metablolite table to select a metabolite\n"),
        h4("2. Select isotope label count\n"),
        h4("3. Click colored peaks to select annotated fragment group\n"),
        h4(
          "4. Select a fragment to view its atom map to molecule, notice there are possibly multiple fragment for one peak(fragment group)\n"
        ),
        h4(
          "5. Molecule Structure show the labeled probability of atom in their filled color, and probability of maping to fragment atom in their border color"
        )

      )
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
            choices = "U"
          )
        ),
        column(
          3,
          selectInput(
            inputId = "select_iso_count",
            label = "Iso count",
            choices = paste0("M", 1:5)
          )
        ),
        column(
          3,
          checkboxInput(
            inputId = "spectra_show_noise",
            label = "show noise"
          )
        ),
        column(
          3,
          verbatimTextOutput(
            outputId = "compound_info"
          )
        )),
        plotlyOutput(outputId = "plotly_ms2_sp", height = "300px")
      ),
      navbarPage(
        title = NULL,
        tabPanel("Atom label probability", fluidRow(#h2("Molecule Structrue"),
          fluidPage(
            column(
              width = 3,
              wellPanel(
                checkboxInput(
                  inputId = "show_atom_id",
                  label = "Show Atom ID"
                ),
                br(),
                tableOutput(outputId = "atom_prob_table")
              )
              ),
            column(
              width = 8,
              #h5("Atom labeled probability"),
              fluidRow(
                align = "right",
                plotOutput(
                  outputId = "atom_prob_legend",
                  inline = F,
                  width = "200px",
                  height  = "70px"
                )
              ),
              fluidRow(visNetworkOutput(outputId = "mol_graph_atom_prob"))
            ),

          ))),
        tabPanel(
          "Atom map",
          fluidRow(
            column(
              width = 6,
              visNetworkOutput(outputId = "mol_graph_atom_map")
            ),
            column(
              width = 6,
              selectInput(
                inputId = "select_fragment_id",
                label = "Possible fragment structure",
                choices = NULL
              ),
              fluidRow(align = "center", textOutput(outputId = "frag_formula")),
              visNetworkOutput(outputId = "frag_graph")
            ))
          ),
        tabPanel(
          "Fragment map",
          fluidRow(
            column(
              3,
              actionButton(
                inputId = "Re_calc_button",
                label = "Re calculate"),
              DTOutput(
                outputId = "include_fragment_group"
              )
            ),
            column(
              9,
              plotOutput(
                  outputId = "heatmap_fg_map",
                  height = "500px",
                  width = "800px"
                )

            )
            )
          )
      )
    )





  )



}




MSIP_shiny_Acq_ui <- function() {

  fluidPage(
    verbatimTextOutput(outputId = "test_info"),
    column(
      width = 4,
      selectInput(
        inputId = "select_polarity",
        label = "Polarity",
        choices = c("Positive", "Negative")
      ),
      DTOutput(outputId = "feature_tab", height  = "700px"),
      fluidRow(
        column(
          width = 2,
          offset = 8,
          actionButton(inputId = "save_button", label = "Save")
        ),
        column(
          width = 2,
          offset = 0,
          actionButton(inputId = "quit_button", label = "Quit")
        )
      )
    ),
    #br(),
    column(
      width = 8,
      h1("chromatograms"),
      plotlyOutput(outputId = "feature_chrom", height  = "400px"),
      fluidRow(
        column(width = 3, (
          plotlyOutput(outputId = "p1", height  = "300px")
        )),
        column(width = 3, (
          plotlyOutput(outputId = "p2", height  = "300px")
        )),
        column(width = 3, (
          plotlyOutput(outputId = "p3", height  = "300px")
        )),
        column(width = 3, (
          plotlyOutput(outputId = "p4", height  = "300px")
        ))
      )
    ),


  )

}
