MSIP_shiny_ui <- function(){


  fluidPage(
  title = "MSIP V1",
  column(
      width = 5,
      #shinythemes::themeSelector(),
      h1("Metabolite table"),
      wellPanel(
        DTOutput(outputId = "metabolite_table")
      ),
      wellPanel(
        h3("Help"),
        h4("1. Click Metablolite table to select a metabolite\n"),
        h4("2. Select isotope label count\n"),
        h4("3. Click colored peaks to select annotated fragment group\n"),
        h4("4. Select a fragment to view its atom map to molecule, notice there are possibly multiple fragment for one peak(fragment group)\n"),
        h4("5. Molecule Structure show the labeled probability of atom in their filled color, and probability of maping to fragment atom in their border color")

      )
    ),
  br(),
    column(
      width = 7,
      verbatimTextOutput(outputId = "test_info"),
      navbarPage(title = NULL,
                 tabPanel("Main result",
                          wellPanel(
                            fluidRow(
                              column(4,offset = 4,
                                selectInput(inputId = "select_sample",
                                            label = "Select sample",
                                            choices = "U" )
                              ),
                              column(4,
                                selectInput(inputId = "select_iso_count",
                                            label = "Select iso-labeled count",
                                            choices = paste0("M",1:5) )
                              )
                            ),
                            plotlyOutput(outputId = "plotly_ms2_sp")
                          ),
                          fluidRow(
                            h2("Molecule Structrue"),
                            fluidPage(
                              column(width = 2,
                                     fluidRow(
                                       align = "right",
                                       actionButton(inputId = "mol.ig.button",
                                                    label = "Show Atom ID",
                                                    height  = "50px")
                                     ),
                                     br(),
                                     column(width = 3,
                                            tableOutput(outputId = "atom_prob_table"))),
                              column(width = 5,
                                     #h5("Atom labeled probability"),
                                     fluidRow(
                                       align = "right",
                                       plotOutput(outputId = "atom_prob_legend",
                                                         inline = F,
                                                         width = "200px",
                                                         height  = "70px")
                                     ),
                                     fluidRow(
                                       visNetworkOutput(
                                                outputId = "mol_graph")
                                     )
                              ),
                              column(width = 5,
                                     selectInput(inputId = "select_fid",
                                                 label = "Possible fragment structure",
                                                 choices = paste0("M",1:5)),
                                     fluidRow(align = "center",
                                       textOutput(outputId = "frag_formula" )
                                     ),
                                     visNetworkOutput(
                                       outputId = "frag_graph")
                              )
                            )
                          )



                 ),
                 tabPanel("Atom map",
                          h1("Devoping...")
                 ),
                 tabPanel("Iso-form map",
                          h1("Devoping...")
                 )
      )
    )





  )



}




MSIP_shiny_Acq_ui <- function(){

  fluidPage(
    column(
      width = 4,
      selectInput(inputId = "select_polarity",
                         label = "Polarity",
                         choices = c("Positive","Negative")),
      DTOutput(outputId = "feature_tab",height  = "600px"),
      fluidRow(
        column(width = 2,offset = 8,
               actionButton(inputId = "save_button",label = "Save")),
        column(width = 2,offset = 0,
               actionButton(inputId = "quit_button",label = "Quit")))
    ),
    br(),
    column(
      width = 8,
      h1("chromatograms"),
      verbatimTextOutput(outputId = "test_info"),
      wellPanel(plotlyOutput(outputId = "feature_chrom",height  = "600px"))
    )

  )

}
