#' Generate Server for shiny
#'
#' @return
#' @export
#' @import visNetwork MSCC DT
#' @examples
get_cfm_shiny_server <- function(){

  server <- function(input, output) {


    ### var trans
    cfm_data <-reactiveVal(shiny_cfm_data_list[[1]])
    peak_mz_selected <- reactiveVal()
    frag_selected <- reactiveVal(
      as.character(shiny_cfm_data_list[[1]]$fragment_define2$fragment_id[1])
    )
    frag_highlight <- reactiveVal(
      as.character(shiny_cfm_data_list[[1]]$fragment_define2$fragment_id)
    )
    atom_root_highlight <- reactiveVal()


    ### click event metabolite table----
    observeEvent(input$metabolite_table_rows_selected,{

      cfm_data( shiny_cfm_data_list[[input$metabolite_table_rows_selected]])
      frag_selected(as.character(shiny_cfm_data_list[[input$metabolite_table_rows_selected]]$fragment_define2$fragment_id[1]))
      peak_mz_selected()
      })

    ### click event spectra----
    observeEvent( event_data("plotly_click",source = "sp_plot"),{

      peak_mz_selected(event_data("plotly_click",source = "sp_plot")$x )
      frag_selected({
        cfm_data()$peak.assignment%>%
          dplyr::filter(mz %in% peak_mz_selected(),
                        !is.na(fragment_id))%>%
          dplyr::ungroup()%>%
          dplyr::slice_max(fragment_score,n = 1,with_ties = T)%>%
          dplyr::pull(fragment_id)%>%
          as.character()
      })
      frag_highlight({
        cfm_data()$peak.assignment%>%
          dplyr::filter(mz %in% peak_mz_selected(),
                        !is.na(fragment_id))%>%
          dplyr::pull(fragment_id)%>%
          as.character()
      })
      atom_root_highlight({
        cfm_data()$tracing_stat%>%
          dplyr::filter(fragment_id %in% frag_highlight())%>%
          dplyr::pull(root_atom)%>%
          unlist()%>%unique()%>%na.omit()
      })
      })

    ### click event fragment transition ----
    observeEvent(input$Fragment_transition_selected,{

      frag.selected <-input$Fragment_transition_selected
      if (frag.selected == as.character(cfm_data()$fragment_define2$fragment_id[1])) {
        frag.selected <- " "
      }
      atom_root_highlight({
        cfm_data()$tracing_stat%>%
          dplyr::filter(fragment_id %in% frag.selected)%>%
          dplyr::pull(root_atom)%>%
          unlist()%>%unique()%>%na.omit()
      })

    })


    ### s1 text output show var------
    output$show_sp_clicked <- renderPrint({
      event_data("plotly_click",source = "sp_plot")
      atom_root_highlight()
    })


    ### s2  ----
    ### s2c1 metabolite table ----
    output$metabolite_table <- renderDT({
      shiny_metabolite_table%>%
        datatable( selection = list(
          mode = 'single', selected =1
        ))},
      server = F

    )

    ### s2c2 spectra ----
    output$sp_plot <- renderPlotly({

      peak_data <-cfm_data()$peak.assignment
      plotly_cfm_spectra(peak_data, peak_mz_selected() )%>%
        event_register("plotly_click")

    })

    ### s2c3 molecular structure
    output$mol_graph <- renderVisNetwork(
      {

        mol.graph <-cfm_data()$fragment_igraph[[1]]
        vis_sdf_igraph(mol.graph,highlight = atom_root_highlight() )

      }
    )
    output$mol_formula <- renderText({
      cfm_data()$fragment_define2$formula[1]
    })

    ### s3c1 Atom tracing ----
    {
      output$tracing_tablle <- renderDT({
        dt.data <- cfm_data()$tracing_stat%>%
          dplyr::mutate(selected = fragment_id %in% frag_highlight() )%>%
          dplyr::arrange(-selected)%>%
          dplyr::select(-selected)
        selected_no <- which(dt.data$fragment_id %in% frag_highlight() )
        dt.data%>%
          dplyr::select(-root_atom)%>%
          datatable(
            extensions = 'RowGroup',
            options = list(rowGroup = list(dataSrc = 1),
                           pageLength = 5),
            selection = list(
              mode = 'multiple', selected = selected_no
            )
          )

      })

    }


    ### s3c2 Fragment_transition----
    output$Fragment_transition <- renderVisNetwork({


      cfm_data()$transition_graph%>%
        vis_fragment_trans(frag.selected = frag_selected(),
                           highligted = frag_highlight())
    })


    ### s3c3 Fragment_sdf----
    output$Fragment_sdf <- renderVisNetwork({
      req(input$Fragment_transition_selected)
      vis_sdf_igraph( cfm_data()$fragment_igraph[[input$Fragment_transition_selected]]  )
    })
    output$frag_formula <- renderText({
      req(input$Fragment_transition_selected)
      cfm_data()$fragment_define2[input$Fragment_transition_selected,]$formula
    })


  }

  return(server)
}

