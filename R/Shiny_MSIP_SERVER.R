MSIP_shiny_server <- function(object){


  iso.msip.list <- object@statData$MSIP$isotopologues_data
  all.sample<- object@sampleInfo%>%
    dplyr::filter(sample.type == "Sample")%>%
    dplyr::pull(sample.source)%>%
    unique()
  function(input, output, session) {




    output$test_info <- renderPrint({

      names(iso.msip()$MSIP_result[[1]])
      #vdata(mol.ig())
    })


    ### reactiveval
    {
      fid.selected <- reactiveVal()
      iso.msip <- reactiveVal()
      iso.count.selected <- reactiveVal()
      sample.selected <- reactiveVal()



      iso.sp.data <- reactiveVal()
      fg.selected <- reactiveVal()
      mol.ig <- reactiveVal()
      mol.atom.map <- reactiveVal()
      frag.ig <- reactiveVal()
      frag.atom.map <- reactiveVal()
      mol.atom.c.prob <- reactiveVal()
      mol.ig.show.label <- reactiveVal(F)

    }

    ### metabolite table
    {
      output$metabolite_table <- renderDT({
        get_iso_cfm_compound_info(iso.msip.list)%>%
          datatable( selection = list(
            mode = 'single', selected =1
          ))}
      )


      observeEvent(input$metabolite_table_rows_selected,{
        fid.selected(
          names(iso.msip.list)[input$metabolite_table_rows_selected]
        )
        iso.msip( iso.msip.list[[fid.selected()]])
        all.iso.count <- names(iso.msip()$MSIP_result)
        iso.count.selected(all.iso.count[1])
        updateSelectInput(inputId = "select_iso_count",
                          choices = all.iso.count,
                          selected = iso.count.selected())
        sample.selected(all.sample[1])
        updateSelectInput(inputId = "select_sample",
                          choices = all.sample,
                          selected = sample.selected())
        #mol.ig.show.label(F)
      })



    }

    ### selectInput, iso.count and sample
    {
    }

    ### Spectra
    {
      #input$select_iso_count
      observe({
        sample.selected(input$select_sample)
        iso.count.selected(input$select_iso_count)
        iso.sp.data(shiny_get_sp_data(iso_msip = iso.msip(),
                                      sample = sample.selected(),
                                      iso_count =  iso.count.selected()))
        #fg.selected(shiny_get_fg(sp.data = iso.sp.data(),
        #                         x = event_data("plotly_click",
        #                                        source = "plotly_ms2_sp",
        #                                        priority  = "event")$x))
        })

      output$plotly_ms2_sp <-  renderPlotly({

        shiny_plotly_iso_msip_spectra(iso.sp.data())

      })




    }

    ### Vis struc
    {

      observeEvent(fg.selected(),{

        #updateSelectInput(inputId = "select_fid",
        #                  choices = shiny_get_fid(iso_msip = iso.msip(),
        #                                          fg = fg.selected()))

      })


      observe({
        #mol.ig(iso.msip()$CFM_annotation@fragment_igraph[[1]])
        #frag.ig(shiny_get_fg_ig(iso_msip =iso.msip(),
        #                        fid = input$select_fid  ))
        #atom_map <- shiny_get_atom_map(iso.msip(),
        #                               input$select_fid ,
        #                               prob = T)
        #mol.atom.map(atom_map[[1]])
        #frag.atom.map(atom_map[[2]])
        #mol.atom.c.prob(iso.msip()$msip_data[[input$select_iso_count]]$c.prob)
      })

      output$mol_graph <- renderVisNetwork(
        {

          #shiny_vis_ig(ig = mol.ig() ,
          #             prob.border = mol.atom.map(),
          #             prob.fill =mol.atom.c.prob())
        }
      )
      output$frag_graph <- renderVisNetwork(
        {
#
          #shiny_vis_ig(frag.ig(),
          #             prob.border  = frag.atom.map() )

        }
      )
      output$frag_formula <- renderText({
        #shiny_get_frag_formula(iso_msip = iso.msip(),
        #                       fid =input$select_fid  )
      })


      output$atom_prob_legend <- renderPlot(
        draw(Legend(col_fun = shiny_col_ramp_atom_prob(),
                    direction = "h",border = T,
                    title = "Atom labeled\nprobability",
                    title_gp = gpar(fontsize = 15)),
                    #legend_height = unit(1,"inch"),
                    #legend_width = unit(2,"inch")),
             #x = unit(0, "npc"), y = unit(0, "npc"),
             just = "c")
      )

      output$atom_prob_table <- renderTable({
#
        #if (!length(mol.atom.c.prob())) {
        #  tb <- NULL
        #}else{
        #  tb <- data.frame(Atom = names(mol.atom.c.prob()),
        #                   Prob = mol.atom.c.prob())%>%
        #    dplyr::mutate(Prob = round(Prob,4))
        #}
        #tb


      })

      observeEvent(input$mol.ig.button,{
        #mol.ig.show.label(!mol.ig.show.label() )
        #updateActionButton(inputId = "mol.ig.button",
        #                   label = ifelse(
        #  mol.ig.show.label(),
        #  "Show Atom ID",
        #  "Hide Atom ID"))
        #igvd <- shiny_change_ig_label(ig = mol.ig(),
        #                      id = mol.ig.show.label())

        #visNetworkProxy(shinyId = "mol_graph")%>%
        #  visUpdateNodes(nodes = igvd,
        #                 updateOptions= F)

      })
    }






  }

}


#' Title
#'
#' @param object MSdev
#' @import DT
#' @return MSdev
MSIP_shiny_Acq_server <- function(object){

  acq.list <- object@statData$iso.acq.list
  function(input, output, session) {


    ### REACTVAL
    {
      acq.list.table <- reactiveVal()
      ### checkbox
      selected <- reactiveVal(
        list(Positive = acq.list$Positive$selected_to_acq,
             Negative = acq.list$Negative$selected_to_acq)
      )
      clicked_row <- reactiveVal(NA)
      xchrom <- reactiveVal()
      feature_id <- reactiveVal()
      fid_selected <- reactiveVal()
      fid_seed <- reactiveVal()


    }

    observeEvent(input$select_polarity,{
      acq.list.table(shiny_format_acq(
        acq.list[[input$select_polarity]],
        selected()[[input$select_polarity]]
      ))
      dataTableProxy("feature_tab")%>%
        showCols(show = 1:5,reset = T)
    })

    observeEvent(input$feature_tab_cell_clicked,{
      if (length(input$feature_tab_cell_clicked$col )) {
        if (input$feature_tab_cell_clicked$col== 5) {
          x <- selected()
          clicked_row(input$feature_tab_cell_clicked$row)
          x[[input$select_polarity]][input$feature_tab_cell_clicked$row] <-
            !x[[input$select_polarity]][input$feature_tab_cell_clicked$row]
          selected(x)

        }
      }



      x <- shiny_format_acq(
        acq.list[[input$select_polarity]],
        selected()[[input$select_polarity]]
      )
      dataTableProxy("feature_tab")%>%
        replaceData(x,resetPaging = F,clearSelection ="none" )%>%
        showCols(show = 1:5,reset = T)


    })

    observeEvent(input$feature_tab_rows_selected,{

      feature_id(
        acq.list.table()$feature_id[input$feature_tab_rows_selected]
      )
      cdf <- featureDefinitions(object@xcmsData[[paste0(input$select_polarity,"_Chromatograms")]])
      id <- match(feature_id(),cdf$feature_id)
      if (is.na(id)){
        xchrom(NA)
      }else{
        xchrom(object@xcmsData[[paste0(input$select_polarity,"_Chromatograms")]][
          id,
        ])
      }


    })

    output$feature_tab <- renderDT({
      acq.list.table()%>%
        datatable(escape = F,selection = "single",
                  extensions = c('RowGroup',"Scroller"),
                  options = list(rowGroup = list(dataSrc = 6),
                                 autoWidth = F,
                                 columnDefs = list(
                                   list(width = "500px",
                                        targets = 0:4)
                                 ),
                                 ordering = F,
                                 searching= F,
                                 deferRender = TRUE,
                                 scrollY = 500,
                                 scroller = TRUE))
    }
    )
    output$feature_chrom <- renderPlotly({

      if (!all(is.na(xchrom()))) {
        shiny_plotly_chrom(xchrom()  )
      }


    })



    output$test_info <- renderPrint({

      head(as.HTML.checkbox.logical( acq.list.table()$selected))
      input$feature_tab_cell_clicked
      input$feature_tab_cell_clicked$col
      head(selected()[[input$select_polarity]])
      dim(xchrom())
      input$feature_tab_rows_selected
      feature_id()
    })

    observeEvent(input$save_button,{
      object@statData$iso.acq.list[[input$select_polarity]]$selected_to_acq <-
        selected()[[input$select_polarity]]
    })

    observeEvent(input$quit_button,{
      object@statData$iso.acq.list$Positive$selected_to_acq <-
        selected()$Positive
      object@statData$iso.acq.list$Negative$selected_to_acq <-
        selected()$Negative
      stopApp(object)
    })


  }

}
