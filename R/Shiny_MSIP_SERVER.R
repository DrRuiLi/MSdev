MSIP_shiny_server <- function(object){



  all.sample <- names(na.omit(.get_MSIP_tracer(object)))
  iso.msip.list.statistic <- object@statData$MSIP$isotopologues_data
  function(input, output, session) {

    output$test_info <- renderPrint({



    })


    ### reactiveval
    {
      fid.selected <- reactiveVal()
      iso.msip.list <- reactiveVal( object@statData$MSIP$isotopologues_data)
      iso.msip <- reactiveVal( object@statData$MSIP$isotopologues_data[[1]])
      msip.core.data <- reactiveVal()
      iso.count.selected <- reactiveVal()
      sample.selected <- reactiveVal()



      iso.sp.data <- reactiveVal()
      fragment_group_selected <- reactiveVal()
      mol.ig <- reactiveVal()
      mol.atom.map <- reactiveVal()
      frag.ig <- reactiveVal()
      frag.atom.map <- reactiveVal()
      mol.atom.c.prob <- reactiveVal()
      mol.ig.show.label <- reactiveVal(F)
      sp.x.clicked <- reactiveVal()

    }


    ### input var
    {

      observe({
      ### metabolite_table
      input$metabolite_table_rows_selected

      ###selectInput
      input$select_sample
      input$select_iso_count
      input$select_fragment_id

      input$show_atom_id

      })


    }

    ### metabolite table
    {
      output$metabolite_table <- renderDT({
        get_iso_cfm_compound_info(iso.msip.list.statistic)%>%
          datatable( selection = list(
            mode = 'single', selected =1
          ))})


      observeEvent(input$metabolite_table_rows_selected,{
        fid.selected(
          names(iso.msip.list())[input$metabolite_table_rows_selected]
        )
        iso.msip( iso.msip.list()[[fid.selected()]])
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

      observeEvent(input$metabolite_table_rows_selected,{



      })


      output$compound_info <- renderText(
        shiny_isotopologues_info(iso_msip = iso.msip(),
                                 sample = sample.selected(),
                                 iso_count =  iso.count.selected())
      )


      observeEvent(
        {iso.msip()
          input$select_sample
          input$select_iso_count
        },
        {
          iso.sp.data(shiny_get_sp_data(iso_msip = iso.msip(),
                                        sample = input$select_sample,
                                        iso_count =  input$select_iso_count))
          sp.x.clicked(NULL)
          #message("input$select_sample")
          fragment_group_selected(
            shiny_get_fg(sp.data = iso.sp.data(),
                         x = sp.x.clicked()))
          msip.core.data(
            iso.msip()$MSIP_result[[input$select_iso_count]][[input$select_sample]]
            )
          mol.atom.c.prob({
            shiny_get_C_prob(iso.msip(),
                             input$select_iso_count,
                             input$select_sample)
          })
        })



    }

    ### Spectra
    {


      output$plotly_ms2_sp <-  renderPlotly({


        shiny_plotly_iso_msip_spectra(iso.sp.data())

      })

      observeEvent({
        event_data("plotly_click",
                   source = "plotly_ms2_sp",
                   priority  = "event")},
        {
          sp.x.clicked(event_data(
            "plotly_click",
            source = "plotly_ms2_sp",
            priority  = "event")$x) })
      observeEvent(sp.x.clicked(),{
        fragment_group_selected(
          shiny_get_fg(sp.data = iso.sp.data(),
                       x = sp.x.clicked()))
      })


      #input$select_iso_count
      #observe({
      #  sample.selected(input$select_sample)
      #  iso.count.selected(input$select_iso_count)
      #  iso.sp.data(shiny_get_sp_data(iso_msip = iso.msip(),
      #                                sample = sample.selected(),
      #                                iso_count =  iso.count.selected()))
      #  fragment_group_selected(shiny_get_fg(sp.data = iso.sp.data(),
      #                           x = event_data("plotly_click",
      #                                          source = "plotly_ms2_sp",
      #                                          priority  = "event")$x))
      #  })
#



    }

    ### atom prob
    {
      observeEvent(iso.msip(),{
        mol.ig(shiny_get_fg_ig(iso_msip =iso.msip(),
                               fid = 1  ))
      })

      observeEvent(fragment_group_selected(),{
        pf.id <- shiny_get_fid(iso_msip = iso.msip(),
                      fg = fragment_group_selected())
        selected <- ifelse(length(pf.id),pf.id[1], character(0))
        updateSelectInput(inputId = "select_fragment_id",
                          label = paste0(length(pf.id)," Possible fragment structure"),
                          choices = pf.id,
                          selected =selected )

      })

      observeEvent(input$select_fragment_id,{
        frag.ig(shiny_get_fg_ig(iso_msip =iso.msip(),
                                fid = input$select_fragment_id  ))
        atom_map <- shiny_get_atom_map(iso.msip(),
                                       input$select_fragment_id ,
                                       prob = T)
        mol.atom.map(atom_map[[1]])
        frag.atom.map(atom_map[[2]])
      })



      output$mol_graph_atom_prob <- renderVisNetwork(
        {
#
         #shiny_vis_ig(ig = mol.ig() ,
                      #prob.border = mol.atom.map(),
         #             prob.fill =mol.atom.c.prob())
         vis_sdf_igraph(mol.ig(),
                        show.id = input$show_atom_id,
                        prob.fill=mol.atom.c.prob())
        }
      )

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
        if (!length(mol.atom.c.prob())|all(is.na(mol.atom.c.prob()))){
          tb <- NULL
        }else{
          tb <- data.frame(Atom = names(mol.atom.c.prob()),
                           Prob = mol.atom.c.prob())%>%
            dplyr::mutate(Prob = round(Prob,4))
        }
        tb
      })

    }


    ### atom map
    {
      output$mol_graph_atom_map <- renderVisNetwork(
        {

          vis_sdf_igraph(mol.ig(),
                         show.id = input$show_atom_id,
                         prob.border =  mol.atom.map(),
                         prob.fill=mol.atom.c.prob())
        }
      )

      output$frag_formula <- renderText({
        shiny_get_frag_formula(
          iso_msip = iso.msip(),
          fid =input$select_fragment_id)
      })


      output$frag_graph <- renderVisNetwork(
        {
          vis_sdf_igraph(frag.ig(),
                         show.id = input$show_atom_id,
                       prob.border  = frag.atom.map() )
        }
      )


    }

    ### fg map
    {

      fg.include <- reactiveVal()
      observeEvent(msip.core.data(),{
        fg.include( shiny_get_fg_include(msip.core.data()) )
      })



      output$heatmap_fg_map <- renderPlot({

        fg.map <- msip.core.data()[["fg.map"]]

        if (all(is.na(fg.map))) {
          plot.new()
          title("no data")
        }else
          heatmap.fg.map(fg.map )

      })


      output$include_fragment_group <- renderDT({

        shiny_DT_fg_include(msip.core.data())%>%
          datatable(escape = F,
                    #rownames = F,
                    colnames = " ",
                    selection = "none",
                    extensions = c("Scroller"),
                    options = list(
                      columnDefs = list(
                        list(width = "100px",
                             targets = 0:1)
                      ),
                      info = F,
                      autoWidth = F,
                      ordering = F,
                      searching= F,
                      deferRender = TRUE,
                      scrollY = 500,
                      scroller = T))

      })

      observeEvent(input$include_fragment_group_cell_clicked,{

        if (length(input$include_fragment_group_cell_clicked$row )) {


          fg.include(shiny_update_fg_include(frag.include = fg.include(),
                                    input$include_fragment_group_cell_clicked$row))
          #print(x)
          x <- shiny_update_msip_core_data(msip.core.data(),fg.include())
          x <- shiny_DT_fg_include(x)
          dataTableProxy("include_fragment_group")%>%
            replaceData(x,resetPaging = F,clearSelection ="none" )
        }

      })

      observeEvent(input$Re_calc_button,{

        x <- shiny_update_msip_core_data(msip.core.data(),fg.include())
        x <- .re_calculate_isotopologues_label_fraction(x)


        iso.msip.list(shiny_update_iso_msip_list(
          iso.msip.list =iso.msip.list(),feature_id = fid.selected(),
          sample = sample.selected(),
          iso_count =  iso.count.selected(),
          msip.core.data  = x
        ))

        iso.msip( iso.msip.list()[[fid.selected()]])

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

  acq.list <- object@statData$MSIP$isotopologues_table
  rt.range <- range(rtime(object@xcmsData$PositiveMS1),
                    rtime(object@xcmsData$NegativeMS1))
  col.sample.source <- make_group_color(object@sampleInfo$sample.source,palette = "npg")
  function(input, output, session) {


    ### REACTVAL
    {
      acq.list.table <- reactiveVal()
      ### checkbox
      acq.selected <- reactiveVal(
        list(Positive = acq.list$Positive%>%
               dplyr::pull(selected_to_acq,name = feature_id),
             Negative = acq.list$Negative%>%
               dplyr::pull(selected_to_acq,name = feature_id))
      )
      clicked_row <- reactiveVal(NA)

      xchrom <- reactiveVal()
      feature_id <- reactiveVal()
      ratio_matrix <- reactiveVal()
      purity_matrix <- reactiveVal()


    }

    observeEvent(input$select_polarity,{
      acq.list.table(shiny_format_acq(
        acq.list[[input$select_polarity]],
        acq.selected()[[input$select_polarity]]
      ))

      dataTableProxy("feature_tab")%>%
        showCols(show = 1:6,reset = T)

      ratio_matrix(object@statData$MSIP$isotopologues_matrix$ratio_to_seed
                   [[input$select_polarity]]    )
      purity_matrix(object@statData$MSIP$isotopologues_matrix$ms1_purity
                   [[input$select_polarity]]    )
    })

    observeEvent(input$feature_tab_cell_clicked,{
      if (length(input$feature_tab_cell_clicked$col )) {
        if (input$feature_tab_cell_clicked$col== 6) {
          x <- acq.selected()
          clicked_row(input$feature_tab_cell_clicked$row)
          clicked_fid <- acq.list.table()$feature_id[ clicked_row() ]
          x[[input$select_polarity]][clicked_fid] <-
            !x[[input$select_polarity]][clicked_fid]
          acq.selected(x)

        }
      }



      x <- shiny_format_acq(
        acq.list[[input$select_polarity]],
        acq.selected()[[input$select_polarity]]
      )
      dataTableProxy("feature_tab")%>%
        replaceData(x,resetPaging = F,clearSelection ="none" )%>%
        showCols(show = 1:6,reset = T)


    })

    observeEvent(input$feature_tab_rows_selected,{

      feature_id(
        acq.list.table()$feature_id[input$feature_tab_rows_selected]
      )

      ### get chrom
      {
        cdf <- featureDefinitions(object@xcmsData[[paste0(input$select_polarity,"_Chromatograms")]])
        id <- match(feature_id(),cdf$feature_id)
        if (is.na(id)){
          xchrom(NA)
        }else{
          xchrom(object@xcmsData[[paste0(input$select_polarity,"_Chromatograms")]][
            id,
          ])
        }
      }

    })

    output$feature_tab <- renderDT({
      acq.list.table()%>%
        datatable(escape = F,
                  extensions = c('RowGroup',"Scroller"),
                  selection = list(
                    mode = 'single', selected = 1,target="row"
                  ),
                  options = list(rowGroup = list(dataSrc = 7),
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
    })
    output$feature_chrom <- renderPlotly({

      if (!all(is.na(xchrom()))) {
        shiny_plotly_chrom(xchrom(),col.map = col.sample.source,rtr = rt.range  )
      }


    })
    output$p1 <- renderPlotly({
      if (!all(is.na(xchrom())))
      shiny_plotly_feature_int(xchrom(),col.map = col.sample.source)

    })
    output$p2 <- renderPlotly({
      if (!all(is.na(xchrom())))
      shiny_plotly_feature_ratio(xchrom(),
                                 col.map = col.sample.source,
                                 ratio_to_seed = ratio_matrix()[feature_id(),] )

    })
    output$p3 <- renderPlotly({
      if (!all(is.na(xchrom())))
        shiny_plotly_feature_purity(xchrom(),
                                   col.map = col.sample.source,
                                     ms1_purity = purity_matrix()[feature_id(),] )
    })
    output$p4 <- renderPlotly({

    })
    output$test_info <- renderPrint({

     #object@statData$MSIP$isotopologues_table$Positive%>%
     #   dplyr::pull(selected_to_acq,name = feature_id)
      #head(ratio_matrix())
    })


    observeEvent(input$save_button,{

      object@statData$MSIP$isotopologues_table$Positive$selected_to_acq <-
        acq.selected()$Positive[object@statData$MSIP$isotopologues_table$Positive$feature_id]
      object@statData$MSIP$isotopologues_table$Negative$selected_to_acq <-
        acq.selected()$Negative[object@statData$MSIP$isotopologues_table$Negative$feature_id]

    })

    observeEvent(input$quit_button,{
      #object@statData$iso.acq.list$Positive$selected_to_acq <-
      #  selected()$Positive
      #object@statData$iso.acq.list$Negative$selected_to_acq <-
      #  selected()$Negative
      object@statData$MSIP$isotopologues_table$Positive$selected_to_acq <-
        acq.selected()$Positive[object@statData$MSIP$isotopologues_table$Positive$feature_id]
      object@statData$MSIP$isotopologues_table$Negative$selected_to_acq <-
        acq.selected()$Negative[object@statData$MSIP$isotopologues_table$Negative$feature_id]

      stopApp(object)
    })


  }

}
