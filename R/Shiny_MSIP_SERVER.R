MSIP_shiny_server <- function(object){


  ### for DT::renderDT, to avoid fresh
  object.temp <- object

  function(input, output, session) {


    # Trigger the pop-up when the app starts
    #session$sendCustomMessage("openPopup", "https://drruili.github.io/MSIP/articles/msip_shiny_demo.html")

    output$test_info <- renderPrint({
      FSIS_selected()
    })

    observeEvent(input$help_act, {
      browseURL("https://drruili.github.io/MSIP/articles/msip_shiny_demo.html")  # Replace with your URL or local file
    })


    ### reactiveval
    {
      fid.selected <- reactiveVal()
      iso.data.list <- reactiveVal( object@statData$MSIP$isotopologues_data)
      iso.data <- reactiveVal()
      msip.core.data <- reactiveVal()



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
      output$metabolite_table <- DT::renderDT({

        message_with_time("metabolite_table")
        get_MSIP_compound_info(object.temp)%>%
          DT::datatable(options = list(columns =list(orderable = F),
                                       pageLength = 25),
                    selection = list(
            mode = 'single', selected =1
          ))})


      observeEvent(input$metabolite_table_rows_selected,{
        message_with_time("metabolite_table_rows_selected: ",input$metabolite_table_rows_selected)
        fid.selected(
          names(iso.data.list())[input$metabolite_table_rows_selected]
        )
        iso.data( iso.data.list()[[fid.selected()]])


        ms2.count <- iso.data()$compound_info$ms2_count


        all.iso.count <- setdiff(rownames(ms2.count),"M0")
        all.iso.count <- levels(groupStringFactor(all.iso.count))
        updateSelectInput(inputId = "select_iso_count",
                          choices = all.iso.count,
                          selected = all.iso.count[1])

        all.sample <- colnames(ms2.count)
        all.sample <- levels(groupStringFactor(all.sample))
        updateSelectInput(inputId = "select_sample",
                          choices = all.sample,
                          selected = all.sample[1])



        iso.sp.data(shiny_get_sp_data(iso_data =  iso.data(),
                                      sample = input$select_sample,
                                      iso_count =  input$select_iso_count))
        fragment_group_selected(
          shiny_get_fg(sp.data = iso.sp.data(),
                       x = sp.x.clicked()))
      })



    }



    ### selectInput, iso_count and sample
    {

      output$compound_info <- renderText(

        shiny_isotopologues_info(iso_data = iso.data(),
                                 sample = input$select_sample,
                                 iso_count =  input$select_iso_count)
      )



      natural.ratio <- reactiveVal(0)

      observeEvent(
        { iso.data()
          input$select_sample
          input$select_iso_count
        },
        {
          message_with_time("shiny_get_sp_data")
          iso.sp.data(shiny_get_sp_data(iso_data =  iso.data(),
                                        sample = input$select_sample,
                                        iso_count =  input$select_iso_count))

          natural.matrix <- iso.data()[["compound_info"]][["natural_matrix"]]
          natural.ratio(
            get_matrix_value_fill_with_NA(natural.matrix,
                                          input$select_iso_count, input$select_sample)

          )
          sp.x.clicked(NULL)
          #message("input$select_sample")
          message_with_time("shiny_get_fg")
          fragment_group_selected(
            shiny_get_fg(sp.data = iso.sp.data(),
                         x = sp.x.clicked()))
          msip.core.data(
            iso.data()$MSIP_result[[input$select_iso_count]][[input$select_sample]]
            )
          ###
          mol.atom.c.prob({
            shiny_get_C_prob(msip.core.data())
          })
        })



    }

    ### Spectra
    {


      output$plotly_ms2_sp <-  plotly::renderPlotly({

        message_with_time("plotly_ms2_sp")

        shiny_plotly_iso_data_spectra(iso.sp.data(),
                                      show.rawdata = input$spectra_show_rawdata)

      })

      observeEvent({
        event_data("plotly_click",
                   source = "plotly_ms2_sp",
                   priority  = "event")},
        {
          message_with_time("plotly_click")
          sp.x.clicked(event_data(
            "plotly_click",
            source = "plotly_ms2_sp",
            priority  = "event")$x) })

      observeEvent(sp.x.clicked(),{
        message_with_time("sp.x.clicked")
        fragment_group_selected(
          shiny_get_fg(sp.data = iso.sp.data(),
                       x = sp.x.clicked()))

      })


      output$plotly_natural_ratio <-  plotly::renderPlotly({
        message_with_time("plotly_natural_ratio")

        shiny_plotly_natural_ratio(natural.ratio())

      })

      output$plotly_fragment_iso_distribution <-  plotly::renderPlotly({
        message_with_time("shiny_get_frag_iso_distribution")
        frag.iso.dis <- shiny_get_frag_iso_distribution(msip.core.data(),
                                                        fragment_group_selected())
        message_with_time("shiny_plotly_iso_distribution")
        shiny_plotly_iso_distribution(frag.iso.dis)

      })



    }

    ### atom prob
    {
      observeEvent({
        iso.data()
        fragment_group_selected()
      },{
        message_with_time("shiny_get_mol_ig :",fragment_group_selected())
        mol.ig(shiny_get_fg_ig(iso_data = iso.data(),
                               fg = fragment_group_selected(),
                               fid = "seed" ))
      })

      observeEvent(fragment_group_selected(),{
        message_with_time("shiny_get_fid")
        possible.fragment.id <- shiny_get_fid(iso_data  = iso.data(),
                      fg = fragment_group_selected())
        selected <- ifelse(length(possible.fragment.id),possible.fragment.id[1], character(0))
        updateSelectInput(inputId = "select_fragment_id",
                          label = paste0(length(possible.fragment.id)," Possible fragment structure"),
                          choices = possible.fragment.id,
                          selected =selected )

      })

      observeEvent(input$select_fragment_id,{
        message_with_time("select_fragment_id")
        frag.ig(shiny_get_fg_ig(iso_data  =iso.data(),
                                fg = fragment_group_selected(),
                                fid = input$select_fragment_id  ))
        atom_map <- shiny_get_atom_map(iso.data(),
                                       fg = fragment_group_selected(),
                                       input$select_fragment_id ,
                                       prob = T)
        mol.atom.map(atom_map[[1]])
        frag.atom.map(atom_map[[2]])
      })



      output$mol_graph_atom_prob <- visNetwork::renderVisNetwork(
        {
          message_with_time("mol_graph_atom_prob")

          shiny_vis_sdf_igraph(mol.ig(),
                          show_id = input$show_atom_id,
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

        message_with_time("atom_prob_table")
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
      output$mol_graph_atom_map <- visNetwork::renderVisNetwork(
        {
          message_with_time("mol_graph_atom_map")

          shiny_vis_sdf_igraph(mol.ig(),
                         show_id = input$show_atom_id,
                         prob.border =  mol.atom.map(),
                         prob.fill= 0)
        }
      )

      output$frag_formula <- renderText({
        message_with_time("frag_formula")
        shiny_get_frag_formula(
          iso_data = iso.data(),
          fg = fragment_group_selected(),
          fid =input$select_fragment_id)
      })


      output$frag_graph <- visNetwork::renderVisNetwork(
        {
          message_with_time("frag_graph")
          shiny_vis_sdf_igraph(frag.ig(),
                         show_id = input$show_atom_id,
                        prob.border  = frag.atom.map() ,
                        prob.fill= 0)
        }
      )


    }

    ### fg map and re-calc
    {

      fg.include <- reactiveVal()
      observeEvent(msip.core.data(),{

        message_with_time("shiny_get_fg_include")
        fg.include( shiny_get_fg_include(msip.core.data()) )

        message_with_time("shiny_get_int_thresh")
        msip.core.int_thresh <- shiny_get_int_thresh(msip.core.data())
        updateSliderInput(inputId = "int_thresh",
                          value = msip.core.int_thresh)
        updateSliderInput(inputId = "certainty_thresh",
                          value = msip.core.data()@solve$certainty_thresh)


      })



      output$heatmap_fg_map <- renderPlot({

        message_with_time("shiny_heatmap_fgmap")
        shiny_heatmap_fgmap(msip.core.data())
      },
      height = reactive(nrow(msip.core.data()@FG_map@FG.atom.matrix)*50+50)
      )


      if(F){
        output$include_fragment_group <- DT::renderDT({

          message_with_time("include_fragment_group")
          shiny_DT_fg_include(msip.core.data())%>%
            DT::datatable(escape = F,
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

          message_with_time("observeEvent316")
          if (length(input$include_fragment_group_cell_clicked$row )) {


            fg.include(shiny_update_fg_include(frag.include = fg.include(),
                                               input$include_fragment_group_cell_clicked$row))
            #print(x)
            x <- shiny_update_msip_core_data(msip.core.data(),fg.include())
            x <- shiny_DT_fg_include(x)
            DT::dataTableProxy("include_fragment_group")%>%
              DT::replaceData(x,resetPaging = F,clearSelection ="none" )
          }

        })
      }


      observeEvent(input$Re_calc_button,{

        message_with_time("observeEvent333")
        shinybusy::show_modal_spinner()

        #x <- shiny_update_msip_core_data(msip.core.data(),
        #                                 fg.include())
        message_with_time("ReSolve")
        wf <- get_MSIP_weight_fun(object)
        msip.core.data(MSIPCore_solve(msip.core.data(),weight_fun = wf,
                                      int_thresh = 10^input$int_thresh,
                            certainty_thresh = input$certainty_thresh))


        message_with_time("update_result")
        iso.data.list(shiny_update_iso_data_list(
          iso.data.list =iso.data.list(),
          feature_id = fid.selected(),
          sample =input$select_sample,
          iso_count =  input$select_iso_count,
          msip.core.data  = msip.core.data()
        ))

        iso.data( iso.data.list()[[fid.selected()]])
        shinybusy::remove_modal_spinner()

      })

    }


    ### predict and natural prob
    {
      output$pred_nat_prob <- plotly::renderPlotly({

        message_with_time("pred_nat_prob")
        plotly_MSIPCore_pred_nature_prob(msip.core.data())


      })

    }



    ### FSIS
    {
      FSIS_table <- reactiveVal()
      message_with_time("shiny_get_FSIS_table")
      observeEvent(msip.core.data(),{

        FSIS_table(shiny_get_FSIS_table(msip.core.data()))
      })

      output$FSIS_table <- DT::renderDT({
        FSIS_table()%>%
          DT::datatable(options = list(
            columns =list(orderable = F),
            dom = 't',   # 't' means only the table body is displayed
            paging = FALSE,  # Disable pagination
            searching = FALSE,  # Disable search box
            info = FALSE  # Disable table information
          ),
                    selection = list(
                      mode = 'single', selected =1
                    ))
      })

      FSIS_selected <- reactiveVal()
      observeEvent(input$FSIS_table_rows_selected,{

        FSIS_selected(
          FSIS_table()$FSIS[input$FSIS_table_rows_selected]
        )


      })


      message_with_time("vis_MSIPcore_isotopoer_set")
      output$Vis_isotopomer_set <- visNetwork::renderVisNetwork({
        vis_MSIPcore_isotopoer_set(msip.core.data(),
                                   mol.ig(),FSIS_selected())
      })




    }


    ### isotopomer
    {
      isotopomer_table <- reactiveVal()
      observeEvent(msip.core.data(),{
        message_with_time("shiny_get_isotopomer_table")
        isotopomer_table(shiny_get_isotopomer_table(msip.core.data()))
      })

      output$isotopomer_table <- DT::renderDT({

        message_with_time("isotopomer_table")
        x <- isotopomer_table()
        colnames(x)[4] <- "Mixed\nisotopomers"
        x%>%
          DT::datatable(
            rownames = FALSE,
            options = list(
            columns =list(orderable = F),
            dom = 't',   # 't' means only the table body is displayed
            paging = FALSE,  # Disable pagination
            searching = FALSE,  # Disable search box
            info = FALSE  # Disable table information
          ),
          selection = list(
            mode = 'single', selected =1
          ))
      })

      isotopomer_selected <- reactiveVal()
      observeEvent(input$isotopomer_table_rows_selected,{

        isotopomer_selected(
          isotopomer_table()$Isotopomer[input$isotopomer_table_rows_selected]
        )


      })


      output$Vis_isotopomer <- visNetwork::renderVisNetwork({
        message_with_time("Vis_isotopomer")
        vis_MSIPcore_isotopomer(msip.core.data(),
                                   mol.ig(),isotopomer_selected())
      })




    }

  }

}


#' MSIP_shiny_Acq_server
#'
#' @param object MSdev
#' @return MSdev
MSIP_shiny_Acq_server <- function(object){

  acq.list <- object@statData$MSIP$isotopologues_table
  rt.range <- range(rtime(object@xcmsData$PositiveMS1),
                    rtime(object@xcmsData$NegativeMS1))
  col.sample.source <- make_group_color(object@sampleInfo$sample.source,palette = "npg")
  xcms.chrom.data <- object@xcmsData[c("Positive_Chromatograms",
                                        "Negative_Chromatograms")]
  message_with_time("onDiskData_retrieve...")
  xcms.chrom.data <- lapply(xcms.chrom.data,
                            onDiskData_retrieve  )
  message_with_time("Done")
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

      xchrom.selected.polarity <- reactiveVal()
      xchrom <- reactiveVal()
      feature_id <- reactiveVal()
      ratio_matrix <- reactiveVal()
      purity_matrix <- reactiveVal()


    }

    observeEvent(input$select_polarity,{

      message_with_time("observeEvent 402")
      acq.list.table(shiny_format_acq(
        acq.list[[input$select_polarity]],
        acq.selected()[[input$select_polarity]]
      ))

      DT::dataTableProxy("feature_tab")%>%
        DT::showCols(show = 1:6,reset = T)

      ratio_matrix(object@statData$MSIP$isotopologues_matrix$ratio_to_seed
                   [[input$select_polarity]]    )
      purity_matrix(object@statData$MSIP$isotopologues_matrix$ms1_purity
                   [[input$select_polarity]]    )

      xchrom.selected.polarity(xcms.chrom.data[[paste0(input$select_polarity,"_Chromatograms")]])
    })

    observeEvent(input$feature_tab_cell_clicked,{

      message_with_time("observeEvent 419")
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
      DT::dataTableProxy("feature_tab")%>%
        DT::replaceData(x,resetPaging = F,clearSelection ="none" )%>%
        DT::showCols(show = 1:6,reset = T)


    })

    observeEvent(input$feature_tab_rows_selected,{

      message_with_time("observeEvent 452")

      feature_id(
        acq.list.table()$feature_id[input$feature_tab_rows_selected]
      )

      ### get chrom
      {
        cdf <- featureDefinitions(xchrom.selected.polarity())
        pdf <- pData(xchrom.selected.polarity())
        id <- match(feature_id(),cdf$feature_id)
        if (is.na(id)){
          xchrom(NA)
        }else{
          xchrom(xchrom.selected.polarity()[id, pdf$sample.type != "Blank" ,drop = F])
        }
      }
      message_with_time("observeEvent 452 done")

    })

    output$feature_tab <- DT::renderDT({


      message_with_time("feature_tab")

      acq.list.table()%>%
        DT::datatable(escape = F,
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
    output$feature_chrom <- plotly::renderPlotly({
      message_with_time("feature_chrom")

      if (!all(is.na(xchrom()))) {
        shiny_plotly_chrom(xchrom(),col.map = col.sample.source,rtr = rt.range  )
      }


    })
    output$p1 <- plotly::renderPlotly({
      message_with_time("p1")
      if (!all(is.na(xchrom())))
      shiny_plotly_feature_int(xchrom(),col.map = col.sample.source)

    })
    output$p2 <- plotly::renderPlotly({
      message_with_time("p2")
      if (!all(is.na(xchrom())))
      shiny_plotly_feature_ratio(xchrom(),
                                 col.map = col.sample.source,
                                 ratio_to_seed = ratio_matrix()[feature_id(),] )

    })
    output$p3 <- plotly::renderPlotly({
      message_with_time("p3")
      if (!all(is.na(xchrom())))
        shiny_plotly_feature_purity(xchrom(),
                                   col.map = col.sample.source,
                                     ms1_purity = purity_matrix()[feature_id(),] )
    })
    output$p4 <- plotly::renderPlotly({
      message_with_time("p4")

    })
    output$test_info <- renderPrint({

     #object@statData$MSIP$isotopologues_table$Positive%>%
     #   dplyr::pull(selected_to_acq,name = feature_id)
      #head(ratio_matrix())
    })


    observeEvent(input$save_button,{
      message_with_time("observeEvent 535")

      object@statData$MSIP$isotopologues_table$Positive$selected_to_acq <-
        acq.selected()$Positive[object@statData$MSIP$isotopologues_table$Positive$feature_id]
      object@statData$MSIP$isotopologues_table$Negative$selected_to_acq <-
        acq.selected()$Negative[object@statData$MSIP$isotopologues_table$Negative$feature_id]

    })

    observeEvent(input$quit_button,{
      message_with_time("observeEvent 545")
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




MFN_manul_Shiny_server <- function(object){


  function(input, output, session){



    ### initial
    {
      mfn.r <- vdata(object)%>%
        dplyr::filter(node.type == "Reaction")
      mfn.c <- vdata(object)%>%
        dplyr::filter(node.type == "Compound")

      mfn.selected.node.type <- reactiveVal("Reaction")
      mfn.selected.reaction <- reactiveVal(mfn.r$id[1])
      mfn.selected.compound <- reactiveVal(mfn.c$id[1])


    }

    ### test
    {
      output$test_output <- renderPrint({


        #edata(mfn_for_update())$id
        #names(input$Metabolic_flux_network_vis_edges)

      })
    }



    output$Metabolic_flux_network_vis <- visNetwork::renderVisNetwork({

      message_with_time("render_Metabolic_flux_network_vis")
      mfn <- mfn_for_vis()
      shinybusy::show_modal_spinner(spin = "circle", color = "#007bff")  # Show spinner
      vis_Metabolic_flux_network(mfn)%>%
        visOptions(width = "100%",height = "100%",
                   #nodesIdSelection = list(enabled=T, selected= "C00631" ),
                   manipulation = list(
                     enabled = TRUE,
                     addEdge = F,
                     deleteEdge = F,
                     editEdge = F,
                     addNode = F,
                     deleteNode  = T,
                     editNode = F
                   ))%>%
        visEvents(
          selectNode = "function(properties) {
        Shiny.onInputChange('Metabolic_flux_network_vis_selected_node', properties.nodes);
      }",
          stabilizationIterationsDone = "function() {
        Shiny.setInputValue('network_ready', true, {priority: 'event'});
      }"
        )%>%
        visInteraction(selectConnectedEdges  = T)

    })


    observeEvent(input$network_ready, {
      shinybusy::remove_modal_spinner()  # Hide the modal spinner when the network is ready
    })
    ### Metabolic_flux_network_vis select node
    {


      observeEvent(input$Metabolic_flux_network_vis_selected_node,{


        ### reactiveVal in
        {
          mfn <- mfn_for_update()
        }

        ### update
        {
          vda <- vdata(mfn)%>%
            dplyr::filter(id == input$Metabolic_flux_network_vis_selected_node)

          mfn.selected.node.type(vda$node.type)
          if (vda$node.type =="Reaction") {
            mfn.selected.reaction(vda$id)
          }
          if (vda$node.type =="Compound") {
            mfn.selected.compound(vda$id)
          }

        }


      })




    }


    ### Reaction_atom_transfer select
    {

      rat_for_vis <- reactiveVal()
      rat_for_update <- reactiveVal()

      observeEvent(mfn.selected.reaction(),{

        ### reactiveVal in
        {
          rid <- mfn.selected.reaction()
          mfn <- mfn_for_update()

        }

        ### reactiveVal update
        {
          rat <- V(mfn@metabolic_network)[[rid]]$Reaction_atom_transfer
          rat_for_vis(rat)
          rat_for_update(rat)
          }

      })

      mig_for_vis <- reactiveVal()
      observeEvent(mfn.selected.compound(),{

        ### reactiveVal in
        {
          nid <- mfn.selected.compound()
          mfn <- mfn_for_update()

        }

        ### reactiveVal update
        {
          mig <- V(mfn@metabolic_network)[[nid]]$Molecule_igraph
          mig_for_vis(mig)
        }

      })




    }


    output$Atom_transfer_vis <- visNetwork::renderVisNetwork({
      message_with_time("render_Atom_transfer_vis")
      ### reactiveVal
      {
        rat <- rat_for_vis()
      }

      vis_Reaction_atom_transfer(rat)%>%
        visOptions(manipulation = list(
          enabled = TRUE,
          addEdge = htmlwidgets::JS(
            "function(data, callback) {
            Shiny.onInputChange('Atom_transfer_vis_edited', new Date().getTime() )
            // Check if edge is valid
            if (data.from === data.to) {
              alert('Edges cannot connect a node to itself.');
              return;
            }
            // Set custom attributes for the edge
            data.color = {color: '#54BFBF'};  // Set edge color
            data.dashes = [10, 20];       // Dashed pattern
            data.label = 'Custom Transfer'; // Edge label
            data.width = 8;
            data.arrows = 'to';
            data.connect_type = 'custom_atom_transfer';
            callback(data);              // Finalize the edge addition
          }"
          ),
          deleteEdge = htmlwidgets::JS(
            "function(data, callback) {
            Shiny.onInputChange('Atom_transfer_vis_edited', new Date().getTime() )
            callback(data);              // Finalize the edge addition
          }"
          ),
          editEdge = htmlwidgets::JS(
            "function(data, callback) {
            Shiny.onInputChange('Atom_transfer_vis_edited', new Date().getTime() )
            callback(data);              // Finalize the edge addition
          }"
          ),
          addNode = F,deleteNode  = F,editNode = F

        ))

    })






    ### Reaction_atom_transfer edit
    {

      observeEvent(input$Atom_transfer_vis_edited,{
        message_with_time("Atom_transfer_update")
        visNetworkProxy("Atom_transfer_vis")%>%
          visGetEdges()

      })

      observeEvent(input$Atom_transfer_vis_edges,{
        message_with_time("Atom_transfer_update_record")

        ### reactiveVal
        {
          mat <- mat_for_update()
          mat.id <- input$Atom_transfer_id
          visGetEdges <- input$Atom_transfer_vis_edges
        }
        mat <- get_mat_from_visGetEdges(mat = mat,id = mat.id,visGetEdges =  visGetEdges )



        mat_for_update(mat)
      })
    }



    ### MFN save and export
    {
      mfn_for_vis <- reactiveVal(object)
      mfn_for_update <- reactiveVal(object)

    }


    ### MFN update
    {
      observeEvent(input$Metabolic_flux_network_vis_edited,{

        message_with_time("Metabolic_flux_network_vis_edited")
        visNetworkProxy("Metabolic_flux_network_vis")%>%
          visGetEdges()

      })


      observeEvent(input$Metabolic_flux_network_vis_edges,{
        message_with_time("Metabolic_flux_network_vis_update_record")

        ### reactiveVal
        {
          mfn <- mfn_for_update()
          visGetEdges <- input$Metabolic_flux_network_vis_edges
        }
        mfn <- Metabolic_flux_network_update_from_visGetEdges(mfn = mfn,
                                                              visGetEdges =  visGetEdges )


        mfn_for_update(mfn)

      })

      observeEvent(input$reverse_edge,{


        ### input
        {
          mfn <- mfn_for_update()
          edge.selected <- Metabolic_flux_network_vis_selected_edge()
        }

        ### output
        {

          mfn <- Metabolic_flux_network_reverse(mfn,edge.selected)
          mfn_for_update(mfn)
          rat_for_vis(
            mat <- E(mfn@metabolic_network)[[edge.selected]]$atom_transfer)
        }

        ### update vis
        {


        }

      })

      observeEvent(mfn_for_update(),{
        #visNetworkProxy("Metabolic_flux_network_vis")%>%
        #  visUpdateEdges(edata(mfn_for_update()))
      })
    }




    ### Metabolic_flux_network_vis select edge
    {
      Metabolic_flux_network_vis_selected_edge <- reactiveVal(1)
      Atom_transfer_id_selected <- reactiveVal(1)

      observeEvent(input$Metabolic_flux_network_vis_selected_edge,{

        message_with_time("Metabolic_flux_network_vis_selected_edge")

        ### update mat_for_update to MFN
        {
          if (!identical( mat_for_update(), rat_for_vis() )) {
            #print("update mat_for_update to MFN")
              mfn_for_update({
                mfn <- mfn_for_update()
                eid <- Metabolic_flux_network_vis_selected_edge()
                mat <- mat_for_update()
                E(mfn@metabolic_network)[[eid]]$atom_transfer <- mat
                mfn
              })
              #shiny_test_fun(mfn)
          }

        }


        Metabolic_flux_network_vis_selected_edge(input$Metabolic_flux_network_vis_selected_edge)


      })

      observeEvent(Metabolic_flux_network_vis_selected_edge(),{
        message_with_time("Metabolic_flux_network_vis_selected_edge_update_MAT")

        ### reactiveVal
        {
         # MFN <- mfn_for_update()
         # edge.selected <- Metabolic_flux_network_vis_selected_edge()
        }
#
    #   if (is.null(edge.selected)) edge.selected <- 1
    #   mat <- E(MFN@metabolic_network)[[edge.selected]]$atom_transfer
    #   mat.id <- rownames(mat@transfer_matrix)
    #   Atom_transfer_id_selected(mat.id[1])
    #   updateSelectInput(inputId = "Atom_transfer_id",
    #                     choices = mat.id,
    #                     selected = mat.id[1])

        ### reactiveVal update
        {
    #      rat_for_vis(mat)
    #      mat_for_update(mat)

        }
      })


    }

    output$Reaction_info <- renderText({

      message_with_time("Reaction_info")

      ### in
      {

        mfn <- mfn_for_update()
        rid <- mfn.selected.reaction()
        vda <- V(mfn@metabolic_network)[rid]

      }



      ### out
      {
        paste0("Reaction ID: ",vda$id,"\n",
               "Reaction Name: ",vda$Name,"\n",
               "Reaction Equation: ",vda$equation,"\n",
               "Reaction Enzyme: ",vda$enzyme,"\n"
        )
      }

    })




    ### update Molecule_igraph
    {

      mol.ig.selected <- reactiveVal(
        V(object@metabolic_network)[[1]]$Molecule_igraph
      )


    }

    output$Metabolite_isotopomer_statu_vis <- visNetwork::renderVisNetwork({

      message_with_time("Metabolite_isotopomer_statu_vis")
      {
        mol.ig <- mig_for_vis()
        iso.select <- selected.isotopomer()
      }

      vis_Molecule_igraph_isotopomer(mol.ig,isotopomer = iso.select)

    })


    output$Metabolite_isotopomer_statu_table <- DT::renderDataTable({

      message_with_time("Metabolite_isotopomer_statu_table")
      {
        mol.ig <- mig_for_vis()
      }


      dt <- mol.ig@isotopomer[,c("isotopomer","isotopologue","label","path"),drop = F]
      dt%>%
        DT::datatable(selection =list(mode = "single", selected = 1))
    })


    ### update selected.isotopomer() when select in Metabolite_isotopomer_statu_table
    {
      selected.isotopomer <- reactiveVal(1)
      observeEvent(input$Metabolite_isotopomer_statu_table_rows_selected,{
        selected.isotopomer(input$Metabolite_isotopomer_statu_table_rows_selected)
      })

    }



    ### save MFN to files
    {

      observeEvent(input$save_buttion,{

        message("save_buttion")

        save.dir <- "C:/Users/91879/OneDrive/Code/R/Projecct/2024.01.11.MSIP/Data/Metabolic_flux_network/"
        saveRDS(mfn_for_update(),
                file = paste0(save.dir,input$save_name,"_",str_time(),".rds"))


        message("MFN_saved")

      })
    }



  }
}
