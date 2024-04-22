MSIP_shiny_server <- function(iso.msip.list){


  function(input, output, session) {



    output$test_info <- renderPrint({


      #vdata(mol.ig())
    })


    ### reactiveval
    {
      iso.msip <- reactiveVal(iso.msip.list[[1]])
      iso.sp.data <- reactiveVal()
      fg.selected <- reactiveVal()
      fid.selected <- reactiveVal()
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
        iso.msip( iso.msip.list[[input$metabolite_table_rows_selected]] )
        updateSelectInput(inputId = "select_iso_count",
                          choices = names(iso.msip()$msip_data),
                          selected = names(iso.msip()$msip_data)[1])
        mol.ig.show.label(F)
      })



    }

    ### Spectra
    {
      #input$select_iso_count
      observe({

        iso.sp.data(shiny_get_sp_data(iso_msip = iso.msip(),
                                      iso_count =  input$select_iso_count))
        fg.selected(shiny_get_fg(sp.data = iso.sp.data(),
                                 x = event_data("plotly_click",
                                                source = "plotly_ms2_sp",
                                                priority  = "event")$x))
        })

      output$plotly_ms2_sp <-  renderPlotly({

        shiny_plotly_iso_msip_spectra(iso.sp.data())

      })




    }

    ### Vis struc
    {

      observeEvent(fg.selected(),{

        updateSelectInput(inputId = "select_fid",
                          choices = shiny_get_fid(iso_msip = iso.msip(),
                                                  fg = fg.selected()))

      })


      observe({
        mol.ig(iso.msip()$CFM_annotation@fragment_igraph[[1]])
        frag.ig(shiny_get_fg_ig(iso_msip =iso.msip(),
                                fid = input$select_fid  ))
        atom_map <- shiny_get_atom_map(iso.msip(),
                                       input$select_fid ,
                                       prob = T)
        mol.atom.map(atom_map[[1]])
        frag.atom.map(atom_map[[2]])
        mol.atom.c.prob(iso.msip()$msip_data[[input$select_iso_count]]$c.prob)
      })

      output$mol_graph <- renderVisNetwork(
        {

          shiny_vis_ig(ig = mol.ig() ,
                       prob.border = mol.atom.map(),
                       prob.fill =mol.atom.c.prob())
        }
      )
      output$frag_graph <- renderVisNetwork(
        {

          shiny_vis_ig(frag.ig(),
                       prob.border  = frag.atom.map() )

        }
      )
      output$frag_formula <- renderText({
        shiny_get_frag_formula(iso_msip = iso.msip(),
                               fid =input$select_fid  )
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

        if (!length(mol.atom.c.prob())) {
          tb <- NULL
        }else{
          tb <- data.frame(Atom = names(mol.atom.c.prob()),
                           Prob = mol.atom.c.prob())%>%
            dplyr::mutate(Prob = round(Prob,4))
        }
        tb


      })

      observeEvent(input$mol.ig.button,{
        mol.ig.show.label(!mol.ig.show.label() )
        #updateActionButton(inputId = "mol.ig.button",
        #                   label = ifelse(
        #  mol.ig.show.label(),
        #  "Show Atom ID",
        #  "Hide Atom ID"))
        igvd <- shiny_change_ig_label(ig = mol.ig(),
                              id = mol.ig.show.label())

        visNetworkProxy(shinyId = "mol_graph")%>%
          visUpdateNodes(nodes = igvd,
                         updateOptions= F)

      })
    }






  }

}
