# Thu Nov 30 18:49:39 2023 ------------------------------
msdev.xch <- MSdev("d:/2023.11.30.XCH/rawData/")
msdev.xch <- load_as_var("d:/2023.11.30.XCH/MSdev_2023_11_30.Rdata")
msdev.xch <- MSdev_msConvert(msdev.xch)
msdev.xch <- MSdev_checkSampleInfo(msdev.xch)
msdev.xch <- MSdev_xcmsProcessing(msdev.xch)
msdev.xch <- MSdev_extract_Spectra(msdev.xch)
msdev.xch <- MSdev_match_Spectra_to_feature(msdev.xch)
msdev.xch <- MSdev_annotation(msdev.xch,
                              db.path = "D:/MSdb.2023.05.30/LipidBlast.rda")
msdev.xch <- MSdev_get_Stat(msdev.xch)
saveMSdev(msdev.xch)
exportMSdev(msdev.xch)

# Fri Dec  1 14:31:37 2023 ------------------------------
install.packages("d:/temp/ChemmineOB_1.40.0.tar.gz",
                 repo = NULL,
                 configure.args = list(OPEN_BABEL_INCDIR="d"))

i <- 4
plot(fragment.sdf[i])
conMA(fragment.sdf[i])

p <- ggplot_sdf(sdf)
g <- ggplotGrob(p)
p.patch <- p.sp+annotation_custom(g,xmin = 500,900,
                       ymin = 1e6,ymax = 1.5e6)
open_ggplot_win(p.patch,5,4)


###
ggdraw()+
  draw_plot(p.sp,0,0,1,1)+
  draw_plot(p,0.8,0.8,.1,.1)->p.patch


f <- function(smile){
  suppressMessages(f_smile_sdf(smile))
}

f_smile_sdf <- function (smiles)
{
  if (!any(class(smiles) %in% c("character", "SMIset"))) {
    stop("input must be SMILES strings stored as \"SMIset\" or \"character\" object")
  }
  if (inherits(smiles, "SMIset"))
    smiles <- as.character(smiles)
  .ensureOB()
  sdf = ChemmineR:::definition2SDFset(convertFormat("SMI", "SDF", paste(paste(smiles,
                                                                  names(smiles), sep = "\t"), collapse = "\n")))
  #cid(sdf) = sdfid(sdf)
  sdf
}




# Mon Dec  4 14:17:07 2023 ------------------------------
library(igraph)
g <- make_graph(~ Alice-Bob:Claire:Frank, Claire-Alice:Dennis:Frank:Esther,
                George-Dennis:Frank, Dennis-Esther)

plot(g,layout = layout_randomly(g))

# Wed Dec  6 09:30:46 2023 ------------------------------
nodes <- data.frame(id = 1:3)
edges <- data.frame(from = c(1,2), to = c(1,3))
visNetwork(nodes, edges) %>%
  visNodes(shape = "square",
           title = "I'm a node", borderWidth = 3)
visNetwork(nodes, edges) %>%
  visNodes(color = list(hover = "green")) %>%
  visInteraction(hover = TRUE)

visNetwork(nodes, edges) %>% visNodes(color = "red")








# Wed Dec  6 12:53:13 2023 ------------------------------
cfm.data <- read_CFM_annotate_result()




temp.path <- "d:/temp/"
node.data <- cfm.data$fragment_define2%>%
  dplyr::mutate(id = fragment_id,
                name = fragment_id,
                label = smile,
                checked = check_smile(smile),
                formula = get_smile_formula(smile),
                image =  paste0(temp.path,"Node",id,".png"),
                shape = "image")%>%
  column_to_rownames("fragment_id")
edge.data <- cfm.data$fragment_transition%>%
  dplyr::mutate(arrows = "to",
                formula = get_smile_formula(smile),
                formula.mz = chemform_mz(formula))%>%
  dplyr::mutate(from.mz = node.data[`from`,"fragment_mz"],
                to.mz = node.data[`to`,"fragment_mz"],
                diff.mz = from.mz-to.mz)


for (i in 1:nrow(node.data)) {

  if (node.data$checked[i]) {

    p <- ggplot_sdf(smiles2sdf(node.data$smile[i])[[1]],show_ele = T)
    export::graph2png(p,
      file = node.data$image[i],
      width = 2,height = 2
    )
  }

}


visNetwork(node.data[1:10,], edge.data) %>%
  visNodes(shapeProperties = list(useBorderWithImage = TRUE),
           image = "file:///D:/temp/Node0.png") %>%
  visLayout(randomSeed = 2)




ig <- graph_from_data_frame(edge.data,vertices = node.data)


# Wed Dec  6 17:59:21 2023 ------------------------------


cfm.data <-read_CFM_annotate_result()

#name.map <- paste0("nd",1:nrow(node.data))
#names(name.map) <- node.data$id
#$emp.path <- "d:/temp/"

node.data <- cfm.data$fragment_define2%>%
  dplyr::mutate(id = fragment_id,
                #name = name.map[id],
                checked = check_smile(smile),
                formula = get_smile_formula(smile),
                label = formula,
                #image =  paste0("/pic/Node",id,".png"),
                #image.name =  paste0("Node",id,".png"),
                shape = "image")
edge.data <- cfm.data$fragment_transition%>%
  dplyr::mutate(arrows = "to",
                formula = get_smile_formula(smile),
                formula.mz = chemform_mz(formula))%>%
  dplyr::mutate(from.mz = node.data[`from`,"fragment_mz"],
                to.mz = node.data[`to`,"fragment_mz"],
                diff.mz = from.mz-to.mz
                #from=name.map[from],
               # to = name.map[to]
               )

visNetwork(node.data[1:10,],edge.data)%>%
  visNodes(image = "aaa:///Node0.png")




load_all()
shinyApp(ui= ui,server = get_server(cfm.data ))
# Thu Dec  7 14:32:55 2023 ------------------------------


msdev.agc <- MSdev("d:/2023.11.MSIP/20231205_MS2_PARAM/Result/")
msdev.agc <- MSdev_msConvert(msdev.agc)
msdev.agc <- MSdev_checkSampleInfo(msdev.agc)
msdev.agc <- MSdev_xcmsProcessing(msdev.agc)
saveMSdev(msdev.agc)


msdev.agc <- load_as_var("d:/2023.11.MSIP/20231205_MS2_PARAM/MSdev_2023_12_08.Rdata")
xcms.xcms <- msdev.agc@xcmsData$NegativeMS1
xcms.xcms <- xcms_get_dda_ms2_assignment(xcms.xcms)
plot_xcms_dda_acquisition(xcms.xcms)
xcms.scam <- get_xcms_scan_Stat(xcms.xcms )%>%
  dplyr::filter(!is.na(precursorMZ))%>%
  dplyr::mutate(groupMz(precursorMZ))
table(sampleNames(xcms.xcms)[xcms.scam$fileIdx],xcms.scam$mz.center)

sp <- get_xcms_Spectra(xcms.xcms)
sp <- filterMsLevel(sp,2)
plot_Spectra_Injection(sp)



files <- unique(sp.data$dataOrigin)
sp.data <- spectraData(sp)%>%as.data.frame()%>%
  dplyr::mutate(condition = basename(dataOrigin),
                condition = factor(condition,level = unique(condition)))%>%
  dplyr::filter(dataOrigin %in% files[c(5)])
sp.data$injectionTime
ggplot(sp.data)+
  geom_point(aes(x = log10(precursorIntensity),
                 y = log10(totIonCurrent),
                 col = injectionTime),
             shape = 19,
             size = 2,
             stroke = 0,
             #col = "transparent",
             alpha = 0.5)+
  #ggsci::scale_fill_gsea()+
  ggsci::scale_color_gsea()+
  theme_bw()



ggplot(sp.data)+
  geom_jitter(aes(x = condition,
                 y = totIonCurrent))+
  scale_y_log10()


xcms.xcms <- load_demo("xcms")
xcms.xcms <- xcms_get_dda_ms2_assignment(xcms.xcms)
plot_xcms_dda_acquisition(xcms.xcms)->p
open_ggplot_win(p,16,9)





