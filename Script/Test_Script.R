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



cfm_data <- read_CFM_annotate_result()


lw <- 5
sw <- 3

case_when(width = case_when(bond_type== 1 ~ lw,
                            bond_type== 2  ) )

bond.data<- bond.data%>%
  dplyr::mutate(bond_id = num2str(1:n()))%>%
  dplyr::slice(rep(1:n(),bond_type) )%>%
  dplyr::group_by( bond_id )%>%
  dplyr::mutate( id = 1:n(),
                 width = id*5,
                 color = (id%%2 +(bond_type+1)%%2)%%2,
                 color = ifelse(color==0,"white","grey"))

nodes <- data.frame(id = 1:4)
edges <- data.frame(from = c(2,4,3,3,2), to = c(1,2,4,2,1))

visNetwork(nodes, edges, width = "100%")



"[O-]CC1OC(O)C(O)C(O)C1O"%>%
  smiles2sdf()%>%
  `[[`(1)%>%
  get_sdf_igraph()%>%
  vis_sdf_igraph(highlight = 1:3)



vis_sdf_igraph(graph.from)
vis_sdf_igraph(graph.to)



###
ggplot(fragment.data,
       aes(x = (dis.to.root),
           y =retrieve.ratio ,
           group =dis.to.root,
           col = dis.to.root))+
  labs(title = "Fragment Atom Tracealibity",
       x = "Distance to root",
       y = "Traciable atom ratio",
       col = "Distance to root")+
  geom_boxplot()+
  ggforce::geom_sina()+
  scale_color_gsea()+
  theme_bw()->p
open_ggplot_win(p,5,3)


a <- shortest_paths(frag.trans.graph,1)$vpath%>%
  sapply(length)



visIgraph(trans.graph)%>%
  visNodes()%>%
  visOptions(
    highlightNearest = T,
    nodesIdSelection = list(
      enabled=T,
      useLabels = F
    ),
    manipulation = TRUE)
# Thu Dec 14 15:34:46 2023 ------------------------------
cfm.data <- read_CFM_annotate_result()
f3 <- cfm.data$fragment_transition
f3$smile[!check_smile(f3$smile)]%>%table()
f3.sdf[[1]]

a <- read.SDFset("d:/temp/Structure2D_COMPOUND_CID_783.sdf"  )
b <- read.SDFset("D:/temp/Structure2D_COMPOUND_CID_962.sdf")
smile_map <- c(a,b)

smile_map@ID <- c("[HH]","O")
c[["H2"]]


a <- peak.assignment%>%
  dplyr::mutate(fragment_id = unname(fragment.id.new[fragment_id]))
a$fragment_id


a <- read_CFM_annotate_result()

cfm.data <- read_CFM_annotate_result()
frag.def <- cfm.data$fragment_define
frag.trans <- cfm.data$fragment_transition
frag.sdf <- get_smile_sdf(cfm.data$fragment_define$smile,
                          smile.id = frag.def$fragment_id)

i <- 13
frag.parent <- frag.sdf[[frag.trans$from[i]]]
frag.product <- frag.sdf[[frag.trans$to[i]]]
frag.fmc


sdf.igraph1 <- get_sdf_igraph(frag.parent)
sdf.igraph2 <- get_sdf_igraph(frag.product)


frag.fmc <- fmcs(frag.parent,frag.product,au = 0,bu = 5)
frag.fmc
i_fmc <- 1
vis_sdf_igraph_compare(sdf.igraph1 ,sdf.igraph2 ,
                       frag.fmc@mcs1$mcs1[[i_fmc]],
                       frag.fmc@mcs2$mcs2[[i_fmc]])


vis_sdf_igraph(sdf.igraph1)

sdf.igraph.highlight <- add_sdf_igraph_highlight(sdf.igraph1)
V(sdf.igraph.highlight)$shape = "icon"
V(sdf.igraph.highlight)$icon.face = 'Ionicons'
V(sdf.igraph.highlight)$icon.code  = "f1005"
sdf.igraph.highlight%>%
  visIgraph(idToLabel = T)%>%
  visNodes(font = list(size = 40,
                       strokeWidth = 10),
           size = 60,
           borderWidth  = 5,
           color = list(background = "transparent",
                        border = "#2B7CE9"))%>%
  visEdges(arrows = list(to = F),
           length = 0.8)%>%
  addIonicons()




CFM_annotate(smiles_or_inchi = smile,
             spectrum_file = paste0(project.dir,"/Result/Spectra.",id,".for.cfm.txt"),
             id = id,
             param_adduct = "[M-H]-",
             output_file =NULL) ->a



glu.cfm.pred <- CFM_predict(prob_thresh = 0.001,param_adduct = "[M-H]-")
glu.cfm.fraggen <- CFM_fraggen(max_depth = 2,param_adduct = "[M-H]-")
glu.cfm.exp <- read_CFM_annotate_result()


peak.pred <- glu.cfm.pred$peak_assignment
peak.exp <- glu.cfm.exp$peak_assignment
frag.def.pred <- glu.cfm.pred$fragment_define
frag.def.exp<- glu.cfm.exp$fragment_define
frag.def.gen <- glu.cfm.fraggen$fragment_define

frag.def.exp$smile %in% frag.def.gen$smile

# Sat Dec 16 16:39:48 2023 ------------------------------
frag.def.exp <- frag.def.exp%>%
  dplyr::mutate(formula = get_smile_formula(smile),
                mz = chemform_mz(formula,-1),
                inS = smile %in% frag.def.pred$smile)

frag.def.pred <- frag.def.pred%>%
  dplyr::mutate(formula = get_smile_formula(smile),
                mz = chemform_mz(formula,-1),
                inS = smile %in% frag.def.exp$smile)

peak.pred <- peak.pred%>%
  dplyr::mutate(mz = frag.def.pred$mz[match(fragment_id,frag.def.pred$fragment_id)],
                smile = frag.def.pred$smile[match(fragment_id,frag.def.pred$fragment_id)])%>%
  dplyr::group_by(energy)%>%
  dplyr::mutate(esum = sum(fragment_score,na.rm =T))%>%
  dplyr::group_by(energy,mz)%>%
  dplyr::mutate(msum = sum(fragment_score))%>%
  dplyr::group_by(energy)%>%
  dplyr::mutate(nmsum = msum/max(msum,na.rm = T)*100,
                peak_ratio = 100*intensity/(sum(unique(intensity))))

peak.exp <- peak.exp%>%
  dplyr::mutate(mz = frag.def.exp$mz[match(fragment_id,frag.def.exp$fragment_id)],
                smile = frag.def.exp$smile[match(fragment_id,frag.def.exp$fragment_id)])%>%
  dplyr::group_by(energy)%>%
  dplyr::mutate(esum = sum(fragment_score,na.rm =T))%>%
  dplyr::group_by(energy,mz)%>%
  dplyr::mutate(msum = sum(fragment_score))%>%
  dplyr::group_by(energy)%>%
  dplyr::mutate(nmsum = msum/max(msum,na.rm = T)*100,
                peak_ratio = 100*intensity/(sum(unique(intensity))))




MSdev.tca <- MSdev("d:/2023.11.MSIP/20231213_TCA/Result/")
MSdev.tca <- MSdev_msConvert(MSdev.tca)
MSdev.tca <- MSdev_checkSampleInfo(MSdev.tca)
MSdev.tca <- MSdev_xcmsProcessing(MSdev.tca)
saveMSdev(MSdev.tca)


xcms.xcms <- MSdev.tca@xcmsData$NegativeMS1
xcms.sp.TCA <- get_xcms_Spectra(xcms.xcms)
xcms.ms2.TCA <- filterMsLevel(xcms.sp.TCA,2)
xcms.ms2.data.TCA <- spectraData(xcms.ms2.TCA)%>%
  as.data.frame( )%>%
  dplyr::mutate(dp = precursorIntensity * injectionTime/1000)


ggplot(xcms.ms2.data)+
  geom_point(aes(x = log10(dp),
                 y = log10(totIonCurrent),
               #  y= injectionTime,
                 col = collisionEnergy))+
  geom_abline(slope = 1)+
  scale_color_gsea()


# Sat Dec 16 20:08:28 2023 ------------------------------
msdev <- load_as_var("d:/2023.11.MSIP/231108_Glc_Tracing_rawdata/MSdev_2023_11_09.Rdata")
xcms.xcms <- msdev@xcmsData$PositiveMS1
xcms.sp <- get_xcms_Spectra(xcms.xcms)
xcms.ms2 <- filterMsLevel(xcms.sp,2)
xcms.ms2.data <- spectraData(xcms.ms2)%>%
  as.data.frame()%>%
  dplyr::mutate(dp = precursorIntensity * injectionTime/1000)

xcms.ms2.data <-xcms.ms2.data[sample(1:nrow(xcms.ms2.data),
                                     3000),]

ggplot(xcms.ms2.data)+
  geom_point(aes(
                  #x = injectionTime,
                 x = log10(dp),
                 y = log10(basePeakIntensity),
                 #  y= injectionTime,
                 col = injectionTime))+
  geom_abline(slope = 1)+
  scale_color_gsea()



# Sun Dec 17 22:14:58 2023 ------------------------------
msdev <- load_as_var("d:/2023.11.MSIP/231108_Glc_Tracing_rawdata/MSdev_2023_11_09.Rdata")
xcms.pos <- msdev@xcmsData$PositiveMS1
xcms.neg <- msdev@xcmsData$NegativeMS1
xcms.pos <- filterFile(xcms.pos,which(!grepl("Blank",sampleNames(xcms.pos))))
xcms.neg <- filterFile(xcms.neg,which(!grepl("Blank",sampleNames(xcms.neg))))

xcms.tune.pos <- get_xcms_Autotuner(xcms.pos)
xcms.tune.neg <- get_xcms_Autotuner(xcms.neg)



msdev.tca <- MSdev("d:/2023.11.MSIP/231108_Glc_Tracing_rawdata/rawData/")
msdev.tca <- MSdev_msConvert(msdev.tca)
msdev.tca <- MSdev_checkSampleInfo(msdev.tca)
msdev.tca@experimentInfo <- MS_Experiment[10]
get_MSdev_param(msdev.tca)
msdev.tca <- MSdev_xcmsProcessing(msdev.tca)


# Mon Dec 18 14:39:49 2023 ------------------------------
msdev.tca <- load_as_var("d:/2023.11.MSIP/20231126Tune/MSdev_2023_11_27.Rdata")
get_MSdev_param(msdev.tca)
xcms.xcms<- msdev.tca@xcmsData$NegativeMS1
xcms.peaks <- get_xcms_peaks_stat(xcms.xcms)
table(sampleNames(xcms.xcms)[xcms.peaks$sample])



# Tue Dec 19 09:24:18 2023 ------------------------------
i <- 7
expSpec <- sp.exp[[i]]
refSpec <- sp.ref[[i]]
scorem <- compareSpectra(expSpec,refSpec)
scorem[is.infinite(scorem)|is.na(scorem )] <- 0
dim(scorem) <- c(length(expSpec),length(refSpec))

apply(scorem,2,max,na.rm=T)


Spectra_database[13997]%>%
  filterSpectra_below_PrecursorMz()


xcms.xcms <- load_demo("xcms")
a <- chromPeaks(xcms.xcms)

groupf

groupFeatures()

xmse <- groupFeatures(xcms.xcms, EicSimilarityParam(threshold = 0.7, n = 2))
xmse <- xcms.xcms
plot(featureDefinitions(xmse)$rtmed, featureDefinitions(xmse)$mzmed,
     xlab = "retention time", ylab = "m/z", main = "features",
     col = "#00000080", pch = 21, bg = "#00000040")
grid()
xmse <- groupFeatures(xmse, param = SimilarRtimeParam(5))
#plotFeatureGroups(xmse)
table(featureGroups(xmse))
xmse <- groupFeatures(
  xmse, AbundanceSimilarityParam(threshold = 0.7, transform = log2),
  filled = TRUE)
table(featureGroups(xmse))


cor_plot <- function(x, y) {
  C <- cor(x, y, use = "pairwise.complete.obs")
  col <- ifelse(C >= 0.7, yes = "#0000ff80", no = "#ff000080")
  points(x, y, pch = 16, col = col)
  grid()
}
fts <- grep("FG.0409.001", featureGroups(xmse))
pairs(t(fvals[fts, ]), gap = 0.1, main = "FG.040", panel = cor_plot)
register(SerialParam())
xmse <- groupFeatures(xmse, EicSimilarityParam(threshold = 0.7, n = 5))
table(featureGroups(xmse))

fidx <- grep("FG.013.001.", featureGroups(xmse))
eics <- featureChromatograms(
  xmse, features = rownames(featureDefinitions(xmse))[fidx],
  filled = TRUE, n = 1)

xcms.grouped <- xcms_get_feature_group(xcms.xcms)




load("d:/temp/xcms.feature.group.temp")


msdev <- load_as_var("d:/2023.11.MSIP/20231126Tune/MSdev_2023_11_27.Rdata")
xcms.xcms <- msdev@xcmsData$PositiveMS1
xcms.grouped <- xcms_get_feature_group(xcms.xcms)



a <- featureDefinitions(xcms.grouped)%>%
  as.data.frame()


get_xcms_feature_group_split <- function(feature_group ){

  a <- str_split(pattern  = "\\.",string =  feature_group)
  b <- do.call("rbind",a)%>%
    as.data.frame()
  table(b$V4)


}


msdev.demo <- MSdev("d:/2023_12_19-Xinchenhao/Data/")
msdev.demo <- MSdev_checkSampleInfo(msdev.demo)
msdev.demo <- MSdev_msConvert(msdev.demo)
msdev.demo <- MSdev_xcmsProcessing(msdev.demo)
msdev.demo <- MSdev_extract_Spectra(msdev.demo)
msdev.demo <- MSdev_match_Spectra_to_feature(msdev.demo)
msdev.demo <- MSdev_annotation(msdev.demo,
                               db.path = "d:/MSdb.2023.05.30/LipidBlast.rda")
msdev.demo <- MSdev_get_Stat(msdev.demo)
saveMSdev(msdev.demo)
MSdev_export(msdev.demo)



xcms.metabolite <- msdev.demo@statData$metabolite.se%>%
  rowData()%>%
  as.data.frame()



# Sat Dec 23 14:21:58 2023 ------------------------------
get_MSdev_param(msdev.demo)



# Mon Dec 25 16:56:17 2023 ------------------------------
ppm.record <- data.frame(ppm = 1:10000,c13.count = 0)
for (i in 1:10000) {

  message(i)
  fdf.edge <- fdf.edge%>%
    rowwise()%>%
    dplyr::mutate(closest.c13.count = which.min(abs(c13.mz-mz.diff)))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(closest.c13.mz = c13.mz[closest.c13.count],
                  mz.error = abs(mz.diff-closest.c13.mz),
                  is.c13 = mz.error/from.mz < i*1e-6)
  ppm.record$c13.count[i] <- sum(fdf.edge$is.c13)
}

plot_ly(ppm.record,x = ~ppm, y = ~c13.count)


CentWavePredIsoParam()
MassifquantParam()







xcms.mf <- filterFile(xcms.xcms,1)
xcms.mf <- findChromPeaks(xcms.mf,param = MassifquantParam())
xcms.peaks <- chromPeaks(xcms.mf)

a <- get_xcms_peaks_stat(xcms.xcms )
b <- a%>%
  dplyr::mutate(into_maxo = into/maxo,
                i_m_t = into_maxo/peakWidth)%>%
  dplyr::slice(sample(1:n(),3000))

plot(b$into_maxo,b$peakWidth)
plot(b$i_m_t)



# Wed Dec 27 08:55:53 2023 ------------------------------
msdev.dcx <- MSdev("d:/2023.12.26.DCX.Metabolomic/Data/")
msdev.dcx <- MSdev_msConvert(msdev.dcx)
msdev.dcx <- MSdev_checkSampleInfo(msdev.dcx)
msdev.dcx <- MSdev_xcmsProcessing(msdev.dcx)
msdev.dcx <- MSdev_extract_Spectra(msdev.dcx)
msdev.dcx <- MSdev_match_Spectra_to_feature(msdev.dcx)
msdev.dcx <- MSdev_annotation(msdev.dcx,
                              expand_adduct= T,
                              db.path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/SpectraDB/HMDB_KEGG_export_human_pathway.rda")
msdev.dcx <- MSdev_get_Stat(msdev.dcx)
saveMSdev(msdev.dcx)
MSdev_export(msdev.dcx)



groupSimilarityMatrix(x, threshold = 0.9)

# Wed Dec 27 19:10:39 2023 ------------------------------
msdev.dcx <- load_as_var("d:/2023.12.26.DCX.Metabolomic/MSdev_2023_12_27.Rdata")
msdev.dcx <- MSdev_checkSampleInfo(msdev.dcx)
c("[M-H]-",
  "[M-H2O-H]-",
  "[2M-H]-" ,
  "[M+FA-H]-" ,
  "[M+H]+" ,
  "[M-H2O+H]+",
  "[M+NH4]+"
)

# Wed Dec 27 19:23:18 2023 ------------------------------
data.se <- get_MSdev_DEP_se(msdev.dcx,"metabolite")
data.se <- data.se[,!data.se$group%in% c("QC","Blank")]
rowData(data.se)$name <- rowData(data.se)$ID
data_filt <- DEP::filter_missval(data.se,
                                 thr = max(table(data.se$group))*0)
data_norm <- DEP::normalize_vsn(data_filt)
plot_normalization(data_filt,data_norm)
data_imp <- DEP::impute(data_norm, fun = "MinProb")
#plot_missval(data_norm)
#p.pca <- DEP_plot_PCA(data.se)
data.diff <- DEP_test_diff(data_imp)
data.diff <- DEP_add_rejections(data.diff,p.adjust = T)
DEP_plot_volcano(data.diff,
                 contrast = DEP_list_contrast(data.diff)[4],
                 show.label = T)
p <- DEP_plot_volcano(data.diff,contrast = "all")



# Thu Dec 28 12:13:08 2023 ------------------------------
x <- rbind(
  c(12, 34, 231, 234, 9, 5, 7),
  c(900, 900, 800, 10, 12, 9, 4),
  c(25, 70, 400, 409, 15, 8, 4),
  c(12, 13, 14, 15, 16, 17, 18),
  c(14, 36, 240, 239, 12, 7, 8),
  c(100, 103, 80, 2, 3, 1, 1)
)
x.cor <- cor(t(x))
res <- groupFeatures(x, AbundanceSimilarityParam())
res



# Thu Dec 28 12:21:15 2023 ------------------------------
x <- c(2, 3, 4, 5, 10, 11, 12, 14, 15)

## Group the values using a `group` function. This will create larger
## groups.
groupFeatures(x, param = SimilarRtimeParam(2, MsCoreUtils::group))

## Group values using the default `groupClosest` function. This creates
## smaller groups in which all elements have a difference smaller than the
## defined `diffRt` with each other.
groupFeatures(x, param = SimilarRtimeParam(5))

## Grouping on a SummarizedExperiment
##
## load the test SummarizedExperiment object
library(SummarizedExperiment)
data(se)

## No feature groups defined yet
featureGroups(se)

## Determine the column that contains retention times
rowData(se)

## Column "rtmed" contains the (median) retention time for each feature
## Group features that are eluting within 10 seconds
res <- groupFeatures(se, SimilarRtimeParam(10), rtime = "rtmed")

featureGroups(res)

## Evaluating differences between retention times within each feature group
rts <- split(rowData(res)$rtmed, featureGroups(res))
lapply(rts, function(z) abs(diff(z)) <= 10)

## One feature group ("FG.053") has elements with a difference larger 10:
rts[["FG.053"]]
abs(diff(rts[["FG.053"]]))

## But the difference between the **sorted** retention times is < 10:
abs(diff(sort(rts[["FG.053"]])))

## Feature grouping with pre-defined feature groups: groupFeatures will
## sub-group the pre-defined feature groups, features with the feature group
## being `NA` are skipped. Below we perform the feature grouping only on
## features 40 to 70
fgs <- rep(NA, nrow(rowData(se)))
fgs[40:70] <- "FG"
featureGroups(se) <- fgs
res <- groupFeatures(se, SimilarRtimeParam(10), rtime = "rtmed")
featureGroups(res)




p+ggrepel::geom_text_repel(data = label_df,
            aes(x= x,y=y,label = label),
            hjust = 0)->x


# Fri Dec 29 18:49:42 2023 ------------------------------
MSdev_annotation()
MSdev_get_Stat()


dn <- density(1:100,bw= 5  )
plot(dn)


plotChromPeakDensity(xcms.xcms,simulate=T,
                     param = PeakDensityParam(pData(xcms.xcms)$group,
                                              bw = 200,
                                              minFraction = 0.5,
                                              binSize = 0.001))


# Sat Dec 30 23:36:58 2023 ------------------------------
msdev.ylf <- load_as_var("D:/2023.11.MSIP/20231221_FS/MSdev_2023_12_23.Rdata")
xcms.xcms <- msdev.ylf@xcmsData$PositiveMS1

xcms.xcms <- groupChromPeaks(xcms.xcms,param =
                               PeakDensityParam(pData(xcms.xcms)$group,
                                                bw = 200,
                                                minFraction = 0.7,
                                                binSize = 0.002))
xcms.xcms <- xcms_get_feature_stat(xcms.xcms)
fdf <- featureDefinitions(xcms.xcms)%>%
  as.data.frame()
ggplot(fdf)+
  geom_density(aes(x = (mzmax-mzmin)/mzmed*1e6))+
  #geom_density(aes(x = Glu))+
  xlim(c(0,5))



# Sun Dec 31 00:06:16 2023 ------------------------------
cox.data%>%
  dplyr::group_by(phase)%>%
  summarise(mean(flare.time))


