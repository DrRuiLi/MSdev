
x <- MS.network[[5]]
##
adduct.candidate <- x[["adduct.candidate"]]
adduct.mz <- adduct.candidate$exact.mz
peak.mz <- xcms.peaks[,"mz"]
adduct.matrix <- matrix(rep(adduct.mz,length(peak.mz)) , nrow = length(adduct.mz))
peak.matrix <- matrix(rep(peak.mz,length(adduct.mz)) ,
                      nrow = length(adduct.mz),
                      byrow = T)
sub.matrix <- adduct.matrix - peak.matrix
tol.matrix <- peak.matrix * ppm.thresh*1e-6
pass.matrix <- abs(sub.matrix) < tol.matrix

matched.id <- which(pass.matrix,arr.ind = T)
matched.id

adduct <- cbind(adduct.candidate[matched.id[,1],],
                matrix(xcms.peaks[matched.id[,2],],nrow = 1))


b <-MS.network.pos[[16]][["adduct"]]
ggplot(b)+
  geom_point(aes(x = peak.rt,y=peak.mz,size = log10(peak.intb),col = adduct),
             alpha = 0.8)+
  ggsci::scale_color_aaas()+
  xlim(c(0,750))
export::graph2png(file = "b.png",width = 10,height = 8)

ggplot(as.data.frame(xcms.pos.peaks))+
  geom_point(aes(x = rt,y=mz,
                 col = log10(intb),
                 alpha = log10(intb)/10,
                 size = rtmax-rtmin),
             )+
  scale_color_gradient2(low = "yellow",mid = "red",high = "blue",midpoint = 5)
export::graph2png(file = "a.png",width = 5,height = 4)

ggplot(as.data.frame(xcms.peaks))+
  geom_point(aes(x = rt,y=mz,
                 col = log10(intb),
                 alpha = log10(intb)/10,
                 size = (rtmax-rtmin)),
  )+
  scale_size_area(max_size = 10)+
  xlim(c(0,800))+
  scale_color_gradient2(low = "yellow",mid = "red",high = "blue",
                        midpoint = 5)
export::graph2png(file = "b.png",width = 10,height = 8)

ggplot(xcms.peaks)+
  geom_segment(aes(x = rtmin,
                   xend = rtmax, y=mzmin,yend = mzmax,
                 col = log10(intb),
                 alpha = log10(intb)/10)
  )+
  scale_size_area(max_size = 20)+
  xlim(c(0,800))+
  scale_color_gradient2(low = "yellow",mid = "red",high = "blue",
                        midpoint = 5)
export::graph2png(file = "b.png",width = 5,height = 4)


spec.neg <- Spectra(unique(compound.record$mzML.negative))

a <-spec.neg%>%
  filterMsLevel(2)%>%
  filterPrecursorMzRange(c(129.00,129.02))
precursorMz(a)
rtime(a)
plotSpectra(a[1])

export_sample_information_from_wiff(raw.data.dir = choose.dir(),project.dir = choose.dir())

peak.density.param <- PeakDensityParam(sampleGroups = sample.info.xcms$sample.type,
                                       minFraction = 0.8,bw = 30,
                                       binSize = 0.015)

xcms.xcms <- groupChromPeaks(xcms.xcms,param = peak.density.param)


openxlsx::write.xlsx(compound.record , file = compound.record.file,
                     sheetName = "aaa")


a <- loadWorkbook(compound.record.file)
sheets(a)
compound.record <- read.xlsx(a)%>%
  select(pubchem.cid)
addWorksheet(wb = a,sheetName = "b")
sheets(a)
writeData(x = compound.record,wb = a,sheet = "b")
saveWorkbook(a , "a.xlsx")

a <- quantify(ms.ana$xcms.positive)%>%
  colData()%>%
  as.data.frame()
b <- featureChromatograms(ms.ana$xcms.positive,expandRt = 700,features = sample(1:10000,10))
plot(b[10],col  = brewer.pal(7,"Set1"))

a <- featureValues(xcms.xcms)%>%as.data.frame()

quantify(xcms.xcms)


xcms.xcms <- xcms.pos
sampleNames(xcms.xcms)
featureValues(xcms.xcms)%>%colnames()


sampleNames(xcms.xcms) <- paste0("S",1:7)
featureValues(xcms.xcms)%>%colnames()

xcms.xcms@phenoData@data$sampleNames<- paste0("S",1:7)
featureValues(xcms.xcms)%>%colnames()


fileNames(xcms.xcms) <-  paste0("S",1:7)

xcms.xcms@processingData@files <- paste0("S",1:7)
fileNames(xcms.xcms)
featureValues(xcms.xcms)%>%colnames()



xcms.xcms <- MsFeatures::groupFeatures(xcms.xcms,
                                       param = MsFeatures::SimilarRtimeParam(10))
plotFeatureGroups(xcms.xcms)
xcms.xcms <- MsFeatures::groupFeatures(xcms.xcms,
                                       param = MsFeatures::AbundanceSimilarityParam(threshold = 0.7,
                                                                                    transform = log2))
library(SummarizedExperiment)
se <- quantify(xcms.pos)
assay(se) <-featureValues(xcms.pos,filled = T,
                          missing = "rowmin_half")
assay(se)
assayNames(se)

tl <- c(1,system.time(a))
for (i in seq(1,10000,200)) {

  system.time(a <- groupFeatures(se[1:i],
                              param = SimilarRtimeParam(10), rtime = "rtmed")

  )->t
  message(i,", time: ",t[1])
  tl<-rbind(tl,c(i,t))
  x <- tl[,1]

  y <- tl[,2]
  #fit2 <- lm(y~poly(x,2))
  plot(x,y)
  #lines(x,predict(fit2))
}
plot(tl[,1],tl[,2])

x <-20000
5e-6 *x*x - 0.0009*x+0.1339


a <-data.frame(x = tl[,1],
               y = tl[,2])


fit <- lm(y~poly(x) , data = a)
fit

b <- data.frame(x = 1:10000 )
c <- predict(fit,b)
head(c)


std.matrix <- featureValues(xcms.pos,missing = "rowmin_half")[,1:5]

plot(std.matrix[6,])
std.matrix <- log2(std.matrix)
cor.list <- cor(t(std.matrix), compound.sample.info$concentration[1:5])
densityplot(cor.list)
sum(cor.list>0.5,na.rm = T)
which(cor.list >0.8)
plot(compound.sample.info$concentration[1:5],std.matrix[2570 ,])











xcms.xcms <- xcms.pos
featureGroups(xcms.xcms) <- NA
xcms.xcms <- MsFeatures::groupFeatures(xcms.xcms,
                                       param = MsFeatures::SimilarRtimeParam(10))
featureGroups(xcms.xcms)
plotFeatureGroups(xcms.xcms)

xcms.xcms <- MsFeatures::groupFeatures(xcms.xcms,
                                       param = MsFeatures::AbundanceSimilarityParam(threshold = 0.7))
featureGroups(xcms.xcms)%>%table()

t.s <- Sys.time()
xcms.xcms <- MsFeatures::groupFeatures(xcms.xcms,
                                       param = xcms::EicSimilarityParam(n = 2))
t.e <- Sys.time()
save(xcms.xcms,file = paste0("xcms.msfeature.temp.Rdata"))


xcms.peaks <- chromPeaks(xcms.pos)%>%as.data.frame()
table(xcms.peaks$sample)

xcms.to.plot <- filterFile(xcms.pos,7)
sampleNames(xcms.to.plot)
plot_xcms_peaks_distribution(xcms.to.plot)

library(SummarizedExperiment)
ms.features <- quantify(xcms.pos ,missing = "rowmin_half")
colData(ms.features)
ms.feature.matrix<- assay(ms.features)

a <- sample(nrow(ms.feature.matrix),100)
toplot <- ms.feature.matrix[a,]%>%max_min_normalize()


ComplexHeatmap::Heatmap(toplot,cluster_columns = F)



ms_condition_to_workbook <- function(ms_condition){

  ms_condition$column[[1]] -> ms_column
  ms_condition$phaseA[[1]] -> ms_phase_A
  ms_condition$phaseB [[1]]-> ms_phase_B
  ms_condition$gradient [[1]]-> ms_gradient

  ms_condition <- ms_condition[-which(colnames(ms_condition) %in% c("column","phaseA","phaseB","gradient"))]

  ms_condition_table <- as.data.frame(ms_condition,row.names = "value")%>%
    t%>%as.data.frame()%>%rownames_to_column("params")

  ### edit in excel
  temp.xlsx <- paste0(tempdir(), "/temp.xlsx")
  openxlsx::write.xlsx(
    list(
      "Exp_Condition" = ms_condition_table,
      "Column" = ms_column,
      "phase A" = ms_phase_A,
      "phase B" = ms_phase_B,
      "Gradient" = ms_gradient
    ),
    file = temp.xlsx
  )
  openxlsx::openXL(temp.xlsx)
  readline("Edit MS condition")
  ms_workbook <- openxlsx::loadWorkbook(file = temp.xlsx)

  ms_column <- openxlsx::read.xlsx(ms_workbook, sheet = "Column")
  ms_phase_A <- openxlsx::read.xlsx(ms_workbook, sheet = "phase A")
  ms_phase_B <- openxlsx::read.xlsx(ms_workbook, sheet = "phase B")
  ms_gradient <-
    openxlsx::read.xlsx(ms_workbook, sheet = "Gradient")


  ### load new data
  {
    ms_condition_table <-
      openxlsx::read.xlsx(ms_workbook, sheet = "Exp_Condition")
    ms_column <- openxlsx::read.xlsx(ms_workbook, sheet = "Column")
    ms_phase_A <- openxlsx::read.xlsx(ms_workbook, sheet = "phase A")
    ms_phase_B <- openxlsx::read.xlsx(ms_workbook, sheet = "phase B")
    ms_gradient <-
      openxlsx::read.xlsx(ms_workbook, sheet = "Gradient")


    ms_condition <- DataFrame(MSC_id = NA)
    for (i in 1:nrow(ms_condition_table)) {
      exprss <-
        paste0(
          "ms_condition$",
          ms_condition_table$params[i],
          " <- \"",
          ms_condition_table$values[i],
          "\""
        )
      eval(parse(text = exprss))
    }
    ms_condition$column <- ms_column %>% list()
    ms_condition$phaseA <- ms_phase_A %>% list()
    ms_condition$phaseB <- ms_phase_B %>% list()
    ms_condition$gradient <- ms_gradient %>% list()
    }
  return(ms_condition)
}

check_chemform(isotopes, chemforms = "[13]C0[15]N0[17]O0[18]O0[2]H5C9H6N1O2")



#49   63
data("chemforms")
for (i in 1:length(chemforms)) {
  formula <- chemforms[i]
  isopat <- isotopes_pattern_enviPat(formula )
  check <- check_chemform(isotopes , isopat$formula)
  if (any(check$warning)) {
    message(i,":",formula,"\n")
  }
}

enviPat::check_chemform(isotopes , "Br1")

export::graph2ppt(file = "temp.xlsx",width = 6,height = 4,append = T)



ms.features <- as.data.frame(xcms.pos.peaks)%>%
  filter( mz > 374 ,mz < 375)

ms.chrom <- chromatogram(xcms.xcms ,mz = c( 280.1605
,280.1621))

plot(ms.chrom)
grid()


peaks.ppm <- data.frame(xcms.peaks)%>%
  mutate(error = (mzmax-mz)/mz * 1e6)


xcms.xcms <- groupChromPeaks(xcms.xcms , param =PeakDensityParam())

xcms.feature <- featureDefinitions(xcms.xcms) %>%as.data.frame()%>%
  mutate(rtdiff = rtmax - rtmin,
         mzdiff = (mzmax-mzmin)/mzmed)

xcms.peaks <- chromPeaks(xcms.xcms)
a <- xcms.peaks[xcms.feature$peakidx[308][[1]],]




adduct <- MS.network[["2"]][["adduct"]]








#####################


extract_featuregroup <- function(i,MS.network){
  x <- MS.network[[i]]
  compound <- x[["compound"]]
  adduct <- x[["adduct"]]%>%
    dplyr::filter(rt.filter)%>%
    mutate(feature_group_id = paste0("FG",sprintf("%03d",i)))

  adduct

}
feature_group <- lapply(1:5, extract_featuregroup,MS.network)%>%
  data.table::rbindlist()%>%
  distinct(feature.id,feature_group_id )

xcms.features[feature_group$feature.id , "feature_group"] <-feature_group$feature_group_id

featureDefinitions(xdata)$feature_group <-xcms.features$feature_group

xcms.eic <- groupFeatures(xdata, EicSimilarityParam(threshold = 0.7, n = 1,onlyPeak = T))


ms.chrom.1 <- featureChromatograms(xcms.xcms , features = "FT0526")
ms.chrom.2 <- featureChromatograms(xcms.xcms , features = "FT0520")

compareChromatograms(ms.chrom.1,ms.chrom.2)

mix <- c(ms.chrom.1,ms.chrom.2)
plotChromatogramsOverlay(normalize(mix),lwd = 2, peakType = "none",col = c("red","blue"))
plotChromatogramsOverlay(normalize(mix),lwd = 2, peakType = "none",col = c("red","blue"))

plot(ms.chrom.1)
plot(ms.chrom.2)

ms.chrom.mix <- featureChromatograms(xcms.xcms  ,features = feature_group$feature.id)
cor.matrix <- compareChromatograms(ms.chrom.mix)


cor.matrix[is.na(cor.matrix)] <- 0




compare_eic <- function( MS.network ,xcms.xcms){

  #x <- MS.network[[2]]
  compound <- x[["compound"]]
  adduct <- x[["adduct"]]%>%
    dplyr::filter(rt.filter)

  adduct
  chrom <- featureChromatograms(xcms.xcms , features = adduct$feature.id)
  cor.matrix <- compareChromatograms(chrom)
  cor.to.main.peak <- cor.matrix[which.max(adduct$feature.intb),]
  cor.to.main.peak[is.na(cor.to.main.peak) ] <- 0
  adduct$cor.to.main.peak <- cor.to.main.peak

  x[["chromatogram"]] <- chrom
  x[["adduct"]] <- adduct
  return(x)

}

MS.network <- lapply(MS.network, compare_eic,xcms.xcms)
export::graph2ppt(file = "temp.xlsx",width = 7,height = 4,append = T)




tra.adduct <- MS.network.pos[["2"]][["adduct"]]%>%filter(cor.to.main.peak>0.5)

isopat <- isotopes_pattern_enviPat("C19H22Cl1N5O1")

tra.adduct <- tra.adduct[c(1,2,4,6),]
isopat <- isopat[c(1,3,4,6),]

isoadduct <- data.frame(adduct = tra.adduct$adduct,
                        abundance = isopat$abundance,
                        intb = tra.adduct$feature.intb)
isoadduct$intb[1] <- 35257012.120*0.6
isoadduct <- isoadduct%>%
  mutate(intb = intb /21154207.272*100)%>%
  pivot_longer( 2:3,names_to = "f",values_to = "v" )

library(randomcoloR)
ggplot(isoadduct)+
  geom_bar(aes(x = adduct , y = v , fill = f),stat = "identity",position = "dodge2")+
  scale_fill_manual(values =randomColor(2),label = c("Relative intensity","Theorical isotope abundance") )+
  labs(fill = "Ion type", x = "Adduct/isotope", y = "Relative intensity/\nTheoretical abundance")+
  theme_bw()

plot_adduct_distribution(MS.network.neg,3)+
  xlim(c(0,100))



Sys.time()
####  2022-08-10 21:04:16 CST

library(tibble)
setClass("MS_exp" ,slots = list(

  name = "character",
  MSC_id = "character",
  Pre_processing = "tbl",
  Column = "tbl",
  PhaseA  ="tbl",
  PhaseB = "tbl",
  Gradient = "tbl",
  MS_info = "tbl"

))

MS_exp <- function(){
  new("MS_exp")
}




use_r("MS_exp-class")





project.dir ="d:/2022_08_08-WYQ/"
raw.data.dir ="d:/2022_08_08-WYQ/Data/"



metabolomic_workflow(project.dir ,raw.data.dir )



prototype = list(
  General = tibble(
    "MSC_id" =  "MSC0001",
    "name" =  "Metabolomics",
    "creat_time" = paste0(Sys.time()),
    "instrument" = "SCIEX TripleTOF 6600",
    "data_aquisition" = "DDA TOP10",
    "link" = "",
    "note" = ""
  ),
  Pre_process = tibble(

  )
)
# Fri Aug 12 15:10:20 2022 ------------------------------


object <- MS_Exp(General )
object
a <- c(object,object)
a
length(a)
.remove_MS_Exp(a , c(1))
a[-2]
a[1] <- object



# Fri Aug 12 15:15:30 2022 ------------------------------
MS_workbook <- openxlsx::createWorkbook()
sapply(slotNames(x), function(y){
  openxlsx::addWorksheet( wb = MS_workbook,sheetName = y)
} )
x@General %>%t %>%as.data.frame()%>%
  select(value = 1)%>%
  openxlsx::writeData(wb = MS_workbook,sheet = 1,rowNames = T)
x@Pre_process %>%as.data.frame()%>%
  openxlsx::writeData(wb = MS_workbook,sheet = 2,rowNames = F)
x@Moblie_phase_A %>%as.data.frame()%>%
  openxlsx::writeData(wb = MS_workbook,sheet = 3,rowNames = F)
x@Moblie_phase_B %>%as.data.frame()%>%
  openxlsx::writeData(wb = MS_workbook,sheet = 4,rowNames = F)
x@Chroma_column%>%unlist%>%as.data.frame()%>%as.data.frame()%>%
  select(value = 1)%>%
  openxlsx::writeData(wb = MS_workbook,sheet = 5,rowNames = T)
x@Chroma_gradient %>%as.data.frame()%>%
  openxlsx::writeData(wb = MS_workbook,sheet = 6,rowNames = F)
x@Mass_Spectrum%>%unlist%>%as.data.frame()%>%
  select(value = 1)%>%
  openxlsx::writeData(wb = MS_workbook,sheet = 7,rowNames = T)

saveWorkbook(wb = MS_workbook, file = "b.xlsx",overwrite = T)



# Fri Aug 12 16:21:35 2022 ------------------------------
MS_Experiment@General$MSE_id %>%
  str_extract(pattern = "[0-9]+")%>%
  as.numeric()%>%
  `+`(1)%>%
  sprintf("%05d",.)%>%
  paste0("MSE",.)


i= 1
min((i+1),length(MS_Experiment)):length(MS_Experiment)



MS_Experiment@General$MSE_id %>%
  str_extract(pattern = "[0-9]+")%>%
  as.numeric()%>%
  setdiff(1:1e5,.)%>%
  min()
  `+`(1)%>%
  sprintf("%05d",.)%>%
  paste0("MSE",.)




  data("starwars")
  colnames(starwars )
  vars <- c("mass", "height")
  mutate(starwars, prod = .data[[vars[[1]]]] * .data[[vars[[2]]]])
  vars <- c("aaa", "bbb")
  mutate(starwars, "{vars[1]}":= 0 ) %>%colnames()



  MS_Experiment@General$MSE_id %>%
    str_extract(pattern = "[0-9]+")%>%
    as.numeric()%>%
    setdiff(1:10,.)%>%
    min

library(stringi)
a <- stri_rand_strings(10,5,pattern = "[A-D]")
b <- stri_rand_strings(4,5,pattern = "[A-D]")

stringdist::stringdisttmatrix(a,b)
a
b

gunzip()


use_r("dev_openxlsx")




df<- starwars

wb <- df_to_wb(df)



lipid_blast <- as.data.frame(spectraData(spectra.database))
hmdb <- as.data.frame(spectraData(hmdb.exp.spectra))



library(devtools)
load_all()
use_r("ms_ana-function")




load("d:/2022_08_08-WYQ/ms.ana.2022-08-12.Rdata")
ms.ana <-edit_sample_info(ms.ana)



for (i in 1:5) {


  plot_adduct_distribution(MS.network.pos,i)->p

  export::graph2ppt(p , file = "MS_Network.pptx",width = 7,height = 4,append = T)
  plot_adduct_chromatogram(MS.network.pos,i,
                           rt.filter = T,cor.thresh = 0.5,
                           norm = F,
                           move = T) ->p
  export::graph2ppt(p , file = "MS_Network.pptx",width = 7,height = 4,append = T)


  plot_adduct_distribution(MS.network.neg,i)->p

  export::graph2ppt(p , file = "MS_Network.pptx",width = 7,height = 4,append = T)
  plot_adduct_chromatogram(MS.network.neg,i,
                           rt.filter = T,cor.thresh = 0.5,
                           norm = F,
                           move = T) ->p
  export::graph2ppt(p , file = "MS_Network.pptx",width = 7,height = 4,append = T)
}




xcms.s1.pos <- readMSData(files = "d:/2022.6.27.STD/mzML/pos/Sample2.mzML",
                          mode = "onDisk")

tic.s1.pos <- xcms::chromatogram(xcms.s1.pos,aggregationFun = "sum")
plot(tic.s1.pos ,main = "")
title( main = "Serum metabolomic TIC\nPositive",
       adj = 0)
export::graph2ppt(file = "tic.pptx",width = 7,height = 4,append = T)
### neg
xcms.s1.neg <- readMSData(files = "d:/2022.6.27.STD/mzML/neg/Sample2.mzML",
                          mode = "onDisk")

tic.s1.neg <- xcms::chromatogram(xcms.s1.neg,aggregationFun = "sum")
plot(tic.s1.neg ,main = "")
title( main = "Serum metabolomic TIC\nnegitive",
       adj = 0)
export::graph2ppt(file = "tic.pptx",width = 7,height = 4,append = T)





### IS
xcms.s1.pos <- readMSData(files = "d:/2022.6.27.STD/mzML/pos/Sample7.mzML",
                          mode = "onDisk")

tic.s1.pos <- xcms::chromatogram(xcms.s1.pos,aggregationFun = "sum")
plot(tic.s1.pos ,main = "")
title( main = "Standard metabolomic TIC\nPositive",
       adj = 0)
export::graph2ppt(file = "tic.pptx",width = 7,height = 4,append = T)
### neg
xcms.s1.neg <- readMSData(files = "d:/2022.6.27.STD/mzML/neg/Sample7.mzML",
                          mode = "onDisk")

tic.s1.neg <- xcms::chromatogram(xcms.s1.neg,aggregationFun = "sum")
plot(tic.s1.neg ,main = "")
title( main = "Standard metabolomic TIC\nnegitive",
       adj = 0)
export::graph2ppt(file = "tic.pptx",width = 7,height = 4,append = T)











# Mon Aug 22 15:21:21 2022 ------------------------------
load("d:/2022_08_08-WYQ/ms.ana.2022-08-12.Rdata")


cbind(SummarizedExperiment::rowData(xcms.sum),
      SummarizedExperiment::assay(xcms.sum))->a

# Mon Aug 22 17:49:20 2022 ------------------------------
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


nutrimouse.sgccda <- wrapper.sgccda(X=data,
                                    Y = Y,
                                    design = design,
                                    keepX = list(gene=c(10,10), lipid=c(15,15)),
                                    ncomp = 2,
                                    scheme = "horst")

circosPlot(nutrimouse.sgccda, cutoff = 0.7)
## links widths based on strength of their similarity
circosPlot(nutrimouse.sgccda, cutoff = 0.7, linkWidth = c(1, 10))
## custom legend
circosPlot(nutrimouse.sgccda, cutoff = 0.7, size.legend = 1.1)

## more customisation
circosPlot(nutrimouse.sgccda, cutoff = 0.7, size.legend = 1.1, color.Y = 1:5,
           color.blocks = c("green","brown"), color.cor = c("magenta", "purple"))

par(mfrow=c(2,2))

circosPlot(nutrimouse.sgccda, cutoff = 0.7, size.legend = 1.1)
## also show intra-block correlations
circosPlot(nutrimouse.sgccda, cutoff = 0.7,
           size.legend = 1.1, showIntraLinks = TRUE)
## show lines
circosPlot(nutrimouse.sgccda, cutoff = 0.7, line = TRUE, ncol.legend = 1,
           size.legend = 1.1, showIntraLinks = TRUE)
## custom line legends
circosPlot(nutrimouse.sgccda, cutoff = 0.7, line = TRUE, ncol.legend = 2,
           size.legend = 1.1, showIntraLinks = TRUE)
par(mfrow=c(1,1))

## adjust feature and block names radially
circosPlot(nutrimouse.sgccda, cutoff = 0.7, size.legend = 1.1)
circosPlot(nutrimouse.sgccda, cutoff = 0.7, size.legend = 1.1,
           var.adj = 0.8, block.labels.adj = -0.5)
## ---  example using breast.TCGA data
data("breast.TCGA")
data = list(mrna = breast.TCGA$data.train$mrna,
            mirna = breast.TCGA$data.train$mirna,
            protein = breast.TCGA$data.train$protein)
list.keepX = list(mrna = rep(20, 2), mirna = rep(10,2), protein = c(10, 2))

TCGA.block.splsda = block.splsda(X = data,
                                 Y =breast.TCGA$data.train$subtype,
                                 ncomp = 2, keepX = list.keepX,
                                 design = 'full')
circosPlot(TCGA.block.splsda, cutoff = 0.7, line=TRUE)
## show only first 2 blocks
circosPlot(TCGA.block.splsda, cutoff = 0.7, line=TRUE, blocks = c(1,2))
## show only correlations including the mrna block features
circosPlot(TCGA.block.splsda, cutoff = 0.7, blocks.link = 'mrna')

data("breast.TCGA")
data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna)
list.keepX = list(mrna = rep(20, 2), mirna = rep(10,2))
list.keepY = c(rep(10, 2))

TCGA.block.spls = block.spls(X = data,
                             Y = breast.TCGA$data.train$protein,
                             ncomp = 2, keepX = list.keepX,
                             keepY = list.keepY, design = 'full')
circosPlot(TCGA.block.spls, group = breast.TCGA$data.train$subtype, cutoff = 0.7,
           Y.name = 'protein')
## only show links including mrna
circosPlot(TCGA.block.spls, group = breast.TCGA$data.train$subtype, cutoff = 0.7,
           Y.name = 'protein', blocks.link = 'mrna')



# Mon Aug 22 17:39:35 2022 ------------------------------
library(devtools)
load_all()


use_r("Statistic-function")
load("d:/2022_08_08-WYQ/ms.ana.2022-08-12.Rdata")

plot(pca.pca,3)
plot(pca.data$p1,pca.data$p2)

library(ggsci)
show_col(pal_futurama()(10))



# Tue Aug 23 20:13:20 2022 ------------------------------
load("d:/2022_08_19-Lirui/ms.ana.2022-08-22.Rdata")
ms.ana<-edit_sample_info(ms.ana)

a <- filter(feature , mz >616,mz < 617)
a <- filter(feature , rt >30,rt < 40,mz >600,mz <620)
plot_feature_intensity_distribution(ms.ana,"FT11071_pos")+
  theme_classic()
export::graph2ppt(file = "RSD.pptx",width = 5,height = 3,append = T)
plot_feature_intensity_distribution(ms.ana,"FT06896_neg")+
  theme_classic()

export::graph2ppt(file = "RSD.pptx",width = 5,height = 3,append = T)
## 5fu
plot_feature_intensity_distribution(ms.ana,"FT00256_neg")+
  theme_classic()



xcms.xcms <- ms.ana$xcms.negative
library(xcms)
fileNames(xcms.xcms)
ms.chrom <-featureChromatograms(xcms.xcms , features = 6896)


chrom.sub <- ms.chrom[1,11]
plot(chrom.sub)

id <- "FT11071_pos"
xcms.pos <- ms.ana$xcms.positive
mz.range <- feature[feature$feature.id == id,c("mzmin","mzmax")]%>%as.numeric()
ms.chrom <-chromatogram(filterFile( xcms.pos , 5) ,mz = mz.range)
plot(ms.chrom,xlim = c(30,50))
plot_XChromatograms(ms.chrom)+``
  xlim(c(35,43))


feature <-feature%>%
  dplyr::mutate(ion_mode = case_when(grepl(pattern = "pos",
                                       x = feature.id)~"pos",
                                    T~"neg"))
table(feature$ion_mode,feature$qc_rsd <0.3)
ggplot(feature )+
  geom_boxplot(aes(x = ion_mode , y = qc_rsd,col = ion_mode))+
  geom_hline(yintercept = 0.3,lty = "dashed",size= 2,col = "grey")+
  ylim(c(0,1))+
  theme_bw()



# Wed Aug 24 17:39:18 2022 ------------------------------




# Fri Aug 26 16:19:58 2022 ------------------------------
library(devtools)
load_all()


qe.xcms<-readMSData("d:/2022_8_23_QE_test/QC_MS2_35000.mzML",mode = "onDisk")

centwave.param <- CentWaveParam(peakwidth = c(5,30),
                                prefilter = c(3,100),
                                snthresh = 10,
                                ppm = 20)
xcms.xcms<-findChromPeaks(qe.xcms,
                          param = centwave.param)


xcms.peaks <- chromPeaks(xcms.xcms)


mz.range <- c(651.42480,
              651.42719
)
chrom <- chromatogram(xcms.xcms , mz  = mz.range)
plot(chrom)

xcms.peaks <- as.data.frame(xcms.peaks)%>%
  mutate(ppm = (mzmax-mzmin)/mz*1e6,
         mz_diff = mzmax-mzmin)

lattice::densityplot(xcms.peaks$ppm)
library(ggplot2)


ggplot(xcms.peaks,aes(x = ppm , y = mz)) +
  stat_density_2d(aes(fill= ..level..),
                  contour = T,
                  geom = "polygon",bins = 100)+
  geom_point(size = 0.1,alpha = 0.1)+
  scale_fill_gradient(low="#00000001",high = "red")+
  theme_bw()


plot_xcms_peaks_distribution(xcms.xcms,type = "l" )

ggplot(luv_colours, aes(u, v)) +
  geom_point(aes(colour = col), size = 3) +
  #scale_color_identity() +
  coord_equal()+
  guides(col = "none")


# Sun Aug 28 13:41:52 2022 ------------------------------
point.size <- data.frame( x = 1:5,
                          y = 3,
                          size = 1:5,
                          col = LETTERS[1:5])
ggplot(point.size)+
  geom_dotplot(aes(x = x,y=y,binwidth = size,col = col),alpha = 0.2)+
  scale_size_ordinal()






ggplot(xcms.peaks)+
  geom_segment(aes(x = rtmin , xend = rtmax , y = mz,yend = mz ))


# Sun Aug 28 15:18:24 2022 ------------------------------

tof.xcms <- readMSData("d:/2022.6.27.STD/mzML/pos/Sample2.mzML",
                       mode = "onDisk")
qe.xcms <- readMSData("d:/2022_8_23_QE_test/QC_MS2_35000.mzML",
                      mode = "onDisk")
centwave.param <- CentWaveParam(ppm = 20,
                                peakwidth = c(5,30),
                                snthresh = 10,
                                prefilter = c(3,100))
tof.xcms <- findChromPeaks(tof.xcms , param = centwave.param)
qe.xcms <- findChromPeaks(qe.xcms , param = centwave.param)



plot_xcms_peaks_distribution(tof.xcms,type = "o",plot.title = "TOF6600 Peaks Distribution")
export::graph2ppt(file = "QE_TOF_Comparasion", append = T,
                  width = 8,height = 6)
plot_xcms_peaks_distribution(qe.xcms,type = "o",plot.title = "QE Plus Peaks Distribution")
export::graph2ppt(file = "QE_TOF_Comparasion", append = T,
                  width = 8,height = 6)




plot_xcms_peaks_mzerror_density(tof.xcms,plot.title =  "TOF6600 Peaks mz Error Distribution")
export::graph2ppt(file = "QE_TOF_Comparasion", append = T,
                  width = 8,height = 6)

plot_xcms_peaks_mzerror_density(qe.xcms,plot.title =  "QE Plus Peaks mz Error Distribution")
export::graph2ppt(file = "QE_TOF_Comparasion", append = T,
                  width = 8,height = 6)





plot_xcms_peaks_ms1_scans(tof.xcms,plot.title = "TOF6600 Peaks MS1 Scan Count Distribution")
export::graph2ppt(file = "QE_TOF_Comparasion", append = T,
                  width = 8,height = 6)


plot_xcms_peaks_ms1_scans(qe.xcms,plot.title = "QE Plus Peaks MS1 Scan Count Distribution")
export::graph2ppt(file = "QE_TOF_Comparasion", append = T,
                  width = 8,height = 6)


plot_xcms_peaks_ms2_scans(tof.xcms,plot.title = "TOF6600 Peaks MS2 Scan Count Distribution")
export::graph2ppt(file = "QE_TOF_Comparasion", append = T,
                  width = 8,height = 6)


plot_xcms_peaks_ms2_scans(qe.xcms,plot.title = "QE Plus Peaks MS2 Scan Count Distribution")
export::graph2ppt(file = "QE_TOF_Comparasion", append = T,
                  width = 8,height = 6)



plot_xcms_peaks_SN_distribution(tof.xcms,plot.title = "TOF6600 Peaks SNR(Signal to Noise Ratio)")
export::graph2ppt(file = "QE_TOF_Comparasion", append = T,
                  width = 8,height = 6)

plot_xcms_peaks_SN_distribution(qe.xcms,plot.title = "QE Plus Peaks SNR(Signal to Noise Ratio)")
export::graph2ppt(file = "QE_TOF_Comparasion", append = T,
                  width = 8,height = 6)

### TIC
tof.tic <- chromatogram(tof.xcms)
qe.tic <- chromatogram(qe.xcms)

plot_XChromatograms(tof.tic,norm = F)+
  labs(title = "TOF6600 TIC")+
  guides(col = "none")+
  theme(text = element_text(size = 8))+
  geom_line(data = MS_Experiment@Chroma_gradient[[1]],aes(x = time*60 , y = Contentration_B*1.5e6 ))


export::graph2ppt(file = "QE_TOF_Comparasion_overview", append = T,
                  width = 6,height = 3)
plot_XChromatograms(qe.tic,norm = F)+
  labs(title = "QE Plus TIC")+
  guides(col = "none")+
  theme(text = element_text(size = 8))+
  geom_line(data = MS_Experiment@Chroma_gradient[[5]],aes(x = time*60 , y = Contentration_B*3.3e7 ))


export::graph2ppt(file = "QE_TOF_Comparasion_overview", append = T,
                  width = 6,height = 3)



### Scans
tof.scans <- fData(tof.xcms)
qe.scans <- fData(qe.xcms)
rbind(table(tof.scans$msLevel)%>%as.data.frame()%>%
  mutate(type = "tof",time = 12),
  table(qe.scans$msLevel)%>%as.data.frame()%>%
    mutate(type = "qe",time = 24)
  )%>%select(mslevel = Var1,everything())-> plot.data
ggplot(plot.data)+
  geom_bar(aes(x = type , y = Freq , fill = mslevel),alpha = 0.3,stat = "identity",position = "dodge")+
  geom_bar(aes(x = type , y = Freq/time*10 , col = mslevel),stat = "identity",position = "dodge")+
  scale_x_discrete(labels = c("QE Plus","TOF6600"))+
  scale_color_manual(values = c("#F8766D","#00BFC4"),label = c("MS1","MS2"))+
  scale_fill_manual(values = c("#F8766D","#00BFC4"),label = c("MS1","MS2"))+
  labs(title = "Spectrum count",
       fill ="Total count",col = "Normalize to time",
       x = "Instrument",y = "Spectrum count")+
  guides(fill = guide_legend(order = 1) , col)+
  theme_bw()+
  theme(text = element_text(size = 8))
export::graph2ppt(file = "QE_TOF_Comparasion_overview", append = T,
                  width = 4,height = 6)

# Mon Aug 29 15:07:23 2022 ------------------------------
metabolomic_workflow_single_file(tof.xcms)
tof.anno <- metabolomic_workflow_single_file(tof.xcms)
qe.anno <- metabolomic_workflow_single_file(qe.xcms)





chromPeaks(xcms.xcms)%>%
  as.data.frame()%>%
  mutate(mz_error = mzmax - mzmin )%>%
  pull(mz_error)%>%
  max


xcms.xcms <- tof.xcms
peak.density.param <- PeakDensityParam(sampleGroups = "a",
                                       minFraction = 0.5,bw = 30,
                                       binSize = 0.5)

xcms.xcms <- groupChromPeaks(xcms.xcms,param = peak.density.param)
nrow(featureDefinitions(xcms.xcms))

# Wed Aug 31 09:45:10 2022 ------------------------------
plot.data <- rbind(data.frame( type= "tof",
                  peaks = nrow(chromPeaks(tof.xcms)),
                  features = nrow(tof.anno),
                  annoated = sum(tof.anno$score >0))%>%
        pivot_longer(2:4,names_to = "n",
                     values_to = "v"),
      data.frame( type= "qe",
                  peaks = nrow(chromPeaks(qe.xcms)),
                  features = nrow(qe.anno),
                  annoated = sum(qe.anno$score >0))%>%
        pivot_longer(2:4,names_to = "n",
                     values_to = "v"))
library(ggsci)
ggplot(plot.data)+
  geom_bar(aes(x = type , y = v , fill = n),alpha = 0.8,col="black",stat = "identity",position = "dodge")+
  scale_x_discrete(labels = c("QE Plus","TOF6600"))+
  scale_fill_aaas()+
  #scale_color_manual(values = c("#F8766D","#00BFC4"),label = c("MS1","MS2"))+
  #scale_fill_manual(values = c("#F8766D","#00BFC4"),label = c("MS1","MS2"))+
  labs(title = "Peaks count",
       fill ="Total count",col = "Normalize to time",
       x = "Instrument",y = "Spectrum count")+
  guides(fill = guide_legend(order = 1) , col)+
  theme_bw()+
  theme(text = element_text(size = 8))
export::graph2ppt(file = "QE_TOF_Comparasion_overview", append = T,
                  width = 4,height = 6)









# Thu Sep  1 11:44:44 2022 ------------------------------
metabolomic_workflow(project.dir = "d:/2022_08_30-Lirui/",
                     raw.data.dir = "d:/2022_08_30-Lirui/Data/")



load("d:/2022_08_30-Lirui/ms.ana.2022-09-01.Rdata")
feature.neg <- ms.ana$annotation.negative$annotation.table
a <- feature.neg%>%
  filter(rt >30&rt <40,mz >600,mz <700)




# Fri Sep  2 10:01:18 2022 ------------------------------
qe.fullsacn.70k <- readMSData("d:/2022.9.1/Fullscan_R70000_pos_neg_Sample_2.mzML",
                              mode = "onDisk")

fData(qe.fullsacn.70k) -> qe.scans
qe.scans.pos <- qe.scans %>%
  filter(polarity ==1,msLevel ==1)


time.diff.70k <- diff(qe.scans.pos$retentionTime)
boxplot(time.diff.70k)


qe.fullsacn.14k <- readMSData("d:/2022.9.1/Sample_1_Fullscan_R140000_pos_neg.mzML",
                              mode = "onDisk")

fData(qe.fullsacn.14k) -> qe.scans
qe.scans.pos <- qe.scans %>%
  filter(polarity ==1,msLevel ==1)


time.diff.14k <- diff(qe.scans.pos$retentionTime)
boxplot(time.diff.14k)


qe.fullsacn.only <- readMSData("d:/2022.9.1/Std_Ace_4_1_fullms.mzML",
                              mode = "onDisk")

fData(qe.fullsacn.only) -> qe.scans
qe.scans.neg <- qe.scans %>%
  filter(polarity ==0,msLevel ==1)


time.diff.neg <- diff(qe.scans.neg$retentionTime)
boxplot(time.diff.neg)



# Fri Sep  2 12:21:14 2022 ------------------------------
load("d:/2022_08_30-Lirui/ms.ana.2022-09-01.Rdata")



ms.ana <- get_feature(ms.ana,qc_rsd_thresh = 999)
ms.features <- ms.ana$feature%>%
  filter(grepl("neg",feature.id))%>%
  filter(rt >30,
         rt <40)%>%
  select(c(1,2,3,23,24,25))




plot_feature_intensity_distribution(ms.ana = ms.ana ,
                                    feature.id.to.plot ="FT09944_neg" )


is.chrom <- featureChromatograms(ms.ana$xcms.negative,
                     features = 291)



plot(is.chrom)


ms.features <- ms.ana$feature%>%
  mutate(ion_mode = case_when(grepl("neg",feature.id)~ 0,
                   grepl("pos",feature.id)~1))

ggplot(ms.features) +
   geom_boxplot(aes(x = as.factor(ion_mode) , y = qc_rsd))+
  geom_hline(yintercept = 0.3)+
  ylim(c(0,0.5))+
  scale_x_discrete(label = c("negative","positive"))+
  labs(title = "2022.8",x = "")




#gout.features <- openxlsx::read.xlsx("../../Projecct/2022.4.5.Gout/Discovery/features.no.svr.xlsx")
ms.features <- gout.features%>%
  mutate(ion_mode = case_when(grepl("neg",feature.id)~ 0,
                              grepl("pos",feature.id)~1))

ggplot(ms.features) +
  geom_boxplot(aes(x = as.factor(ion_mode) , y = qc.rsd))+
  geom_hline(yintercept = 0.3)+
  ylim(c(0,0.5))+
  scale_x_discrete(label = c("negative","positive"))+
  labs(title = "2021.6",x = "")


y <- rnorm(1000, 150, 10)

cutoff <- quantile(y, probs = 0.95)

hist.y <- density(y, from = 100, to = 200) %$%
  data.frame(x = x, y = y) %>%
  mutate(area = x >= cutoff)

the.plot <- ggplot(data = hist.y, aes(x = x, ymin = 0, ymax = y, fill = area)) +
  geom_ribbon() +
  geom_line(aes(y = y)) +
  geom_vline(xintercept = cutoff, color = 'red') +
  annotate(geom = 'text', x = cutoff, y = 0.025, color = 'red', label = 'VaR 95%', hjust = -0.1)
the.plot


# Fri Sep  2 14:30:49 2022 ------------------------------
library(xcms)
qe.fullscan <-readMSData("d:/2022.9.1/Fullscan_R70000_pos_neg_Sample_1.mzML",
                         mode = "onDisk")
qe.fullscan.pos  <- filterPolarity(qe.fullscan,1)
qe.fullscan.neg  <- filterPolarity(qe.fullscan,0)

qe.fullscan.pos <- findChromPeaks(qe.fullscan.pos,param = CentWaveParam(ppm = 5))
plot_xcms_peaks_distribution(qe.fullscan.pos,type = "l")
plot_xcms_peaks_ms1_scans(qe.fullscan.pos)
plot_xcms_peaks_SN_distribution(qe.fullscan.pos)
nrow(chromPeakData(qe.fullscan.pos))

plot_xcms_peaks_Chromatogram(qe.fullscan.pos , 1949)

qe.peaks <- chromPeaks(qe.fullscan.pos) %>%
  as.data.frame()%>%
  mutate(pw  = rtmax - rtmin)




use_r("MSdev-class")# Fri Sep  2 18:57:29 2022 ------------------------------
load_all()
object <- MSdev()
#use_r("MSdev-function")

use_r()

load("../../Projecct/2022.1.8_MS.demo/Demo3/MSdev_2022_09_03.Rdata")



test.dir <- "d:/test.dir"
dir.create(test.dir)
test.files <- paste0(test.dir , "/fullscan_pos_neg_",1:9,".raw")
file.create(test.files)


test.msdev<- MSdev(rawDataDir = "d:/test.dir/raw.data/")


raw.files <- dir(path = projectInfo$rawDataDir,
                 pattern = paste0(projectInfo$rawDataFormat,"$"),
                 full.names = T)
raw.files <- test.files

.select_char <- function(char_vector){
  if (all(is.na(char_vector))) {
    return(NA)
  }
  max(char_vector,na.rm = T)
}
sample.info <- data.frame(raw.files = raw.files)%>%
  dplyr::filter(!grepl(pattern = "condition",x = raw.files))%>%
  dplyr::mutate(raw.file.positive = case_when(grepl(pattern = "pos", x = raw.files,ignore.case = T)~raw.files),
                raw.file.negative = case_when(grepl(pattern = "neg", x = raw.files,ignore.case = T)~raw.files))%>%
  dplyr::mutate(sample.abbreviation= gsub(pattern = paste0("pos|neg|",projectInfo$rawDataFormat,"$"),
                                          x = basename(raw.files) ,
                                          ignore.case = T,
                                          replacement = ""),
                sample.abbreviation = tolower(sample.abbreviation))%>%
  dplyr::group_by(sample.abbreviation)%>%
  dplyr::mutate(raw.file.positive = .select_char (raw.file.positive),
                raw.file.negative = .select_char (raw.file.negative))%>%
  dplyr::ungroup()%>%
  dplyr::distinct(sample.abbreviation,.keep_all = T)%>%
  dplyr::mutate(analysis.time.positive =as.character( file.info(raw.file.positive)$mtime),
                analysis.time.negative = as.character(file.info(raw.file.negative)$mtime))%>%
  dplyr::mutate(sample.type = case_when(grepl(pattern = "QC",x = sample.abbreviation,ignore.case = T)~ "QC",
                                        grepl(pattern = "blank|blk",x = sample.abbreviation,ignore.case = T)~ "Blank",
                                        T~"Sample"),
                .before = sample.abbreviation)%>%
  dplyr::mutate(sample.name = paste0(sample.type,str_pad(1:nrow(.), ceiling(log10(nrow(.))),pad = "0")),
                .before = sample.type)%>%
  dplyr::mutate(msData.file.positive = case_when(is.na(raw.file.positive)~raw.file.positive,
                                                 T~paste0(msData.dir,"pos/",sample.name,".mzML")),
                msData.file.negative = case_when(is.na(raw.file.negative)~raw.file.positive,
                                                 T~paste0(msData.dir,"neg/",sample.name,".mzML")))%>%
  dplyr::select(sample.name,sample.type,sample.abbreviation,
                raw.file.positive,raw.file.negative,
                analysis.time.positive,analysis.time.negative,
                msData.file.positive,msData.file.negative)









sample.info <- data.frame(raw.files = raw.files)%>%
  dplyr::filter(!grepl(pattern = "condition",x = raw.files))%>%
  dplyr::mutate(polarity = case_when(grepl(pattern = "pos", x = raw.files,ignore.case = T)~"raw.file.positive",
                                     grepl(pattern = "neg", x = raw.files,ignore.case = T)~"raw.file.negative",

                                     T~"error"),
                sample.abbreviation= gsub(pattern = paste0("pos|neg|",projectInfo$rawDataFormat,"$"),
                                          x = basename(raw.files) ,
                                          ignore.case = T,
                                          replacement = ""),
                sample.abbreviation = tolower(sample.abbreviation))%>%
  pivot_wider(names_from = polarity , values_from = raw.files)%>%
  dplyr::mutate(analysis.time.positive =as.character( file.info(raw.file.positive)$mtime),
                analysis.time.negative = as.character(file.info(raw.file.negative)$mtime))%>%
  dplyr::mutate(sample.type = case_when(grepl(pattern = "QC",x = sample.abbreviation,ignore.case = T)~ "QC",
                                        grepl(pattern = "blank|blk",x = sample.abbreviation,ignore.case = T)~ "Blank",
                                        T~"Sample"),
                .before = sample.abbreviation)%>%
  dplyr::mutate(sample.name = paste0(sample.type,str_pad(1:nrow(.), ceiling(log10(nrow(.))),pad = "0")),
                .before = sample.type)%>%
  dplyr::mutate(msData.file.positive = case_when(is.na(raw.file.positive)~raw.file.positive,
                                                 T~paste0(msData.dir,"pos/",sample.name,".mzML")),
                msData.file.negative = case_when(is.na(raw.file.negative)~raw.file.positive,
                                                 T~paste0(msData.dir,"neg/",sample.name,".mzML")))%>%
  dplyr::select(sample.name,sample.type,sample.abbreviation,
                raw.file.positive,raw.file.negative,
                analysis.time.positive,analysis.time.negative,
                msData.file.positive,msData.file.negative)


# Sun Sep  4 23:19:09 2022 ------------------------------


test.msdev<- MSdev(rawDataDir = "d:/test.dir/raw.data/")
checkSampleInfo(test.msdev)











# Mon Sep  5 10:28:08 2022 ------------------------------
library(devtools)
load_all()
use_r("dev_proteoWizard")
MS_dev_QE <-checkSampleInfo(MSdev)



xcmsProcessingMS1(msDataFiles = MS_dev_QE@sampleInfo$msData.file.positive,
                  ion_mode = 1,
                  peaksGroup =MS_dev_QE@sampleInfo$sample.type,
                  centWaveParam =xcms::CentWaveParam(ppm = 20,
                                                     peakwidth = c(5,50))
)->a.ppm.20
xcmsProcessingMS1(msDataFiles = MS_dev_QE@sampleInfo$msData.file.positive,
                  ion_mode = 1,
                  peaksGroup =MS_dev_QE@sampleInfo$sample.type,
                  centWaveParam =xcms::CentWaveParam(ppm = 5,
                                                     peakwidth = c(5,50))
)->a.ppm.5

b <- featureDefinitions(a)%>%as.data.frame()

nrow(chromPeaks(a.ppm.5))
nrow(featureDefinitions(a.ppm.5))




library(enviPat)
a <- isotopes
b <- adducts

openxlsx::write.xlsx(list(a,b),file = "a.xlsx")



hist(sapply(object@spectra$positiveFeatureMS2Map, length))
a[a > 5] <- 5

hist(a,xlim = c(0,5),breaks = 0:5,ylim = c(0,300),
     main = "Total 3654 features",
     xlab = "MS2 count")


ggplot(a)+
  geom_density_2d_filled(aes( x= precursorMz, y = precursorMz))



a <- spectra.pos%>%
  addProcessing(normalizeSpectra)%>%
  applyProcessing()
a$maxIntensity







library(devtools)
load_all()
load("d:/2022.9.1_QE.Test/MSdev_2022_09_05.Rdata")
MS_dev_QE <- MSdev
rm(MSdev)
sum(sapply(MS_dev_QE@annotation$negativeCandidate, length)>0)
sum(sapply(MS_dev_QE@annotation$positiveCandidate, length)>0)



load("d:/2022_08_30-Lirui/ms.ana.2022-09-01.Rdata")


MS_dev_QE@annotation[["positiveCandidate"]][[31]] ->refSpec
MS_dev_QE@spectra[["positiveFeatureMS2"]][[31]] ->expSpec


files <- dir("f:/YHY.Cloud/LR/MSdemo/",pattern = ".wiff$",full.names = T)
files.mzml <- sub(x = files,pattern = "wiff",replacement = "mzML")
msconvert_wiff2mzML(files,files.mzml)






expSpec <- object@spectra$positiveFeatureMS2[[58]]
refSpec <- object@annotation$positiveCandidate[[58]]


featurematrix <- assay(featurePos)%>%t

toplot <- featurematrix%>%
  scale()

ComplexHeatmap::Heatmap(t(toplot),
                        cluster_columns = F,
                        show_row_names = F,
                        show_row_dend = F)


MS_dev_tof <- MSdev(rawDataDir = "d:/2022.9.1_TOF.test/raw.data",
                   projectDir = "d:/2022.9.1_TOF.test",
                   experimentInfo = MS_Experiment[1])
MS_dev_tof <- checkSampleInfo(MS_dev_tof)
MS_dev_tof <- msConvert(MS_dev_tof)
MS_dev_tof <- xcmsProcessing_fullscan_DDA(MS_dev_tof)
MS_dev_tof <- extractSpectra_fullscan_DDA(MS_dev_tof)
MS_dev_tof <- featureSpectra_fullscan_DDA(MS_dev_tof)
MS_dev_tof <- featureCandidate(MS_dev_tof,mz.ppm = 20)
MS_dev_tof <- annotateMSdev(MS_dev_tof)
MS_dev_tof <- dropSpectra(MS_dev_tof)
saveMSdev(MS_dev_tof)




nrow(MS_dev_QE@xcmsData$positiveMS1%>%featureDefinitions())
nrow(MS_dev_QE@xcmsData$negativeMS1%>%featureDefinitions())
nrow(MS_dev_tof@xcmsData$negativeMS1%>%featureDefinitions())
nrow(MS_dev_tof@xcmsData$positiveMS1%>%featureDefinitions())


sum(sapply(MS_dev_QE@annotation$positiveAnnotation, `[`,"score")>0)
sum(sapply(MS_dev_QE@annotation$negativeAnnotation, `[`,"score")>0)

sum(sapply(MS_dev_tof@annotation$positiveAnnotation, `[`,"score")>0)
sum(sapply(MS_dev_tof@annotation$negativeAnnotation, `[`,"score")>0)


MS_dev_QE@spectra[["positiveMS2"]]->a
plot(a$rtime,a$precursorMz,main = "QE MS2 precursor mz ")



MS_dev_tof@spectra[["positiveMS2"]]->a
plot(a$rtime,a$precursorMz,main = "TOF MS2 precursor mz ")













qe.dda <- readMSData("d:/2022_8_23_QE_test/QC_MS2_35000.mzML",mode = "onDisk")
qe.scans <- fData(qe.dda)%>%
  dplyr::filter(msLevel == 2)

plot(qe.scans$retentionTime,qe.scans$precursorMZ,ylim = c(416.371,416.373))
x.chrom <- chromatogram(qe.dda,mz =  c(416.371,416.373))
plot(x.chrom)



sum(sapply(MS_dev_QE@annotation[["positiveCandidate"]], length)>0)
sum(sapply(MS_dev_QE@annotation[["negativeCandidate"]], length)>0)

sum(sapply(MS_dev_tof@annotation[["positiveCandidate"]], length)>0)
sum(sapply(MS_dev_tof@annotation[["negativeCandidate"]], length)>0)




library(xcms)
qe.dda <- readMSData("d:/2022.9.8.QE_from.FSH/T2_UP.mzML",mode = "onDisk")
qe.scans <- fData(qe.dda)%>%
  dplyr::filter(msLevel == 2)
plot(qe.scans$retentionTime,qe.scans$precursorMZ)
x.chrom <- chromatogram(qe.dda,mz =  c(416.371,416.373))
plot(x.chrom)


mzlist <- a
clusterMz <-function(mzlist){


}



plot_xcms_peaks_distribution(MS_dev_tof@xcmsData$positiveMS1)





library(devtools)

load_all()

plot_feature_intensity_distribution()


load("d:/2022_07_05-Lirui/ms.ana.2022-07-13.Rdata")
library(xcms)
plotAdjustedRtime(ms.ana$xcms.positive)
export::graph2png(file = "a.png",width = 10,height = 5)

load("d:/2021.10.09.LNN/mzML/ms.peaks.pos.Rdata")
plotAdjustedRtime(ms.peaks,col = randomcoloR::randomColor(length(sampleNames(ms.peaks))))
export::graph2png(file = "a.png",width = 10,height = 5)








mzfile <- dir("d:/LR/","mzML",full.names = T)
xcms.xcms <- readMSData(mzfile,mode = "onDisk")

header(xcms.xcms)

a <- openMSfile(mzfile)
runInfo(a)
close(a)




MS_dev_QE <- MSdev(rawDataDir = "d:/2022.9.16.LR.Lipidomic.Co.FuDan/rawData",
                   projectDir = "d:/2022.9.16.LR.Lipidomic.Co.FuDan",
                   experimentInfo = MS_Experiment[5])
MS_dev_QE <- checkSampleInfo(MS_dev_QE)
MS_dev_QE <- msConvert(MS_dev_QE)
MS_dev_QE <- checkSampleInfo(MS_dev_QE)
MS_dev_QE
MS_dev_QE <- xcmsProcessing_fullscan_DDA(MS_dev_QE)
MS_dev_QE <- extractSpectra_fullscan_DDA(MS_dev_QE)
MS_dev_QE <- featureSpectra_fullscan_DDA(MS_dev_QE)
MS_dev_QE <- featureCandidate(MS_dev_QE,mz.ppm = 10,spectraDatabase = "d:/MSdb/MSdb_LipidBlast_from_MSDIAL.Rdata")
MS_dev_QE <- annotateMSdev(MS_dev_QE)
MS_dev_QE <- getStaData(MS_dev_QE)
saveMSdev(MS_dev_QE)
MS_dev_QE <- dropSpectra(MS_dev_QE)



object@projectInfo$MSDB_path <-"d:/MSdb/MSdb_LipidBlast_from_MSDIAL.Rdata"

a <- readMSData("d:/2022.9.16.LR.Lipidomic.Co.FuDan/msData/pos/Sample02.mzML",
                mode = "onDisk")


experimentData(a)->b


expinfo(a)
msInfo(a)->b


sum(sapply(MS_dev_QE@annotation$positiveCandidate, length)>0)
sum(sapply(MS_dev_QE@annotation$negativeCandidate, length)>0)






testfun <- function(...){

  args <- ...
  args

}

testfun(featurePos,a)->b


commandArgs()


use_r("MSdev-Sta_function")


use_r("MSdev_tools")



MS_dev_QE <- checkSampleInfo(MS_dev_QE)

