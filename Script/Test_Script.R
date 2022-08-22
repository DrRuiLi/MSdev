
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



