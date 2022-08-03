
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

