# Wed Oct 18 21:18:36 2023 ------------------------------
### create demo

msdev.demo <- MSdev()
msdev.demo <- checkSampleInfo(msdev.demo)
msdev.demo <- msConvert_MSdev(msdev.demo)
msdev.demo <- xcmsProcessingMSdev(msdev.demo)

xcms.xcms <- msdev.demo@xcmsData$PositiveMS1
save(xcms.xcms,file = paste0(msdev.demo@projectInfo$projectDir,"/XCMSnExp_2023_10_19.rda"))

# Thu Oct 19 10:36:46 2023 ------------------------------
xcms.xcms <- filterFile(xcms.xcms , 2)

xcms.peaks.centwave.iso <- findChromPeaks(xcms.xcms,
                             param = CentWavePredIsoParam(),
                             BPPARAM = SerialParam())%>%
  groupChromPeaks(param = PeakDensityParam(1))

xcms.peaks.centwave <- findChromPeaks(xcms.xcms,
                                      param = CentWaveParam(),
                                      BPPARAM = SerialParam())%>%
  groupChromPeaks(param = PeakDensityParam(1))


a <- chromPeaks(xcms.peaks.centwave)%>%
  as.data.frame()%>%
  dplyr::mutate(cluster_ion(mz,rt))

b <- chromPeaks(xcms.peaks.centwave.iso)%>%
  as.data.frame()%>%
  dplyr::mutate(cluster_ion(mz,rt))


c <- featureDefinitions_PeakSta(xcms.peaks.centwave)



xcms.xcms %>%
  filterMz( mz.range.ppm(828.55040,25)) %>%
  filterRt(c(150,200))%>%
  plot(type = "XIC")




msdev.demo <- MSdev()
msdev.demo <- checkSampleInfo(msdev.demo)
msdev.demo <- msConvert_MSdev(msdev.demo)
msdev.demo <- xcmsProcessingMSdev(msdev.demo)
msdev.demo <- extract_Spectra_MSdev(msdev.demo)
msdev.demo <- match_Spectra_to_feature_MSdev(msdev.demo)


xcms.xcms <- load_demo("XCMS")
get_xcms_MS_report()



a <- readxl::read_excel("aaa.xlsx")

plot.data <- a%>%
  dplyr::mutate(date = as.Date(date))%>%
  dplyr::group_by(name)%>%
  dplyr::mutate(group = case_when(date <max(date)~"HUA",
                                  T~"Gout"))%>%
  dplyr::select(!id)%>%
  pivot_wider(names_from = "group",
              values_from = "date")

ggplot(plot.data)+
  geom_segment(aes(x = HUA,xend = Gout,
                   y= name,yend= name))+
  geom_point(aes(x = HUA,
                   y= name,
                 col = "#91D1C2"),show.legend = F)+
  geom_point(aes(x = Gout,
                 y= name,
                 col = "#E64B35"),show.legend = F)+
  labs( x = "Date", y = "", title = "Sample Collection before and after Gout")+
  theme_bw()+
  theme(axis.text.y = element_blank())->p

open_ggplot_win(p,6,8)


QE.POS <- load_as_var("d:/2023.09.02.DDA.mine/Pos/MSdev_2023_09_02.back.Rdata")
QE.POS.list <-QE.POS@statData[["DDA.mine.list.Positive"]][["DDA.mine.list001"]]




# Mon Oct 23 14:23:59 2023 ------------------------------
data.dir <- "d:/2023.10.22.YLF/"
for (i in c("Pos","Neg")) {
  files <- dir(paste0(data.dir,i),full.names = T)
  file.new.name <- basename(files)%>%
    paste0(data.dir,i,"/",i,"_",.)
  file.rename(files,file.new.name)

}


MS_instrument = "AB6600"
LC_condition = "Metabolomics"
library(MSdev)
MS_dev_obj <- MSdev:::MSdev("d:/2023_10_20-MHR/Data/")
MS_dev_obj <-  MSdev:::xcmsProcessingMSdev(MS_dev_obj,
                                            xcms.findpeak.param = switch(MS_instrument,
                                                                         "AB6600" = MSdev_param_set$xcms.param$findpeakparam$AB6600,
                                                                         "QEplus" = MSdev_param_set$xcms.param$findpeakparam$QEplus)
)


MS_dev_obj <-  MSdev:::extractSpectra_fullscan_DDA(MS_dev_obj)
MS_dev_obj <-  MSdev:::featureSpectra_fullscan_DDA(MS_dev_obj)
MS_dev_obj <-  MSdev:::featureCandidate(MS_dev_obj,
                               mz.ppm = 20,
                               spectraDatabase = switch(LC_condition,
                                                        "Metabolomics" = "d:/MSdb/msdb.HMDB.Rdata",
                                                        "Lipidomics" = "d:/MSdb/MSdb_LipidBlast_from_MSDIAL.Rdata",
                               ))
MS_dev_obj <-  MSdev:::annotateMSdev(MS_dev_obj)
MS_dev_obj <-  MSdev:::getStaDataMSdev(MS_dev_obj,missing = "rowmin_half",
                              MSDB.keys = switch(LC_condition,
                                                 "Metabolomics" = c("Compound_name","adduct","formula","inchikey" ,"database_origin"),
                                                 "Lipidomics" = c("Compound_name","adduct","formula","inchikey","Lipid_subclass" ,"database_origin")
                              )
)


MSdev:::saveMSdev(MS_dev_obj)
MSdev:::exportMSdev(MS_dev_obj)

MS_dev_obj <- load_as_var("d:/2023_10_20-MHR/MSdev_2023_10_23.Rdata")
MS_dev_obj <- checkSampleInfo(MS_dev_obj)


# Tue Oct 24 09:36:00 2023 ------------------------------
#load data
data(gse16873.d)
data(demo.paths)

#KEGG view: gene data only
i <- 1
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id =
                     demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873",
                   kegg.native = TRUE)
str(pv.out)
head(pv.out$plot.data.gene)
#result PNG file in current directory

#Graphviz view: gene data only
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id =
                     demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873",
                   kegg.native = FALSE, sign.pos = demo.paths$spos[i])
#result PDF file in current directory

#KEGG view: both gene and compound data
sim.cpd.data=sim.mol.data(mol.type="cpd", nmol=3000)
i <- 3
print(demo.paths$sel.paths[i])
pv.out <- pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data,
                   pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix =
                     "gse16873.cpd", keys.align = "y", kegg.native = TRUE, key.pos = demo.paths$kpos1[i])
str(pv.out)
head(pv.out$plot.data.cpd)

#multiple states in one graph
set.seed(10)
sim.cpd.data2 = matrix(sample(sim.cpd.data, 18000,
                              replace = TRUE), ncol = 6)
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                   cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.cpd.3-2s", keys.align = "y",
                   kegg.native = TRUE, match.data = FALSE, multi.state = TRUE, same.layer = TRUE)
str(pv.out)
head(pv.out$plot.data.cpd)

#result PNG file in current directory

##more examples of pathview usages are shown in the vignette.



library(CAMERA)
file <- system.file('mzML/MM14.mzML', package = "CAMERA")
xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
an   <- xsAnnotate(xs)
an   <- groupFWHM(an)
an   <- findIsotopes(an)

library(CAMERA)
file <- system.file('mzML/MM14.mzML', package = "CAMERA")
xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
an   <- xsAnnotate(xs)
an   <- groupFWHM(an)
an <- findIsotopes(an)  # optional but recommended.
#an <- groupCorr(an) # optional but very recommended step
an <- findAdducts(an,polarity="positive")
peaklist <- getPeaklist(an) # get the annotated peak list

xcms.xcms <- load_demo("XCMS")
xcms.xcms <- xcms.xcms%>%
  filterFile(1)%>%
  filterMsLevel(1)
xcms.xset <- as(xcms.xcms, "xcmsSet")

an   <- xsAnnotate(xcms.xset)
an   <- groupFWHM(an)
an <- findIsotopes(an)  # optional but recommended.
#an <- groupCorr(an) # optional but very recommended step
an <- findAdducts(an,polarity="positive")
peaklist <- getPeaklist(an) # get the annotated peak list






xcms.scan <- get_xcms_scan_Stat(xcms.xcms)%>%
  dplyr::filter(msLevel == 2)

sp.scan <- Spectra(pData(xcms.xcms)$files)
sp.scan.data <- spectraData(sp.scan)%>%as.data.frame()

mz1 <- xcms.scan$precursorMZ
rt1 <- xcms.scan$retentionTime

mz2 <- precursorMz(Spectra_database)
rt2 <- rtime(Spectra_database)


match.df <- match_mz_rt(featuredef$mzmed,featuredef$rtmed,
                        xcms.scan$precursorMZ,
                        xcms.scan$retentionTime)

match.df <- match.df%>%
  dplyr::mutate(feature_id = featuredef$feature_id[ion1],
                sp_id = xcms.scan$sp_id[ion2])


xcms.fdf <- xcms.fdf %>%
  dplyr::rowwise()%>%
  dplyr::mutate(ms2_sp = sp.ms2.data$sp_id[grep(pattern = feature_id ,
                                                x = sp.ms2.data$ms2_matched_feature
  )])%>%
  dplyr::ungroup()

xcms.fdf$ms2_id <- sapply(xcms.fdf$feature_id,
  function(i){
    sp_id <- sp.ms2.data$sp_id[which(sp.ms2.data$ms2_matched_feature==i)]
    return(sp_id)
  }
)

b <- lapply(a, function(x){
  x[c( "mz"    ,
       "rt" ,
       "ref.mz" ,
       "ref.rt" ,
       "score" ,
       "MSDB_id" )]
})%>%data.table::rbindlist()



pos.def <- featureDefinitions(msdev.demo@xcmsData$PositiveMS1)%>%
  as.data.frame()

neg.def <- featureDefinitions(msdev.demo@xcmsData$NegativeMS1)%>%
  as.data.frame()

xcms.def <- rbind(pos.def,neg.def)













