# Sat Aug  5 15:42:02 2023 ------------------------------
library(xcms)
msdev.pseudo.mrm <- MSdev("d:/2023.08.04.Pseudo/Result/")
msConvert_MSdev(msdev.pseudo.mrm)

sample.info <- msdev.pseudo.mrm@sampleInfo


# pos
{
  sample.info.polarity <- sample.info%>%
    dplyr::filter(!is.na(sample.info$raw.file.positive),
                  sample.type == "QC")
  xcms.xcms <- readSRMData(sample.info.polarity$msData.file.positive)
  xcms.info <- fData(xcms.xcms)
  pdf("D:/2023.08.04.Pseudo/Pseudo.MRM.Chrom.pos.pdf",width = 5,height =3)
  for (i in 1:nrow(xcms.xcms)) {
    message(i,"...")
    plot_XChromatograms(xcms.xcms[i,],F,F)+
      labs(title = xcms.info$chromatogramId[i])+
      theme(legend.position = "none")->p
    plot(p)
  }
  dev.off()


}



# neg
{
  sample.info.polarity <- sample.info%>%
    dplyr::filter(!is.na(sample.info$raw.file.negative),
                  sample.type == "QC")
  xcms.xcms <- readSRMData(sample.info.polarity$msData.file.negative)
  xcms.info <- fData(xcms.xcms)
  pdf("D:/2023.08.04.Pseudo/Pseudo.MRM.Chrom.neg.pdf",width = 5,height =3)
  for (i in 1:nrow(xcms.xcms)) {
    message(i,"...")
    plot_XChromatograms(xcms.xcms[i,],F,F)+
      labs(title = xcms.info$chromatogramId[i])+
      theme(legend.position = "none")->p
    plot(p)
  }
  dev.off()


}



# Tue Aug 15 13:38:27 2023 ------------------------------
library(ComplexHeatmap)
library(SummarizedExperiment)
pseudo.msdev <- load_as_var("d:/2023.08.03.ms.dev/MSdev_2023_08_05.Rdata")

#pos
{
  xcms.xcms <- pseudo.msdev@xcmsData$positiveMS1
  xcms.se <- get_features_from_xcms(xcms.xcms,missing = "rowmin_half" )
  xcms.rowdata <- rowData(xcms.se)%>%as.data.frame()
  xcms.coldata <- colData(xcms.se)%>%as.data.frame()
  xcms.data <- assay(xcms.se)

  ggplot(xcms.rowdata)+
    geom_histogram(aes(x = qc_rsd),
                   fill = "#F32C04",
                   col = "black",binwidth = 0.1)

  heatmap.data <- xcms.data[xcms.rowdata$qc_rsd < 0.6 ,]%>%
    log%>%t%>%scale%>%t
  png("temp.png")
  Heatmap(heatmap.data,
          use_raster=F,
          show_row_names = F,
          show_row_dend = F)
  dev.off()


}


xcms.test <- xcms.xcms[1:10,1:5]

f <- function(i,xchroms,unit.mulit){

  i.position <- which(matrix(1:length(xcms.test),
                             nrow = dim(xcms.test)[1],
                             byrow = T)==i,arr.ind = T)
  x <- i.position[1]
  y <- i.position[2]
  xchroms[x,y]@rtime <- rtime(xchroms[x,y])*unit.mulit
  xchroms[x,y]

}

xcms.test.trans <- bplapply(1:length(xcms.test),FUN = f,
              xchroms = xcms.test,BPPARAM = SerialParam())

xcms.test.trans <- XChromatograms(xcms.test.trans,nrow = dim(xcms.test)[1],byrow = T)
dim(xcms.test.trans)

plot(xcms.test.trans[2,5])
plot(xcms.test[2,5])




xcms.test.trans <- bplapply(1:length(xchroms),FUN = f,
                            xchroms = xchroms,
                            unit.mulit = unit.mulit,
                            BPPARAM = BatchtoolsParam(progressbar = T,
                                                      registryargs = batchtoolsRegistryargs(packages = c("MSnbase"))))


xcms.test <- xcms.xcms[1:100,1:5]

a <- XChromatograms_rt_unit(xcms.test)



xcms.peaks<-findChromPeaks(clean(xcms.xcms),
               param = CentWaveParam(),
               BPPARAM = SerialParam(progressbar = T))
plot(xcms.peaks)






xcms.cleaned <- clean(xcms.chroms)
xcms.empty <- sapply(unlist(xcms.cleaned),isEmpty)%>%which%>%
  lapply( get_matrix_idx,nrow = dim(xcms.chroms)[1],
       ncol = dim(xcms.chroms)[2])%>%do.call("rbind",.)





a <- xcms.xcms[5,1]
a@rtime <- a@rtime[1:2]
a@intensity <- c(100,200)
plot(a)

findChromPeaks(a,param = CentWaveParam())
rt.num <- sapply(1:length(xcms.xcms), function(x){
  x.position <- get_matrix_idx(x,dim(xcms.xcms)[1],dim(xcms.xcms)[2])
  length(rtime(xcms.xcms[x.position[1],x.position[2]]))

})

get_matrix_idx(which(rt.num==1)[1],2132,5)



arrayInd(which(rt.num==1) , .dim = dim(xcms.xcms))




a <- xcms.peaks[11,]

findChromPeaks(a,
               param = CentWaveParam(peakwidth = c(3,20),
                                     prefilter = c(2,10),
                                     fitgauss = T,
                                     snthresh = 0),
               BPPARAM = SerialParam())%>%plot()


plot(xcms.peaks[413,3])
plot(xcms.peaks[6,5])
x = 50:150
y = gaussian_functioin(x , a = 43,
                       b = 100,c = 1.23)
points(x,y)



x = order(xcms.peaks.data[,"egauss"])
pdf("d:/temp/temp.pdf")
counter = 1
for (i in x){
  message("row:",counter,"...")
  counter = counter+1
  plot(xcms.peaks[xcms.peaks.data[i,"row"],xcms.peaks.data[i,"column"]],
      main = xcms.peaks.data[i,"egauss"])
  x = rtime(xcms.peaks[xcms.peaks.data[i,"row"],
                       xcms.peaks.data[i,"column"]])
  h = xcms.peaks.data[i,"h"]
  mu = x[xcms.peaks.data[i,"scpos"]]
  sigma = xcms.peaks.data[i,"sigma"]
  y = SSgauss(x , h = h ,mu = mu,sigma = sigma)
  points(x , y)
 # abline(v = c(mu,mu+sigma*2.35,mu-sigma*2.35))



}
dev.off()






### SCPOS parse
{

  a <- xcms.peaks.data%>%
    as.data.frame()%>%
    dplyr::mutate(scwidth = scmax-scmin,
                  right = scwidth == scale*2)%>%
    dplyr::filter(scale != -1)


  ggplot(a)+
    geom_jitter(aes(x = scale, y = scwidth),width = 0.1,height = 0.1)+
    geom_abline(slope = 2)

}

xcms.test <- readMSData("d:/2023.07.19.ms.dev/Results/mzML/FS_NEG_SV3000.mzML",
                        mode = "onDisk")

xcms.test <- findChromPeaks(xcms.test,
                            param = CentWaveParam(
                              ppm = 5,
                              prefilter = c(5,10000),
                              fitgauss = T,
                              verboseColumns = T))
a <- chromPeaks(xcms.test)




sapply(xcms.peaks.xchroms[1:10])



extract_chrom(xcms.xcms ,
              mzr = peaks.data[1:10,c("mzmin","mzmax")],
              rtr = peaks.data[1:10,c("rtmin","rtmax")],
              sample = peaks.data[1:10,c("sample")])
par.list <- list()
for (i in unique(peaks.data[,"sample"])) {

  par.list[[i]] <- list(f = extract_chrom,
       xcms.xcms =filterFile(xcms.xcms,i),
       mzr = peaks.data[peaks.data[,"sample"]==i,c("mzmin","mzmax")],
       rtr = peaks.data[peaks.data[,"sample"]==i,c("rtmin","rtmax")],
       sample = i)
}


bplapply(par.list,function(x){
  x$f(xcms.xcms = x$xcms.xcms,rtr = x$rtr,mzr = x$mzr,sample = x$sample)
},BPPARAM = SnowParam(progressbar = T))->a



f.snow <- function(x){

  bplapply(1:10,sqrt,BPPARAM = SnowParam())

}
f.serial <-  function(x){

  bplapply(1:10,sqrt,BPPARAM = SerialParam())

}


bplapply(1:10,f.serial , BPPARAM = SnowParam())->a



extract_chrom <- function(xcms.a,
                          rtr,
                          mzr,
                          sample){
  xcms.a <- MSnbase::filterFile(xcms.a,sample)
  xcms::chromatogram(xcms.a,
                     rt = rtr,
                     mz = mzr,
                     BPPARAM = BiocParallel::SerialParam()
  )
}

xcms::chromatogram(xcms.a,
                   rt = rtr,
                   mz = mzr,
                   BPPARAM = BiocParallel::SnowParam(progressbar = T)
)


x.chrom.split <- lapply(peaks.id.split[1:2],
                          FUN = function(x,xcms.x,peaks.data){

                            x <- peaks.data[x,]
                            extract_chrom(xcms.x,
                                          mzr = x[,c("mzmin","mzmax")],
                                          rtr = x[,c("rtmin","rtmax")],
                                          sample = x[,"sample"])
                          },xcms.x=xcms.xcms ,peaks.data=peaks.data)


n <- 100
xcms.features <- featureDefinitions(xcms.xcms)[1:n,]

system.time(xcms.features.peaks <- featureChromatograms(xcms.xcms,
                                      features = 1:n,
                                      BPPARAM = SerialParam()))
xcms.features.peaks.info <- chromPeaks(xcms.features.peaks)

a <- sapply(xcms.features$peakidx,function(x){
  peaks.data[x, ]
})%>%do.call("rbind",.)


system.time(xcms.features.peaks <-
              featureChromatograms(xcms.xcms,
                                   features = 1:10000,
                                   n=1 ,include = "apex_within",
                                   BPPARAM = SerialParam(progressbar = T)))





