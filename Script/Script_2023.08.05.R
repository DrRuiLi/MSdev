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







