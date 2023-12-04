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
ggdraw()+draw_plot(p.sp,0,0,1,1)+
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




