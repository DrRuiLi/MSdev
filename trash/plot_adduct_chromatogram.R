plot_adduct_chromatogram <- function(MS.network,
                                     index,
                                     rt.filter = T,
                                     cor.thresh = 0.7,
                                     norm = F,move = T){

  compound <- MS.network[[index]][["compound"]]
  adduct.table <-MS.network[[index]][["adduct"]]
  if (nrow(adduct.table)== 0 ) {
    message("No adduct detected")
    return()
  }
 #adduct.table <- dplyr::filter(adduct.table ,
 #                              rt.filter,
 #                              cor.to.main.peak >0)

  idx.to.show <-which((adduct.table$rt.filter |!rt.filter )&adduct.table$cor.to.main.peak > cor.thresh)
  chrom <-MS.network[[index]][["chromatogram"]]
  if (isEmpty(idx.to.show)) {
    message("No adduct after filter")
    return()
  }
  adduct.table <- adduct.table[idx.to.show,]
  if (length(idx.to.show)==1) {
    chrom.col <- randomcoloR::distinctColorPalette(1)
    chrom <- chrom[c(0,idx.to.show)]
  }else{
    chrom.col <- randomcoloR::distinctColorPalette(nrow(chrom))
    chrom <- chrom[idx.to.show]
  }




  plot_XChromatograms(chrom,norm = norm,move = move)+
    scale_color_manual(values = chrom.col,label = adduct.table$adduct)+
    labs(title = paste0(compound["name"]),
      x = "Retention time",y = "Intensity",col = "Adduct and isotope")




}

