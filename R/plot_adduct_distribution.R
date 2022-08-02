plot_adduct_distribution <- function(MS.network,index){

###test
  compound <- MS.network[[index]][["compound"]]
  adduct.table <-MS.network[[index]][["adduct"]]
###
  ggplot(adduct.table)+
    geom_point(aes(x = peak.rt,y=peak.mz,size = log10(peak.intb),col = adduct),
               alpha = 0.8)+
    ggsci::scale_color_npg()+
    ggrepel::geom_text_repel(aes(x = peak.rt,y=peak.mz ,
                                 col = adduct,
                                 label = paste0(adduct,
                                                "\nmz = ",sprintf("%.4f",peak.mz))),
                             show.legend = F,direction = "x",size = 2.3)+
    scale_size(range = c(3,5),breaks = c(2,4.5,5))+
    xlim(c(0,750))+
    labs(title = compound["name"],
         subtitle = paste0("exact mass = ",compound["exact.mass"]),
         x = "Retention time",
         y = "mz",
         col = "Adduct form",
         size = "Log10(intensity)")+
    theme(text = element_text(size = 8))->adduct.plot
  #adduct.plot
  return(adduct.plot)
  #export::graph2png(file = "b.png",width = 5,height = 4)

}
