plot_adduct_distribution <- function(MS.network,index,rt.filter = T){


  compound <- MS.network[[index]][["compound"]]
  adduct.table <-MS.network[[index]][["adduct"]]
  if (nrow(adduct.table)== 0 ) {
    message("No adduct detected")
    return()
  }
  if (rt.filter) {
    adduct.table <- dplyr::filter(adduct.table , rt.filter)
  }
###
  ggplot(adduct.table)+
    geom_point(aes(x = feature.rt,y=feature.mz,col = adduct),
               size =5,
               alpha = 0.3)+
    #ggsci::scale_color_npg()+
    ggrepel::geom_text_repel(aes(x = feature.rt,y=feature.mz ,
                                 col = adduct,
                                 label = paste0(adduct,
                                                "\nmz = ",sprintf("%.4f",feature.mz),
                                                "\nrt = ",sprintf("%.2f",feature.rt),
                                                "\nint = ",sprintf("%.3g",feature.intb))),
                             show.legend = F,force = 100,
                             #label.size = 0,fill = "transparent",
                             direction = "both",size = 2.3,hjust = "left")+
    #scale_size(range = c(3,5),breaks = c(2,4.5,5))+
    xlim(c(0,750))+
    labs(title = paste0(compound["name"] ),
         subtitle = paste0(unique(adduct.table$ion_mode),"\nExact mass = ",compound["exact.mass"]),
         x = "Retention time",
         y = "mz",
         col = "Adduct form",
         size = "Log10(intensity)")+
    theme_bw()+
    theme(text = element_text(size = 8))->adduct.plot
  adduct.plot
  return(adduct.plot)
  #export::graph2png(file = "b.png",width = 5,height = 4)

}
