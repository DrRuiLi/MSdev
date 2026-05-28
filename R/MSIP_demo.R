
MSIP_demo <- function(object){

  object <- msdev.Threegroup
# Tue May 21 15:44:49 2024 ------------------------------
### MS1 purity
  {

    xcms.xcms <- object@xcmsData$NegativeMS1
    xcms.fdf <- get_xcms_feature_definitions(xcms.xcms)
    idx <- 453
    xcms.temp <- xcms::filterFile(xcms.xcms, which.max(xcms::featureValues(xcms.xcms, value = "maxo")[idx, ]))
    xcms.sp <- get_xcms_Spectra(xcms.temp)
    xcms.ms1.sp <- ProtGenerics::filterMsLevel(xcms.sp, msLevel = 1)
    xcms.sp.demo <- xcms.ms1.sp[which.min(abs(rtime(xcms.ms1.sp)-xcms.fdf$rtmed[idx]))]
    xcms::plotSpec(xcms.sp.demo)+
      labs(title = 'MS1')->p1
    idx.mz <- xcms.fdf$mzmed[idx]
    xcms.sp.demo%>%
      Spectra::filterMzRange(c(idx.mz-0.5,idx.mz+0.5))%>%
      xcms::plotSpec()+
      labs(title = "Acquisition window")+
      geom_vline(xintercept = c(idx.mz-0.2,idx.mz+0.2),
                 color = "grey",size = 2, linetype = 1)+
      theme_bw()->p2
    p <- p1+p2+plot_layout(widths = c(0.8,0.2))
    open_plot_win(p,width = 10,height = 4)

  }

### purity and intensity distribution
  {

    xcms.xcms <- object@xcmsData$NegativeMS1
    xcms.fdf <- get_xcms_feature_definitions(xcms.xcms)%>%
      dplyr::filter(!is.na(C13_seed))%>%
      dplyr::mutate(acq = log10(peakMaxo) > 4.5&ms1_purity>0.8)
    p1 <- ggplot(xcms.fdf,aes(x = log10(peakMaxo), y = ms1_purity))+
      geom_point(aes(color = acq),size = 0.5,alpha = 0.5) +
      geom_density2d(color = "grey")+
      theme_bw()+
      theme(legend.position = "none")
    p2 <-   ggplot(xcms.fdf)+
      geom_histogram(aes(x = log10(peakMaxo)))+
      theme_void()
    p3 <-   ggplot(xcms.fdf)+
      geom_histogram(aes(y = ms1_purity))+
      theme_void()
    p2+plot_spacer()+p1+p3+
      plot_layout(widths = c(0.8,0.2),heights = c(0.2,0.8))+
      plot_annotation(title = paste0(sum(xcms.fdf$acq),"/",nrow(xcms.fdf)))->p
    open_plot_win(p,5,5)

  }

### Natural iso
  {
    x.iso.cfm <- msdev.combine@advancedAna$MSIP$MSIP_result[[66]]
    ratio_matrix <- x.iso.cfm$compound_info$ratio_matrix
    natural.ratio.matrix <- get_iso_natural_ratio(
      formula = x.iso.cfm$compound_info$formula,
      iso_ele = "[13]C",
      ratio_matrix = x.iso.cfm$compound_info$ratio_matrix)
    f <- function(i){
      ndf <- data.frame(iso = rownames(ratio_matrix),
                        ratio_raw = ratio_matrix[,i],
                        ratio_nature = ratio_matrix[,i]*natural.ratio.matrix[,i])%>%
        dplyr::mutate(ratio_iso = ratio_raw - ratio_nature)%>%
        pivot_longer(c(ratio_iso,ratio_nature),
                     names_to = "component",values_to = "ratio")
      ggplot(ndf)+
        geom_bar(aes(x = iso,y = ratio,
                     fill = component),stat = "identity")+
        labs(title = colnames(ratio_matrix)[i])+
        theme_classic()
    }
    p <- f(1)/f(2)/f(3)
    open_plot_win(p,10,6)
    Heatmap(natural.ratio.matrix,
            col = colramp(),show_heatmap_legend = F,
            cluster_rows = F,cluster_columns = F)%>%
      open_plot_win(5,5)

    cfmd <- x.iso.cfm$CFM_annotation
    sp.iso <- x.iso.cfm$Spectra$M1$U
    iso_count = 1
    sp.frag.data <- CFM_annotate_isotopologues(sp.iso,
                                               cfmd  = cfmd,
                                               ppm = 20,
                                               iso_count = iso_count)
    sp.frag.data <- CFM_spectra_data_int_weight(sp.frag.data,iso_count)
    fg.map <- get_frag_group_map(sp.frag.data,cfmd,iso_count = iso_count)
    heatmap.fg.map(fg.map)%>%
      open_plot_win(height = 6)
    if.map <- get_iso_form_map(fg.map)
    Heatmap(if.map$iso.form.map,
            col = colramp(colors = c("white","white","black")),
            show_column_names = F,show_heatmap_legend = F,
            cluster_rows = F,cluster_columns = F
            )%>%
      open_plot_win(height = 8)
    sp.frag.data2 <- CFM_spectra_data_remove_natural(sp.frag.data,
                                                     natural.ratio.matrix[iso_count+1,3],
                                                     if.map)
    spd1 <- data.frame(sp.frag.data,f = "nature")%>%
      dplyr::mutate(intensity =sp.frag.data$intensity-0
                      #sp.frag.data2$intensity
                      )
    sp.data <- rbind(spd1,
                     data.frame(sp.frag.data2,f = "iso"))%>%
      dplyr::filter(sp.id=="combined_sp")
    ggplot(sp.data)+
      geom_segment(aes(x = mz,xend = mz,
                       y = 0,yend = intensity,
                       colour = f),size = 1)+
      labs(y = "intensity",col = "Source")+
      theme_classic()->p
    open_plot_win(p,width = 10)
    p <- p+xlim(c(75.5,77.5))+ylim(c(0,220000))
    open_plot_win(p,width = 3)
    vis_cfm_data_atom_map(cfmd ,"Fragment129",show.label = T)

    spd1 <- data.frame(sp.frag.data,f = "nature")%>%
      dplyr::mutate(intensity =sp.frag.data$intensity-
                    sp.frag.data2$intensity
      )
    sp.data <- rbind(spd1,
                     data.frame(sp.frag.data2,f = "iso"))%>%
      dplyr::filter(sp.id=="combined_sp")
    fg.case <- sp.data%>%
      dplyr::filter(fragment_group=="FG009")%>%
      dplyr::mutate(label = paste0(f,"_",iso),
                    r = round(intensity/sum(intensity)*100,1))%>%
      dplyr::arrange(rev(label))
    ggplot(fg.case)+
      geom_bar(aes(x = fragment_group,
                   y= intensity,fill =  label),
               color = "black",
               stat = "identity")+
      scale_fill_manual(values = c(nature_M0="#00BFC4" ,
                                   nature_M1="#00BFC4" ,
                                   iso_M0 ="#F8766D" ,
                                   iso_M1 ="#F8766D"))+
      geom_text(aes(x =fragment_group ,
                    y = intensity,
                    label = r),
                position = position_stack(vjust = 0.5))+
      theme_void()->p
  open_plot_win(p,2,3)
  fg.map <- get_frag_group_map(sp.frag.data2,cfmd,iso_count = iso_count)
  heatmap.fg.map(fg.map)%>%
    open_plot_win(height = 6)
}


}
