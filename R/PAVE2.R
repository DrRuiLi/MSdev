PAVE2 <- function(object){


  rt.tol = 10
  ppm = 5

  for (i.pol in 0:1) {

    ### data
    {
      pol <- ifelse(i.pol==0,"Negative","Positive")

      xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
      xcms.net <- get_xcms_feature_connect(xcms.xcms,rt.tol = rt.tol)
      xcms.val <- featureValues(xcms.xcms,missing = 0,value = "maxo")
      xcms.pave.sample <- pData(xcms.xcms)%>%
        dplyr::filter(sample.type %in% c("S12C14N","S12C15N","S13C14N","S13C15N"))

    }

    ### theoritical mass diff match
    {

      cn.mass.diff <- get_CN_mass_diff_table(N_max = 10)[,type :="CN_label"]
      ad.mass.diff <- get_adduct_mass_diff(unique(polarity(xcms.xcms)))[,type := "adduct"]
      is.mass.diff <- get_iso_mass_diff()[,type := "isotope"]

      mass.diff.range <- range(cn.mass.diff$mass_diff,
                               ad.mass.diff$mass_diff,
                               is.mass.diff$mass_diff)
      xcms.net <- xcms.net[between.range(mz.diff,  mass.diff.range)]

      cn.match <- match_mz_foverlaps( xcms.net$mz.diff,cn.mass.diff$mass_diff,
                                ppm.base = xcms.net$mz.mean,ppm = ppm)

      ad.match <- match_mz_foverlaps( xcms.net$mz.diff,ad.mass.diff$mass_diff,
                                ppm.base = xcms.net$mz.mean,ppm = ppm)

      is.match <- match_mz_foverlaps( xcms.net$mz.diff,is.mass.diff$mass_diff,
                                ppm.base = xcms.net$mz.mean,ppm = ppm)



      cn.net <- cbind(xcms.net[cn.match$ion1,],cn.match[,c("mz.ppm","ion1") ],
                      cn.mass.diff[cn.match$ion2,])[mass_diff > 0]
      ad.net <- cbind(xcms.net[ad.match$ion1,],ad.match[,c("mz.ppm","ion1") ],
                      ad.mass.diff[ad.match$ion2,])[mass_diff > 0]
      is.net <- cbind(xcms.net[is.match$ion1,],is.match[,c("mz.ppm","ion1") ],
                      is.mass.diff[is.match$ion2,])[mass_diff > 0]

    }

    ### filter CN label pattern
    {
      cn.net <- cn.net%>%
        dplyr::mutate(pave_pattern = paste0("C",C_count ,"N",N_count  ))
      cn.net.list <- split(cn.net,cn.net$from)
      prefilt <- sapply(cn.net.list,function(x.cn){ 0 %in% x.cn$N_count   })
      cn.net.list <- cn.net.list[prefilt]
      message_with_time("Find CN label pattern...")
      cn.net.list.hit <- bplapply(names(cn.net.list),function(x){

        #message(x)
        x.cn <- cn.net.list[[x]]
        possible.c.count <- unique(x.cn$C_count)
        possible.n.count <- unique(x.cn$N_count)
        c.max <- x.cn$from.mz[1]/14
        possible.c.count <- possible.c.count[possible.c.count < c.max&possible.c.count > 0]
        cn.comb <- expand.grid(C = possible.c.count,
                               N = possible.n.count,
                               p.cor = NA)

        if (!nrow(cn.comb)) return(NULL)
        for (i.cn in 1:nrow(cn.comb)) {
          this.c <- cn.comb$C[i.cn]
          this.n <- cn.comb$N[i.cn]
          all.form <- c(paste0("C0N",this.n,""),paste0("C",this.c,"N0"),paste0("C",this.c,"N",this.n,""))
          all.form <- setdiff(all.form,"C0N0")
          if (!all(all.form %in% x.cn$pave_pattern ) ) next

          #message(x)
          to.id <- x.cn$to[match(all.form,x.cn$pave_pattern)]
          m.detected <- xcms.val[c(x.cn$from[1],to.id),  xcms.pave.sample$sampleNames]
          colnames(m.detected) <- xcms.pave.sample$sample.type
          rownames(m.detected) <- c("C0N0",all.form)
          m.detected <- m.detected/m.detected[1,1]
          m.ideal <- get_ideal_CN_ratio(this.c,this.n)%>%t
          m.ideal <- m.ideal[rownames(m.detected),colnames(m.detected)]
          p.cor <- cor(as.vector(m.detected),as.vector(m.ideal))
          cn.comb$p.cor[i.cn] <- p.cor
        }

        p.cor.max <- max(cn.comb$p.cor,na.rm = T)
        #message(p.cor.max)
        if(p.cor.max< 0.75) return(NULL)
        cn.comb <- cn.comb%>%dplyr::slice_max(p.cor)
        all.form <- c(paste0("C0N",cn.comb$N,""),paste0("C",cn.comb$C,"N0"),
                      paste0("C",cn.comb$C,"N",cn.comb$N,""))
        all.form <- setdiff(all.form,"C0N0")
        x.cn <- x.cn[match(all.form,x.cn$pave_pattern),]
        x.cn$pave_cor <- p.cor.max
        x.cn$pave_formula <-  paste0("C",cn.comb$C,"N",cn.comb$N,"")
        return(x.cn)

      },BPPARAM = SerialParam(progressbar = T))
      names(cn.net.list.hit) <- names(cn.net.list)
      cn.net.list.hit <- cn.net.list.hit[!sapply(cn.net.list.hit,is.null)]
      cn.net.hit <- data.table::rbindlist(cn.net.list.hit)
    }


    ### RT and mz error evaluation
    if(T){
      cn.net.eval <- cn.net%>%
        dplyr::mutate(cn.hit = ion1%in% cn.net.hit$ion1)%>%
        dplyr::arrange(cn.hit)
      #cn.net.eval <- cn.net.eval[1:1000000,]
      cols <- c("TRUE" = "red","FALSE" = "#888888")
      p <- ggplot() +
        geom_point(data = cn.net.eval,
                   aes(x = mz.ppm, y = rt.diff,
                       col = cn.hit),
                   pch = 16,alpha = 0.5,size = 0.2)+
        scale_color_manual(values = cols)+
        labs(x = "mz error (ppm)",y = "rt shift (s)")+
        coord_fixed(ppm/rt.tol)+
        theme_bw()+
        theme(legend.position = "none")
      p.r <- ggplot(cn.net.eval)+
        geom_histogram(aes(y = rt.diff,x = after_stat(density), fill = cn.hit),
                       position = "dodge",#stat = "density",
                       binwidth = 1,col = "white")+
        scale_fill_manual(values =cols)+
        scale_x_continuous(expand = c(0,0))+
        labs(x = NULL, y = NULL)+
        theme_classic()+
        theme(axis.text.y = element_blank(),
              legend.position = "none",
              axis.ticks = element_blank())

      p.u <- ggplot(cn.net.eval)+
        geom_histogram(aes(x = mz.ppm,y = after_stat(density), fill = cn.hit),
                       position = "dodge",#stat = "density",
                       bins = 10,col = "white")+
        scale_fill_manual(values = cols)+
        scale_y_continuous(expand = c(0,0))+
        labs(x = NULL, y = NULL,fill = "CN labeled")+
        theme_classic()+
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank())
      #p.u
      p.all <- p.u+plot_spacer()+p+p.r+
        plot_layout(heights  = c(0.2,0.8),widths = c(0.8,0.2),guides = "collect")
      open_plot_win(p.all)



      ### filter with dynamic error range
      {
        ppm <- quantile(cn.net.hit$mz.ppm,0.95)
        rt.rol <- quantile(cn.net.hit$rt.diff,0.95)

        cn.net.hit <- cn.net.hit[mz.ppm < ppm&rt.diff < rt.tol]
        ad.net <- ad.net[mz.ppm < ppm&rt.diff < rt.tol]
        is.net <- is.net[mz.ppm < ppm&rt.diff < rt.tol]
      }


    }



    ### integration
    {

      xcms.net.candidate <- bind_rows(ad.net,is.net,cn.net.hit)
      xcms.net.candidate <- split(xcms.net.candidate,xcms.net.candidate$ion1)
      xcms.net.matched <- xcms.net[as.numeric(names(xcms.net.candidate)),]

    }


    ### ig
    {
      xcms.net.matched$label <- sapply(
        xcms.net.candidate,
        function(x){
          x%>%
            dplyr::mutate(
              temp = case_when(
                type == "isotope"~ element,
                type == "CN_label" ~ pave_pattern,
                type == "adduct" ~ paste0( adduct.from ," to ", adduct.to)
              ),
              label = paste0(type,": ",temp))%>%
            dplyr::pull(label)%>%
            paste0(collapse = "\n")
        })
      xcms.net.matched$is_CN <- sapply(
        xcms.net.candidate,
        function(x){
          "CN_label" %in%  x$type
        })
      xcms.ig <- igraph::graph_from_data_frame(xcms.net.matched)
      node.group <- igraph::components(xcms.ig)$membership

      vda <- vdata(xcms.ig)%>%
        dplyr::mutate(color = case_when(
          name %in% cn.net.hit$from ~ "#E64B35",
          T~"#97C2FC"
        ))
      vda -> vdata(xcms.ig)
      ngcn <- data.frame(
        ig = unique(node.group),
        cn = 0,
        e = 0
      )

      for (i in unique(node.group)) {

        sub.ig <- igraph_filter_vertex(xcms.ig , node.group==i)
        eda <- edata(sub.ig)
        cn.count  <- sum(eda$is_CN)

        visNetwork::visIgraph(sub.ig)
        ngcn$cn[i] <- cn.count
        ngcn$e[i] <- nrow(eda)
        ngcn$cne[i] <- cn.count/nrow(eda)

      }




    }

    x.to.h <- ad.net%>%
      dplyr::filter( from %in% cn.net.hit$from,
                     to %in% cn.net.hit$from)
    xcms.ig%>%
      igraph_filter_distance(c("3162","3202"),1)%>%
      visNetwork::visIgraph()

  }







}



PAVE2_find_xcms_CN <- function(){





}



get_xcms_feature_connect <- function(xcms.xcms,rt.tol = 5){


  xcms.fdf <- featureDefinitions(xcms.xcms)
  xcms.net <- expand.grid(
    from = 1:nrow(xcms.fdf),
    to = 1:nrow(xcms.fdf)
  )
  xcms.net <- data.table::as.data.table(xcms.net)

  {
    xcms.net <- xcms.net[from < to ][
      , rt.diff := abs(xcms.fdf$rtmed[to]-xcms.fdf$rtmed[from]) ][
        rt.diff < rt.tol,][
          , c("from.mz","to.mz") := .( xcms.fdf$mzmed[from], xcms.fdf$mzmed[to])][
            ,c("mz.diff","mz.mean") := .(to.mz-from.mz,(from.mz+to.mz)/2)]
  }


  return(xcms.net)



}



chemform_simplify <- function(chemform){

  ele.matrix <- MSCC::chemform_parse(chemform)
  MSCC:::chemform_from_ele_matrix(ele.matrix)

}
