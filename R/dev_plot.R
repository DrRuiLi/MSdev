#' @title capture_base_plot
#' @description
#' Evaluate base graphics in a null device and return a \code{recordedplot}
#' object for replay or export (e.g. with \code{\link{open_plot_win}}).
#'
#' @param expr Base graphics commands, typically passed as a braced expression.
#' @param envir Environment in which \code{expr} is evaluated.
#'
#' @return A \code{recordedplot} object (from \code{\link[grDevices]{recordPlot}}).
#' @export
#'

capture_base_plot <- function(expr, envir = parent.frame()) {
  grDevices::pdf(nullfile())
  on.exit(grDevices::dev.off(), add = TRUE)
  grDevices::dev.control(displaylist = "enable")
  eval(substitute(expr), envir = envir)
  grDevices::recordPlot()
}


#' @title open_plot_win
#' @description
#' create a temp.png file and open in Windows
#'
#'
#' @param p ggplot/Complexheatmap
#' @param width num
#' @param height num
#'
#' @return null
#' @export
#'

open_plot_win <- function(p,width = 5,height = 4){

  temp.file <- tempfile(fileext = ".png")
  if (any(c("Heatmap","HeatmapList")%in%class(p))) {
    export::graph2png(ComplexHeatmap::draw(p,padding = unit(c(1, 1, 1, 3), "mm")),
                     file =temp.file,
                     width = width,height= height
                     )
  }else if(any(c("ggplot")%in%class(p))){
    ggplot2::ggsave(filename = temp.file,plot = p,
                    width = width,height= height,dpi = 600)
  }else{
    export::graph2png(p,
                      file =temp.file,
                      width = width,height= height
    )
  }
  open_file(temp.file)

}


open_plot_pdf <- function(p,width = 5,height = 4){

  temp.file <- tempfile(fileext = ".pdf")
  if (any(c("Heatmap","HeatmapList")%in%class(p))) {
    export::graph2png(ComplexHeatmap::draw(p),
                      file =temp.file,
                      width = width,height= height
    )
  }else if(any(c("ggplot")%in%class(p))){
    ggplot2::ggsave(filename = temp.file,plot = p,
                    width = width,height= height,dpi = 600)
  }else{
    export::graph2pdf(p,
                      file =temp.file,
                      width = width,height= height
    )
  }
  open_file(temp.file)

}


#' @title open_plot_ppt
#' @description
#' create a temp.pptx file and open in Windows
#'
#' @param p ggplot/Complexheatmap/base plot/recordedplot
#' @param width num
#' @param height num
#'
#' @return null
#' @export
#'
open_plot_ppt <- function(p,width = 5,height = 4){

  temp.file <- tempfile(fileext = ".pptx")
  if (any(c("Heatmap","HeatmapList")%in%class(p))) {
    export::graph2ppt(ComplexHeatmap::draw(p,padding = unit(c(1, 1, 1, 3), "mm")),
                      file =temp.file,
                      width = width,height= height
    )
  }else if(any(c("ggplot")%in%class(p))){
    export::graph2ppt(p,
                      file =temp.file,
                      width = width,height= height
    )
  }else{
    export::graph2ppt(p,
                      file =temp.file,
                      width = width,height= height
    )
  }
  open_file(temp.file)

}



#' @title ggplot_sum_patchwork
#' @description
#' add all ggplot by patchwork
#'
#' @param ggplot.list a list with all item as ggplot objective
#'
#' @return null
#' @export
#'

ggplot_sum_patchwork <- function(ggplot.list){
  x <- ggplot.list
  x.len <- length(x)
  sum.exp <- 1
  x.exp <- paste0( "x.sum <- ",paste0(paste0("x[[",1:x.len,"]]"),collapse = " + "),
                   "+patchwork::plot_annotation(tag_levels=\"A\")")%>%
    str2expression()
  eval(x.exp)
  return(x.sum)
}



ggplot_ggsci <- function(pal = "npg"){

  ggsci_db <- ggsci:::ggsci_db
  ggsci_pal <- names(ggsci_db)
  pal <- match.arg(pal,ggsci_pal)
  pal.col <- ggsci_db[[pal]][[1]]
  message("Show color palette from ggsci: ",
          crayon::magenta(pal))
  scales::show_col(pal.col)

  return(pal.col)

}

colored_text <- function(x , color = "#E64B35"){

  crayon::make_style(color)(x)


}



#' export plot to pdf file
#' pdf() and export::graph2pdf() not support `append` arg
#' using qpdf::pdf_combine() to realize that function
#'
#' @title Export Graph2pdf
#' @description Graph2pdf.
#' @param p ggplot
#' @param file_path file path
#' @param append logic
#' @param ... additional arguments passed to export functions
#'
#' @return null
#' @export
#'

export_graph2pdf <- function(p ,
                             file_path ,
                             append = F,
                             ...
                             ){
  dir.exists(dirname(file_path))
  if("list" %in% class(p)){
    if((!append) & file.exists(file_path)){
      pdf(file_path)
      dev.off()
    }
    for(i in seq_along(p)){

      export_graph2pdf(p[[i]],file_path = file_path ,append = T,...)


    }
    return(invisible())


  }

  file_path <- normalizePath(file_path,mustWork = F)
  tpf1 <- paste0(tempfile(),".pdf")
  tpf2 <- paste0(tempfile(),".pdf")
  suppressMessages(export::graph2pdf(plot_multi_formate(p),
                                     file = tpf1,...))

  if (append & file.exists(file_path)) {
    qpdf::pdf_combine(input = c(file_path,tpf1),
                      output = tpf2)
  }else{
    tpf2 <- tpf1
  }
  message("Exported graph as ",file_path)
  file.copy(tpf2,file_path,overwrite = T)
  suppressWarnings(file.remove(tpf1,tpf2))
  return(invisible())

}

plot_multi_formate <- function(p){


  if ("Heatmap" %in% class(p)) {

    ComplexHeatmap::draw(p)

  }else{

    p
  }


}



ggplot_from_img <- function(img_path,coord = F,...){

#
  #p.jpg <- readJPEG(img_path,native = T)
  #p.frame <- ggplot()+
  #  xlim(c(0,10))+
  #  ylim(c(0,10))+
  #  theme_void()
  #p <- p.frame + inset_element(p = p.jpg,
  #                             position[1],
  #                             position[2],
  #                             position[3],
  #                             position[4]
  #                             )
  #return(p)


  p <- ggplot() +
    #xlim(c(0,10))+
    #ylim(c(0,10))+
    ggimage::geom_image(aes(x=0,y=0,image = img_path),...)+
    theme_void()
  if (coord) {
    p <- p+theme_classic()
  }

  return(p)

}



ggplot_km <- function(km.km,legend_tile = "group",
                      legend_label = names(km.km$strata) ){

  #km.km <- surv_fit(Surv(time, status) ~ 1, data = colon)
  km.sum <- summary(km.km)
  km.pval <- survminer::surv_pvalue(km.km)
  plot.data <- data.frame(
    time = km.sum$time,
    strata= km.sum$strata,
    surv=km.sum$surv,
    lower=km.sum$lower,
    upper=km.sum$upper
  )

  ggplot(plot.data)+
    geom_line(aes(x = time , y = surv,
                  col = strata))+
    geom_ribbon(aes(x = time ,ymin = lower,
                    ymax = upper,fill = strata),alpha = 0.3)+
    annotate(geom="text",x = quantile(range(plot.data$time),0.2),
             y = 0.2,
             label =paste0("p = ", format(km.pval[2],scientific = F,digits =3)))+
    scale_color_manual(values = ggsci::pal_npg()(8),
                       labels = legend_label)+
    scale_fill_manual(values = ggsci::pal_npg()(8),
                      labels = legend_label)+
    ylim(c(0,1))+
    labs(x = "Time",y = "Survival Probability",
         col = legend_tile,
         fill = legend_tile)+
    theme_bw()

}


colramp<- function(breaks = c(0,0.5,1),
                   colors = c("white","#F7844F","#B20C26"),
                   na.col = "#AAAAAA",
                   ...){

  circlize::colorRamp2(breaks =breaks,
                       colors = colors,
                       ...)
}


make_group_color <- function(x,palette = "random",verbose=F){

  palette <- match.arg(palette,c(names(ggsci:::ggsci_db),"random"))
  x <- unique(x)%>%sort()
  if (palette == "random") {
    col  <- randomcoloR::distinctColorPalette(length(x))
  }else{
    col <- ggsci:::ggsci_db[[palette]][[1]][1:length(x)]

  }
  if (verbose) {
    message(paste0(mapply(x,col,FUN = colored_text),collapse = " "))

  }
  names(col)<-x
  return(col)
}


scale_color_random <- function(...) {
  ggplot2::discrete_scale("colour", "random", function(n) randomcoloR::distinctColorPalette(n), ...)
}



scale_fill_random <- function(...) {
  ggplot2::discrete_scale("fill", "random", function(n) randomcoloR::distinctColorPalette(n), ...)
}

# Function to format legend labels with subscripts, with custom color and labels
scale_subscript_legend <- function(values = NULL, labels = NULL) {

  # Helper function to format labels with subscripts
  format_labels <- function(labels) {
    lapply(labels, function(label) {
      # Check if label contains a number
      if (grepl("\\d", label)) {
        # If number found, split into text and number
        split_label <- strsplit(label, "(?<=\\D)(?=\\d)", perl = TRUE)[[1]]
        # Create an expression with subscript
        expression_text <- bquote(.(split_label[1])[.(split_label[2])])
      } else {
        # If no number, return label as is
        expression_text <- label
      }
      return(expression_text)
    })
  }

  # Create a custom scale for color and fill
  scale_custom_labels <- list(
    scale_color_manual(values = values, labels = format_labels(labels)),
    scale_fill_manual(values = values, labels = format_labels(labels))
  )

  # Return the list of scale components to add to the plot
  return(scale_custom_labels)
}



heatmap_set_size <- function(hm,width  = 5,height= 5){

  hm@heatmap_param$width <- unit(width,"cm")
  hm@heatmap_param$height <- unit(height,"cm")

  return(hm)
}


get_ggplot_from_heatmap <- function(hm){
  wrap_elements(grid::grid.grabExpr(draw(hm)))
}

ggplot_irange <- function(IR, scale = 1e-6){

  ir.data <- as.data.frame(IR)
  ir.data <- (ir.data *  scale)%>%
    dplyr::mutate(no = 1:n())

  ggplot(ir.data)+
    geom_segment( aes(x = start,xend = end, y = no,yend = no) )

}


#' Mirror EIC comparison grid for feature groups
#'
#' Pairwise chromatogram comparison for features in one or more xcms feature
#' groups on an \code{XcmsExperiment} (or compatible) object. Panel \eqn{(i,j)}
#' shows the column feature EIC above \eqn{y = 0} and the row feature EIC below
#' (each cropped to \code{rtmed +/- rt_half} and normalized to its own max).
#' Cell labels are EIC similarities from
#' \code{MsExperiment::otherData(xcms)$EIC_Similarity} when available (e.g.
#' after \code{\link{xcms_group_feature_EIC}}). Strip backgrounds are colored
#' by feature group when \pkg{ggh4x} is installed.
#'
#' @param xcms An \code{XcmsExperiment} / \code{MsExperiment} (or
#'   \code{XCMSnExp}) with feature definitions and \code{featureGroups}.
#' @param feature_group Character vector of feature group id(s) to compare.
#' @param chroms Optional chromatograms for the selected features (rownames =
#'   feature ids). If \code{NULL}, extracted via
#'   \code{\link{get_xcms_feature_chromatogram}}.
#' @param rt_half Half-width (seconds) of the RT window around each feature
#'   \code{rtmed} (default \code{20}).
#' @param max_features Maximum features kept per group (ordered by
#'   \code{rtmed}). \code{NULL} keeps all members (default).
#' @param sample_index Sample index into
#'   \code{otherData(xcms)$EIC_Similarity} when several samples are stored
#'   (default \code{1L}).
#' @param title Optional plot title. Default is
#'   \code{"Feature group: <id> (<rt center>), ..."}, where the RT center is the
#'   median \code{rtmed} of each group's member features.
#' @return A \code{ggplot} object.
#' @export
plot_xcms_feature_group_EIC_comparasion <- function(xcms,
                                               feature_group,
                                               chroms = NULL,
                                               rt_half = 20,
                                               max_features = NULL,
                                               sample_index = 1L,
                                               title = NULL) {
  if (!(inherits(xcms, "XcmsExperiment") ||
        inherits(xcms, "MsExperiment") ||
        inherits(xcms, "XCMSnExp"))) {
    stop("`xcms` must be an XcmsExperiment, MsExperiment, or XCMSnExp")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_xcms_feature_group_EIC_comparasion()")
  }
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required for plot_xcms_feature_group_EIC_comparasion()")
  }
  if (missing(feature_group) || !length(feature_group)) {
    stop("`feature_group` must be a non-empty character vector")
  }
  feature_group <- as.character(feature_group)
  feature_group <- feature_group[!is.na(feature_group) & nzchar(feature_group)]
  if (!length(feature_group)) {
    stop("`feature_group` must be a non-empty character vector")
  }
  if (!is.numeric(rt_half) || length(rt_half) != 1L || rt_half <= 0) {
    stop("`rt_half` must be a positive numeric(1)")
  }

  fdf <- as.data.frame(xcms::featureDefinitions(xcms))
  fg_all <- as.character(xcms::featureGroups(xcms))
  if ("feature_id" %in% names(fdf)) {
    fids_all <- as.character(fdf$feature_id)
  } else {
    fids_all <- rownames(fdf)
  }
  rt_all <- as.numeric(fdf$rtmed)
  mz_all <- as.numeric(fdf$mzmed)
  if ("feature_group_rt" %in% names(fdf)) {
    fgrt_all <- as.numeric(fdf$feature_group_rt)
  } else {
    fgrt_all <- rep(NA_real_, length(fids_all))
  }
  names(rt_all) <- names(mz_all) <- names(fg_all) <- names(fgrt_all) <- fids_all

  miss <- setdiff(feature_group, unique(stats::na.omit(fg_all)))
  if (length(miss)) {
    stop(
      "Unknown feature_group: ", paste(miss, collapse = ", "),
      ". Available examples: ",
      paste(utils::head(unique(stats::na.omit(fg_all)), 5L), collapse = ", ")
    )
  }

  # Preserve caller order of groups; features within group by rtmed
  sel_fg <- feature_group
  gfids_list <- lapply(sel_fg, function(g) {
    gfids <- fids_all[fg_all == g]
    gfids <- gfids[order(rt_all[gfids], na.last = TRUE)]
    if (!is.null(max_features) && length(gfids) > as.integer(max_features)) {
      gfids <- gfids[seq_len(as.integer(max_features))]
    }
    gfids
  })
  names(gfids_list) <- sel_fg
  sel_fids <- unlist(gfids_list, use.names = FALSE)
  if (length(sel_fids) < 2L) {
    stop("Need at least 2 features across the selected feature_group(s)")
  }
  n <- length(sel_fids)
  sel_fg_of <- fg_all[sel_fids]

  sim_mat <- NULL
  if (requireNamespace("MsExperiment", quietly = TRUE) &&
      (inherits(xcms, "MsExperiment") || inherits(xcms, "XcmsExperiment"))) {
    od <- tryCatch(MsExperiment::otherData(xcms), error = function(e) NULL)
    sim_store <- if (!is.null(od)) od$EIC_Similarity else NULL
    if (!is.null(sim_store) && length(sim_store)) {
      si <- as.integer(sample_index)[[1]]
      if (si < 1L || si > length(sim_store)) si <- 1L
      sim_mat <- sim_store[[si]]
    }
  }

  if (is.null(chroms) || identical(chroms, NA)) {
    message("Extracting chromatograms for ", length(sel_fids), " features...")
    bpp <- if (requireNamespace("BiocParallel", quietly = TRUE)) {
      BiocParallel::SerialParam()
    } else {
      NULL
    }
    chrom_args <- list(
      xcms.xcms = xcms,
      feature.id = sel_fids,
      selected_sample = "maxo",
      rt = "expand",
      expandRt = rt_half,
      aggregationFun = "max",
      attachPeaks = FALSE
    )
    if (!is.null(bpp)) chrom_args$BPPARAM <- bpp
    chroms <- do.call(get_xcms_feature_chromatogram, chrom_args)
  }

  .eic_df <- function(fid) {
    rn <- rownames(chroms)
    ii <- match(fid, rn)
    rt_med <- rt_all[[fid]]
    if (is.na(ii) || !is.finite(rt_med)) {
      return(NULL)
    }
    ch <- chroms[ii, 1L]
    if (!requireNamespace("MSnbase", quietly = TRUE)) {
      stop("Package 'MSnbase' is required to read chromatogram rtime/intensity")
    }
    rt <- as.numeric(MSnbase::rtime(ch))
    inten <- as.numeric(MSnbase::intensity(ch))
    keep <- is.finite(rt) & rt >= (rt_med - rt_half) & rt <= (rt_med + rt_half)
    if (!any(keep)) {
      return(NULL)
    }
    data.frame(rt = rt[keep], intensity = inten[keep], stringsAsFactors = FALSE)
  }

  eic_cache <- lapply(sel_fids, .eic_df)
  names(eic_cache) <- sel_fids

  mz_lab <- sprintf("%.3f", mz_all[sel_fids])
  names(mz_lab) <- sel_fids
  if (anyDuplicated(mz_lab)) {
    mz_lab <- make.unique(mz_lab, sep = "_")
    names(mz_lab) <- sel_fids
  }
  feat_lab <- sprintf("%s\n%s", sel_fids, mz_lab)
  names(feat_lab) <- sel_fids
  lab_levels <- unname(feat_lab[sel_fids])

  fg_pal <- grDevices::hcl.colors(max(3L, length(sel_fg)), "Dark 3")[seq_along(sel_fg)]
  names(fg_pal) <- sel_fg
  strip_fills <- unname(fg_pal[sel_fg_of])
  strip_bg <- lapply(strip_fills, function(col) {
    ggplot2::element_rect(fill = col, colour = "grey40", linewidth = 0.2)
  })

  .sim_ij <- function(fi, fj) {
    if (is.null(sim_mat)) return(NA_real_)
    if (identical(fi, fj)) return(1)
    ii <- match(fi, rownames(sim_mat))
    jj <- match(fj, colnames(sim_mat))
    if (is.na(ii) || is.na(jj)) return(NA_real_)
    as.numeric(sim_mat[ii, jj])
  }

  parts <- vector("list", n * n)
  labs_list <- vector("list", n * n)
  k <- 0L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      k <- k + 1L
      fi <- sel_fids[[i]]
      fj <- sel_fids[[j]]
      d_col <- eic_cache[[fj]]
      d_row <- eic_cache[[fi]]
      if (is.null(d_col) || is.null(d_row)) next
      mc <- max(d_col$intensity, na.rm = TRUE)
      mr <- max(d_row$intensity, na.rm = TRUE)
      if (!is.finite(mc) || mc <= 0) mc <- 1
      if (!is.finite(mr) || mr <= 0) mr <- 1
      s <- .sim_ij(fi, fj)
      s_lab <- if (is.finite(s)) sprintf("%.2f", s) else "NA"
      parts[[k]] <- rbind(
        data.frame(
          rt = d_col$rt,
          intensity = d_col$intensity / mc,
          role = "column",
          feature_group = unname(sel_fg_of[[fj]]),
          row_lab = feat_lab[[fi]],
          col_lab = feat_lab[[fj]],
          stringsAsFactors = FALSE
        ),
        data.frame(
          rt = d_row$rt,
          intensity = -d_row$intensity / mr,
          role = "row",
          feature_group = unname(sel_fg_of[[fi]]),
          row_lab = feat_lab[[fi]],
          col_lab = feat_lab[[fj]],
          stringsAsFactors = FALSE
        )
      )
      labs_list[[k]] <- data.frame(
        rt = mean(range(c(d_col$rt, d_row$rt), na.rm = TRUE)),
        intensity = 0.95,
        row_lab = feat_lab[[fi]],
        col_lab = feat_lab[[fj]],
        sim_lab = s_lab,
        stringsAsFactors = FALSE
      )
    }
  }

  plot_df <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
  lab_df <- do.call(rbind, labs_list[!vapply(labs_list, is.null, logical(1))])
  if (is.null(plot_df) || !nrow(plot_df)) {
    stop("No chromatogram data available for the selected feature_group(s)")
  }
  plot_df$row_lab <- factor(plot_df$row_lab, levels = lab_levels)
  plot_df$col_lab <- factor(plot_df$col_lab, levels = lab_levels)
  plot_df$feature_group <- factor(plot_df$feature_group, levels = sel_fg)
  lab_df$row_lab <- factor(lab_df$row_lab, levels = lab_levels)
  lab_df$col_lab <- factor(lab_df$col_lab, levels = lab_levels)

  if (is.null(title)) {
    fg_rt <- vapply(sel_fg, function(g) {
      gfids <- gfids_list[[g]]
      as.numeric(fgrt_all[gfids[[1L]]])
    }, numeric(1))
    fg_lab <- ifelse(
      is.finite(fg_rt),
      sprintf("%s (%s)", sel_fg, round(fg_rt)),
      sel_fg
    )
    title <- sprintf(
      "Feature group: %s",
      paste(fg_lab, collapse = ", ")
    )
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$rt,
      y = .data$intensity,
      color = .data$feature_group,
      fill = .data$feature_group,
      group = .data$role
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, color = "grey40") +
    ggplot2::geom_area(alpha = 0.35, position = "identity", linewidth = 0) +
    ggplot2::geom_line(linewidth = 0.3, na.rm = TRUE) +
    ggplot2::geom_text(
      data = lab_df,
      ggplot2::aes(x = .data$rt, y = .data$intensity, label = .data$sim_lab),
      inherit.aes = FALSE,
      size = 1.6,
      color = "black",
      fontface = "bold"
    ) +
    ggplot2::scale_color_manual(name = "feature group", values = fg_pal) +
    ggplot2::scale_fill_manual(name = "feature group", values = fg_pal) +
    ggplot2::coord_cartesian(ylim = c(-1.05, 1.05)) +
    ggplot2::labs(
      x = "Retention time (s)",
      y = "Normalized intensity (mirror)",
      title = title
    ) +
    ggplot2::theme_bw(base_size = 6) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 3.2, angle = 0, lineheight = 0.85),
      strip.text.y = ggplot2::element_text(size = 3.2, angle = -90, lineheight = 0.85),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.05, "lines"),
      plot.title = ggplot2::element_text(size = 9, face = "bold")
    )

  if (requireNamespace("ggh4x", quietly = TRUE)) {
    p <- p + ggh4x::facet_grid2(
      rows = ggplot2::vars(row_lab),
      cols = ggplot2::vars(col_lab),
      scales = "free_x",
      strip = ggh4x::strip_themed(
        background_x = strip_bg,
        background_y = strip_bg,
        by_layer_x = FALSE,
        by_layer_y = FALSE
      )
    )
  } else {
    message("Package 'ggh4x' not installed; strip fills by feature group disabled.")
    p <- p + ggplot2::facet_grid(
      rows = ggplot2::vars(row_lab),
      cols = ggplot2::vars(col_lab),
      scales = "free_x"
    )
  }

  attr(p, "feature_group") <- sel_fg
  attr(p, "feature_id") <- sel_fids
  p
}


#' Mirror plot of two chromatograms
#'
#' Overlay two chromatograms on a shared RT axis: \code{chrom1} above
#' \eqn{y = 0} and \code{chrom2} mirrored below. Intensities are optionally
#' normalized to each chromatogram's own maximum.
#'
#' @param chrom1,chrom2 A \code{Chromatogram} / \code{XChromatogram}, a
#'   one-cell \code{MChromatograms} / \code{XChromatograms}, or a
#'   \code{data.frame} with columns \code{rt} and \code{intensity}.
#' @param labels Character length-2 labels for the legend (default
#'   \code{c("chrom1", "chrom2")}).
#' @param normalize Logical; if \code{TRUE} (default), scale each trace to
#'   its own max intensity.
#' @param colors Optional length-2 color vector for the two traces.
#' @param title Optional plot title.
#' @return A \code{ggplot} object.
#' @export
plot_Chromatograph_mirror <- function(chrom1,
                                      chrom2,
                                      labels = c("chrom1", "chrom2"),
                                      normalize = TRUE,
                                      colors = NULL,
                                      title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_Chromatograph_mirror()")
  }
  labels <- as.character(labels)
  if (length(labels) < 2L) {
    labels <- c(labels, paste0(labels[[1]], "_2"))[seq_len(2L)]
  }
  labels <- labels[seq_len(2L)]

  .chrom_to_df <- function(chrom, lab) {
    if (is.null(chrom)) {
      stop("Chromatogram input cannot be NULL")
    }
    if (is.data.frame(chrom)) {
      nm <- names(chrom)
      if (!all(c("rt", "intensity") %in% nm)) {
        stop("data.frame chrom must have columns `rt` and `intensity`")
      }
      df <- data.frame(
        rt = as.numeric(chrom$rt),
        intensity = as.numeric(chrom$intensity),
        stringsAsFactors = FALSE
      )
    } else if (inherits(chrom, c("MChromatograms", "XChromatograms"))) {
      if (!requireNamespace("MSnbase", quietly = TRUE)) {
        stop("Package 'MSnbase' is required to read chromatogram objects")
      }
      ch <- chrom[1L, 1L]
      df <- data.frame(
        rt = as.numeric(MSnbase::rtime(ch)),
        intensity = as.numeric(MSnbase::intensity(ch)),
        stringsAsFactors = FALSE
      )
    } else if (inherits(chrom, c("Chromatogram", "XChromatogram"))) {
      if (!requireNamespace("MSnbase", quietly = TRUE)) {
        stop("Package 'MSnbase' is required to read chromatogram objects")
      }
      df <- data.frame(
        rt = as.numeric(MSnbase::rtime(chrom)),
        intensity = as.numeric(MSnbase::intensity(chrom)),
        stringsAsFactors = FALSE
      )
    } else {
      stop(
        "Unsupported chrom type: ", paste(class(chrom), collapse = "/"),
        ". Expect Chromatogram, MChromatograms, or data.frame(rt, intensity)."
      )
    }
    keep <- is.finite(df$rt) & is.finite(df$intensity)
    df <- df[keep, , drop = FALSE]
    if (!nrow(df)) {
      stop("No finite rt/intensity points in chromatogram: ", lab)
    }
    df$label <- lab
    df
  }

  d1 <- .chrom_to_df(chrom1, labels[[1]])
  d2 <- .chrom_to_df(chrom2, labels[[2]])

  if (isTRUE(normalize)) {
    m1 <- max(d1$intensity, na.rm = TRUE)
    m2 <- max(d2$intensity, na.rm = TRUE)
    if (!is.finite(m1) || m1 <= 0) m1 <- 1
    if (!is.finite(m2) || m2 <= 0) m2 <- 1
    d1$intensity <- d1$intensity / m1
    d2$intensity <- d2$intensity / m2
  }
  d2$intensity <- -d2$intensity

  plot_df <- rbind(d1, d2)
  plot_df$label <- factor(plot_df$label, levels = labels)

  if (is.null(colors)) {
    colors <- grDevices::hcl.colors(3L, "Dark 3")[c(1L, 3L)]
  }
  if (length(colors) < 2L) {
    colors <- rep(colors, length.out = 2L)
  }
  colors <- unname(colors[seq_len(2L)])
  names(colors) <- labels

  if (is.null(title)) {
    title <- sprintf("Mirror chromatogram: %s vs %s", labels[[1]], labels[[2]])
  }
  y_lab <- if (isTRUE(normalize)) {
    "Normalized intensity (mirror)"
  } else {
    "Intensity (mirror)"
  }

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$rt,
      y = .data$intensity,
      color = .data$label,
      fill = .data$label,
      group = .data$label
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
    ggplot2::geom_area(alpha = 0.35, position = "identity", linewidth = 0) +
    ggplot2::geom_line(linewidth = 0.45, na.rm = TRUE) +
    ggplot2::scale_color_manual(name = NULL, values = colors) +
    ggplot2::scale_fill_manual(name = NULL, values = colors) +
    ggplot2::labs(
      x = "Retention time (s)",
      y = y_lab,
      title = title
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = 11, face = "bold")
    )
}
