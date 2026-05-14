#' @title Plot isotopologue ratios (circular bar plot)
#' @description
#' Plot isotopologue ratio assay from a \code{MSIPIsotopologueData} object as a
#' circular stacked bar plot. Bars are colored by isotopologue.
#'
#' @param object A \code{MSIPIsotopologueData} (inherits \code{SummarizedExperiment}).
#' @param assay Character. Which assay to plot. Default \code{"ratio"}; if missing,
#'   falls back to \code{"ratio.negative"} then \code{"ratio.positive"}.
#' @param isotopologue_label Character. RowData column used as fill label when present.
#'   Default \code{"label.isotopologue"}, then \code{"isotopologue_form"}, then rownames.
#' @param sample_order Optional character vector of sample.source order.
#' @param isotopologue_order Optional character vector of isotopologue order (fill levels).
#' @param min_ratio Numeric threshold for isotopologue average ratio. Isotopologues
#'   with average ratio (across all samples) below \code{min_ratio} are grouped
#'   into one bar labeled \code{"other"}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_MSIPIsotopologueData_ratio <- function(object,
                                           assay = "ratio",
                                           isotopologue_label = c("label.isotopologue", "isotopologue_form"),
                                           sample_order = NULL,
                                           isotopologue_order = NULL,
                                           min_ratio = 0.01) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  if (is.null(object) || !methods::is(object, "SummarizedExperiment")) {
    stop("object must be a MSIPIsotopologueData (SummarizedExperiment).")
  }

  available <- names(SummarizedExperiment::assays(object))
  assay_pick <- assay
  if (!(assay_pick %in% available)) {
    assay_pick <- intersect(c("ratio.negative", "ratio.positive"), available)[1]
  }
  if (is.na(assay_pick) || is.null(assay_pick) || !(assay_pick %in% available)) {
    stop("No ratio assay found in object. Expected one of: ratio, ratio.negative, ratio.positive.")
  }

  mat <- SummarizedExperiment::assay(object, assay_pick)
  if (is.null(mat) || !nrow(mat) || !ncol(mat)) {
    stop("Selected assay has no data.")
  }

  rda <- tryCatch(as.data.frame(SummarizedExperiment::rowData(object)), error = function(e) NULL)
  if (is.null(rda) || !nrow(rda)) {
    rda <- data.frame(row_id = rownames(mat) %||% paste0("row_", seq_len(nrow(mat))),
                      stringsAsFactors = FALSE)
    rownames(rda) <- rda$row_id
  }

  # Ensure rownames exist for mapping.
  row_id <- rownames(rda)
  if (is.null(row_id) || any(!nzchar(row_id))) {
    if ("isotopologue_id" %in% colnames(rda)) {
      row_id <- as.character(rda$isotopologue_id)
      row_id[is.na(row_id) | !nzchar(row_id)] <- paste0("row_", seq_len(nrow(rda)))[is.na(row_id) | !nzchar(row_id)]
    } else {
      row_id <- paste0("row_", seq_len(nrow(rda)))
    }
    rownames(rda) <- row_id
  }

  # Align rda rows with matrix rows.
  if (is.null(rownames(mat)) || any(!nzchar(rownames(mat)))) {
    if (nrow(mat) == nrow(rda)) {
      rownames(mat) <- rownames(rda)
    } else {
      rownames(mat) <- paste0("row_", seq_len(nrow(mat)))
    }
  }
  if (!all(rownames(mat) %in% rownames(rda))) {
    # fall back to positional alignment if needed
    if (nrow(mat) == nrow(rda)) rownames(mat) <- rownames(rda)
  }

  # Pick isotopologue label column.
  isotopologue_label <- as.character(isotopologue_label)
  fill_col <- intersect(isotopologue_label, colnames(rda))[1]
  if (is.na(fill_col) || is.null(fill_col)) {
    fill_col <- if ("isotopologue_id" %in% colnames(rda)) "isotopologue_id" else NULL
  }
  fill_val <- if (!is.null(fill_col)) as.character(rda[rownames(mat), fill_col]) else rownames(mat)
  fill_val[is.na(fill_val) | !nzchar(fill_val)] <- rownames(mat)[is.na(fill_val) | !nzchar(fill_val)]
  .as_mplus <- function(x) {
    x <- as.character(x)
    out <- x
    has_mplus <- grepl("M\\+[0-9]+", x)
    out[has_mplus] <- sub(".*(M\\+[0-9]+).*", "\\1", x[has_mplus])
    miss <- !has_mplus
    if (any(miss)) {
      num <- suppressWarnings(as.integer(gsub(".*?([0-9]+)$", "\\1", x[miss])))
      ok <- !is.na(num)
      idx <- which(miss)[ok]
      out[idx] <- paste0("M+", num[ok])
    }
    out
  }
  fill_val <- .as_mplus(fill_val)

  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  df$row_id <- rownames(mat)
  df$isotopologue <- fill_val

  df_long <- tidyr::pivot_longer(
    df,
    cols = setdiff(colnames(df), c("row_id", "isotopologue")),
    names_to = "sample.source",
    values_to = "ratio"
  )
  df_long$ratio <- as.numeric(df_long$ratio)
  df_long$ratio[is.na(df_long$ratio)] <- 0

  # Group low-average isotopologues into "other".
  min_ratio <- suppressWarnings(as.numeric(min_ratio))
  if (!is.finite(min_ratio) || is.na(min_ratio)) min_ratio <- 0
  if (min_ratio > 0) {
    iso_avg <- stats::aggregate(ratio ~ isotopologue, data = df_long, FUN = mean, na.rm = TRUE)
    low_iso <- as.character(iso_avg$isotopologue[is.finite(iso_avg$ratio) & iso_avg$ratio < min_ratio])
    if (length(low_iso)) {
      df_long$isotopologue <- as.character(df_long$isotopologue)
      df_long$isotopologue[df_long$isotopologue %in% low_iso] <- "other"
    }
  }

  # Factor orders.
  if (!is.null(sample_order)) {
    df_long$sample.source <- factor(df_long$sample.source, levels = sample_order)
  } else {
    df_long$sample.source <- factor(df_long$sample.source, levels = unique(df_long$sample.source))
  }
  if (!is.null(isotopologue_order)) {
    if ("other" %in% unique(as.character(df_long$isotopologue)) && !("other" %in% isotopologue_order)) {
      isotopologue_order <- c(as.character(isotopologue_order), "other")
    }
    df_long$isotopologue <- factor(df_long$isotopologue, levels = isotopologue_order)
  } else {
    df_long$isotopologue <- factor(df_long$isotopologue, levels = unique(df_long$isotopologue))
  }
  # Keep "other" at the bottom of legend (last factor level).
  if ("other" %in% levels(df_long$isotopologue)) {
    lev <- levels(df_long$isotopologue)
    lev <- c(setdiff(lev, "other"), "other")
    df_long$isotopologue <- factor(as.character(df_long$isotopologue), levels = lev)
  }

  # Colors: keep NPG palette for isotopologues; force "other" to grey.
  iso_levels <- levels(df_long$isotopologue)
  iso_non_other <- iso_levels[iso_levels != "other"]
  base_cols <- if (length(iso_non_other)) ggsci::pal_npg("nrc")(length(iso_non_other)) else character(0)
  fill_vals <- stats::setNames(base_cols, iso_non_other)
  if ("other" %in% iso_levels) fill_vals <- c(fill_vals, other = "grey")

  # Derive compound name for the global title.
  compound_name <- NULL
  if (!is.null(rda) && nrow(rda) && "compound_name" %in% colnames(rda)) {
    cn <- unique(stats::na.omit(as.character(rda$compound_name)))
    if (length(cn)) compound_name <- cn[[1]]
  }
  if (is.null(compound_name) && !is.null(rda) && nrow(rda) && "compound_id" %in% colnames(rda)) {
    cid <- unique(stats::na.omit(as.character(rda$compound_id)))
    if (length(cid)) compound_name <- cid[[1]]
  }
  if (is.null(compound_name)) compound_name <- ""

  .pick_first <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA)
    x[[1]]
  }
  .fmt_num <- function(x, digits = 3) {
    if (is.na(x) || !is.finite(x)) return("NA")
    format(signif(x, digits = digits), scientific = FALSE, trim = TRUE)
  }
  .fmt_sci <- function(x, digits = 2) {
    if (is.na(x) || !is.finite(x)) return("NA")
    formatC(x, format = "e", digits = digits)
  }

  formula_txt <- NA_character_
  if (!is.null(rda) && nrow(rda) && "formula" %in% colnames(rda)) {
    formula_txt <- .pick_first(as.character(rda$formula))
  }
  rt_txt <- NA_real_
  if (!is.null(rda) && nrow(rda) && "rt" %in% colnames(rda)) {
    rt_txt <- .pick_first(as.numeric(rda$rt)/60)
  }
  intensity_txt <- NA_real_
  if (!is.null(rda) && nrow(rda) && "intensity" %in% colnames(rda)) {
    intensity_txt <- .pick_first(as.numeric(rda$intensity))
  }
  subtitle_txt <- paste0(
    "rt: ", .fmt_num(rt_txt),
    " | intensity: ", .fmt_sci(intensity_txt),
    " | formula: ", ifelse(is.na(formula_txt) || !nzchar(formula_txt), "NA", formula_txt)
  )

  # Build one circular bar plot per sample.source, then patchwork them together.
  samples <- levels(df_long$sample.source)
  make_one <- function(ss) {
    d <- df_long[df_long$sample.source %in% ss, , drop = FALSE]
    # x axis is isotopologue (angle), y axis is ratio (radius)
    ggplot2::ggplot(d, ggplot2::aes(x = 1, y = .data$ratio, fill = .data$isotopologue)) +
      ggplot2::geom_col(width = 1, color = "black", linewidth = 0.25) +
      ggplot2::coord_polar(start = 0,theta = "y",direction = -1) +
      ggplot2::scale_fill_manual(values = fill_vals, breaks = iso_levels, drop = FALSE) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
        legend.position = "right"
      ) +
      ggplot2::labs(title = as.character(ss), fill = "Isotopologue")
  }
  plots <- lapply(samples, make_one)
  names(plots) <- as.character(samples)

  patchwork::wrap_plots(plots) +
    patchwork::plot_layout(guides = "collect")  +
    patchwork::plot_annotation(title = compound_name, subtitle = subtitle_txt)
}


#' @title Report isotopologue ratio as a PDF
#' @description
#' Create a PDF report that contains isotopologue ratio circular plots for each
#' compound in \code{object@advancedAna$MSIP$isotopologue_data}.
#'
#' @param object A \code{MSdev} object.
#' @param file Output pdf path. Default:
#'   \code{file.path(object@projectInfo$projectDir, "report", "MSIP_report_isotopologue_ratio.pdf")}.
#' @param assay Assay name passed to \code{\link{plot_MSIPIsotopologueData_ratio}}.
#' @param ... Additional arguments passed to \code{\link{get_MSIPIsotopologueData}} when
#'   isotopologue data needs to be built.
#'
#' @return Invisibly returns \code{file}.
#' @export
MSIP_report_isotopologue_ratio <- function(object,
                                          file = NULL,
                                          assay = "ratio",
                                          ...) {
  message_with_time("MSIP_report_isotopologue_ratio: checking input object")
  if (is.null(object) || !methods::is(object, "MSdev")) {
    stop("object must be a MSdev object.")
  }

  message_with_time("MSIP_report_isotopologue_ratio: resolving output file")
  if (is.null(file)) {
    proj_dir <- tryCatch(object@projectInfo$projectDir, error = function(e) NULL)
    if (is.null(proj_dir) || is.na(proj_dir) || !nzchar(proj_dir)) {
      stop("object@projectInfo$projectDir is missing; please provide `file` explicitly.")
    }
    file <- file.path(proj_dir, "report", "MSIP_report_isotopologue_ratio.pdf")
  }
  message_with_time("MSIP_report_isotopologue_ratio: creating output directory")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  message_with_time("MSIP_report_isotopologue_ratio: retrieving isotopologue_data")
  iso.list <- tryCatch(object@advancedAna$MSIP$isotopologue_data, error = function(e) NULL)
  if (is.null(iso.list) || !is.list(iso.list) || !length(iso.list)) {
    message_with_time("MSIP_report_isotopologue_ratio: isotopologue_data missing, building with get_MSIPIsotopologueData()")
    iso.list <- get_MSIPIsotopologueData(object, ...)
  }
  if (!is.list(iso.list) || !length(iso.list)) {
    stop("No isotopologue_data found/built for this MSdev object.")
  }

  message_with_time("MSIP_report_isotopologue_ratio: start, total compounds = ", length(iso.list))
  message_with_time("Output file: ", normalizePath(file, winslash = "/", mustWork = FALSE))

  message_with_time("MSIP_report_isotopologue_ratio: use export_graph2pdf() for report writing")

  i <- 0L
  n <- length(iso.list)
  append_pdf <- FALSE
  n_written <- 0L
  for (nm in names(iso.list)) {
    i <- i + 1L
    sei <- iso.list[[nm]]
    if (is.null(sei) || !methods::is(sei, "SummarizedExperiment")) {
      message_with_time("[", i, "/", n, "] skip ", nm, " (not SummarizedExperiment)")
      next
    }
    message_with_time("[", i, "/", n, "] plotting ", nm, " (build plot)")
    p <- plot_MSIPIsotopologueData_ratio(sei, assay = assay)
    message_with_time("[", i, "/", n, "] plotting ", nm, " (append to pdf)")
    export_graph2pdf(p, file_path = file, append = append_pdf)
    append_pdf <- TRUE
    n_written <- n_written + 1L
  }

  if (n_written == 0L) {
    message_with_time("MSIP_report_isotopologue_ratio: no valid plot written")
  }
  message_with_time("MSIP_report_isotopologue_ratio: done")
  invisible(file)
}

