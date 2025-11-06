# Tue Oct 14 16:08:40 2025 PAVE------------------------------
{

  library(xcms)


  xcms.xcms <- msdev.pave.demo@xcmsData$NegativeMS1
  xcms.peaks <- chromPeaks(xcms.xcms)%>%as.data.frame()
  xcms.regf <- groupChromPeaks(
    xcms.xcms,
    param = PeakDensityParam(
      minFraction = 0.8,
      sampleGroups = pData(xcms.xcms)$sample.type,
      binSize = 0.002,
      ppm = 10
    )
  )

  featureDefinitions(xcms.regf)%>%nrow()

}
# Thu Oct 16 19:28:44 2025 ------------------------------
{
  files <- dir("d:/data/2025.10.10.PVAE/data7600/rawdata/",full.names = T)
  files.new <- sub(" ",x = files,"")
  file.rename(files,files.new)

  fdf <- featureDefinitions(object@xcmsData$NegativeMS1)
  sum(fdf$pave_cor>0,na.rm = T)
  sum(fdf$pave_seed!="",na.rm = T)
  nrow(fdf)


  chemform_adduct("C10H16N5O13P3","+")
  chemform_adduct("[13]C10H16[15]N5O13P3","+")
  chemform_adduct("[13]C10H16[15]N5O13P3","+") + c(-3:3)*0.00631988

  chemform_adduct("C11H20[15]N4O11P2","+")
  chemform_adduct("[13]C11H20N4O11P2","+")
  chemform_adduct("[13]C11H20[15]N4O11P2","+")



}

# Tue Oct 21 17:13:05 2025 ------------------------------
{

  cfm <- CFM_predict()
  cfm.fg <- CFM_fraggen()
  cfm.ano <- CFM_annotate()
  cfm.ano <- CFM_annotate_by_predict()
  MSdev:::CFM


}
# Thu Oct 23 18:41:56 2025 ------------------------------
{

  library(gt)
  library(dplyr)

  # Example data (simplified)
  df <- tribble(
    ~Category, ~Type, ~`>1E3_neg`, ~`>1E4_neg`, ~`>1E5_neg`, ~`>1E6_neg`, ~`>1E3_pos`, ~`>1E4_pos`, ~`>1E5_pos`, ~`>1E6_pos`,
    "ATOMCOUNT", "Peaks in procedure blank", 10147, 8979, 3051, 447, 20600, 16273, 6975, 1024,
    "ATOMCOUNT", "Other peaks without labeling", 3422, 2097, 320, 35, 6312, 3583, 722, 79,
    "ATOMCOUNT", "Labeling but ρ < 0.75", 186, 153, 33, 5, 339, 247, 81, 13,
    "ATOMCOUNT", "Logical labeling (i.e. biological)", 1151, 998, 363, 81, 1510, 1217, 397, 75,
    "JUNKREMOVER", "Isotopes", 404, 346, 89, 14, 467, 360, 108, 22,
    "JUNKREMOVER", "Dimer or double charge", 45, 40, 16, 3, 25, 19, 7, 0,
    "JUNKREMOVER", "Adducts (same polarity mode)", 278, 234, 81, 10, 431, 367, 102, 9,
    "JUNKREMOVER", "Adducts (opposite polarity mode)", 11, 9, 2, 0, 55, 43, 19, 6
  )

  gt_tbl <- df %>%
    gt(rowname_col = "Type", groupname_col = "Category") %>%
    tab_spanner(label = "Negative mode", columns = 3:6) %>%
    tab_spanner(label = "Positive mode", columns = 7:10) %>%
    cols_label(
      `>1E3_neg` = ">1E3",
      `>1E4_neg` = ">1E4",
      `>1E5_neg` = ">1E5",
      `>1E6_neg` = ">1E6",
      `>1E3_pos` = ">1E3",
      `>1E4_pos` = ">1E4",
      `>1E5_pos` = ">1E5",
      `>1E6_pos` = ">1E6"
    ) %>%
    tab_header(
      title = md("**Table S3. Peak annotation by PAVE as a function of peak height (*E. coli*)**"),
      subtitle = md("*E. coli*")
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_title(groups = "title")
    ) %>%
    fmt_number(columns = 3:10, decimals = 0)

  gt_tbl

  gtsave(gt_tbl, "D:/temp/a.pdf")

}

# Mon Oct 27 16:09:48 2025 ------------------------------
{


  x <- "C6H12O6"
  x.adduct <- chemform_adduct(x,MSCC::adduct.table$Name,value = "all")


  data("adduct.table")
  data("pave_adduct")

  pave.adduct.missing <- MSCC::chemform_adduct_check(
    pave_adduct$Adduct )

}

# Tue Oct 28 13:56:17 2025 ------------------------------
{

  library(data.table)

  match_mz_dt <- function(obs, theo, ppm = 10) {
    dt_obs <- data.table(
      mz = obs,
      idx_obs = seq_along(obs),
      key = "mz"
    )
    dt_theo <- data.table(
      mz = theo,
      idx_theo = seq_along(theo),
      key = "mz"
    )

    # Add tolerance bounds
    dt_obs[, `:=`(
      lower = mz * (1 - ppm * 1e-6),
      upper = mz * (1 + ppm * 1e-6)
    )]

    # Set keys for interval join
    setkey(dt_theo, mz)
    setkey(dt_obs, lower, upper)

    # Overlap join
    overlaps <- foverlaps(dt_theo, dt_obs, type = "within", nomatch = NULL)

    # Calculate ppm
    overlaps[, ppm_error := abs(mz - i.mz) / mz * 1e6]
    overlaps <- overlaps[ppm_error <= ppm]

    result <- overlaps[, .(obs_idx = idx_obs, theo_idx = idx_theo,
                           obs_mz = i.mz, theo_mz = mz, ppm_error)]

    return(result)
  }

  match_mz(xcms.fdf$mzmed,adduct.diff$mass_diff)->a


}

# Tue Oct 28 14:16:17 2025 ------------------------------
{
  data <- edit_df_in_excel(rowname = T)

  plot.data <- data%>%
    tibble::column_to_rownames("X1")
  plot.data <- scale(t(plot.data))%>%t
  ComplexHeatmap::Heatmap(
    plot.data,
    name = "Relative\nAbundance",
    col = circlize::colorRamp2(breaks = c(-2,0,2),
                               colors = c("#1C73AE","white","#E06C55")),
    row_names_side  = "left",
    column_split = paste0("Time",rep(1:3,each = 3)),
    show_row_dend  = F,
    show_column_names = F,
    top_annotation = columnAnnotation(
      df = data.frame(
        Group = paste0("Time",rep(1:3,each = 3))
      ),
      col = list(Group = c("Time1" = "#E64B35",
                           "Time2" = "#4DBBD5",
                           "Time3" = "#00A087"))
    )
  )-> hm

  open_plot_win(hm,8,5)

  p <- plot_PCA(t(plot.data),
           pca.group = paste0("Time",rep(1:3,each = 3)))
  open_plot_win(p,3,3)
}

# Tue Oct 28 15:19:14 2025 ------------------------------
{

  x <- c(100,mz.range.ppm(100,3),200)%>%sort()
  y <- c(50,mz.range.ppm(100,2),300)%>%sort()
  join(
    x,
    y,
    tolerance = 0,
    ppm = 10
  )


  ggplot(xcms.net.matched)+
    #geom_density2d(aes(x = mz.ppm, y = abs(rt.diff) ) )
    geom_point(aes(x = mz.ppm, y = abs(rt.diff) ),alpha = 0.05 )


}

# Tue Oct 28 17:16:29 2025 ------------------------------
{

  match_ppm_sorted <- function(x, y, ppm_tol = 10) {
    x <- sort(x)
    y <- sort(y)

    i <- 1
    j <- 1
    nx <- length(x)
    ny <- length(y)

    matches <- vector("list", nx)

    while (i <= nx && j <= ny) {
      ppm_diff <- abs((y[j] - x[i]) / x[i] * 1e6)

      if (ppm_diff <= ppm_tol) {
        # 找到匹配，存下所有符合ppm范围的y[j]
        k <- j
        tmp <- c()
        while (k <= ny && abs((y[k] - x[i]) / x[i] * 1e6) <= ppm_tol) {
          tmp <- c(tmp, k)
          k <- k + 1
        }
        matches[[i]] <- tmp
        i <- i + 1
      } else if (y[j] < x[i]) {
        j <- j + 1
      } else {
        i <- i + 1
      }
    }

    matches
  }


  x <- sort(runif(1e6, 100, 1000))
  y <- sort(runif(1e6, 100, 1000))
  system.time({
    idx <- match_ppm_sorted(x, y, ppm_tol = 10)
  })


  match_ppm_all <- function(x, y, ppm_tol = 10) {

    x <- sort(x)
    y <- sort(y)

    # 每个x对应的匹配上下限
    xlow <- x * (1 - ppm_tol / 1e6)
    xhigh <- x * (1 + ppm_tol / 1e6)

    # 在y中找到对应范围的索引边界
    left_idx <- findInterval(xlow, y)
    right_idx <- findInterval(xhigh, y)

    res_list <- vector("list", length(x))

    for (i in seq_along(x)) {
      message(i)
      if (left_idx[i] < right_idx[i]) {
        # 匹配的y索引范围
        idx_range <- seq.int(left_idx[i] + 1, right_idx[i])
        # 计算实际ppm误差
        ppm_diff <- abs((y[idx_range] - x[i]) / x[i] * 1e6)
        res_list[[i]] <- data.frame(
          x = x[i],
          y = y[idx_range],
          ppm_diff = ppm_diff
        )
      }
    }

    # 合并所有结果
    res <- do.call(rbind, res_list)
    rownames(res) <- NULL
    res
  }


  set.seed(1)
  x <- sort(runif(1e6, 100, 101))
  y <- sort(runif(3e3, 100, 101))

  res <- match_ppm_all(x, y, ppm_tol = 50)  # 50 ppm
  print(res)

}

# Tue Oct 28 17:22:45 2025 ------------------------------
{
  library(data.table)


  match_ppm_foverlaps <- function(x, y,ppm.base = x, ppm = 10 ) {

    dt_x <- data.table(x = x,
                       xid = seq_along(x))
    dt_x[, `:=`(
      xmin = x - ppm.base * ppm / 1e6,
      xmax = x + ppm.base * ppm / 1e6
    )]
    setkey(dt_x, xmin, xmax)

    dt_y <- data.table(y = y,
                       yid = seq_along(y),
                       y_start = y, y_end = y)
    setkey(dt_y, y_start, y_end)

    res <- foverlaps(dt_y, dt_x, by.x = c("y_start", "y_end"),
                     type = "within", nomatch = 0L)


    res[, ppm := abs((y - x) / ppm.base[xid] * 1e6)]



    res[]
  }



  set.seed(123)
  x <- sort(runif(1e6, 50, 900))
  y <- sort(runif(1e3, 50, 900))

  system.time(resx <- match_mz_foverlaps(x, y, ppm = 10))
  system.time(resy <- match_mz_foverlaps(y, x, ppm = 10))
  print(res)


  ggplot(eda)+
    geom_density_2d(aes(x = mz.ppm , y  = abs(rt.diff)))

  plot_density(eda$mz.ppm)
  plot_density(eda$rt.diff)

}

# Wed Oct 29 11:48:47 2025 ------------------------------
{

  f.raw <- dir("d:/data/2025.10.10.PVAE/data7600/20251029_TOFMS_neg/",recursive = T,full.names = T)
  f.new <- paste0(dirname(f.raw ),"/neg_",basename(f.raw))
  file.rename(f.raw,f.new)


  f.raw <- dir("d:/data/2025.10.10.PVAE/data7600/rawdata/",recursive = T,full.names = T)
  f.new <- gsub(pattern = " ",x = f.raw ,replacement = "_" )
  file.rename(f.raw,f.new)

  msdev.8600.dda <- MSdev("d:/data/2025.10.10.PVAE/data7600/rawdata/")
  msdev.8600.dda <- MSdev_msConvert(msdev.8600.dda)
  msdev.8600.dda <- MSdev_checkSampleInfo(msdev.8600.dda)
  MSdev_save(msdev.8600.dda)

}

# Wed Oct 29 13:46:50 2025 ------------------------------
{
  library(data.table)
  DT = as.data.table(iris)

  # FROM[WHERE, SELECT, GROUP BY]
  # DT  [i,     j,      by]

  DT[Petal.Width > 1.0, mean(Petal.Length), by = Species]
  #      Species       V1
  #1: versicolor 4.362791
  #2:  virginica 5.552000

  DT = data.table(
    ID = c("b","b","b","a","a","c"),
    a = 1:6,
    b = 7:12,
    c = 13:18
  )
  DT
  #        ID     a     b     c
  #    <char> <int> <int> <int>
  # 1:      b     1     7    13
  # 2:      b     2     8    14
  # 3:      b     3     9    15
  # 4:      a     4    10    16
  # 5:      a     5    11    17
  # 6:      c     6    12    18
  class(DT$ID)
  # [1] "character"
}


# Wed Oct 29 15:38:51 2025 ------------------------------
{

  a <- xcms.net.matched%>%
    dplyr::filter(mz.ppm < 5,rt.diff < 5)%>%
    dplyr::mutate(label =  paste0(from.adduct," to ",to.adduct,"\n",
                                  "ppm=",str_digit(mz.ppm),
                                  ",rt=",str_digit(rt.diff)))

  aig <- igraph::graph_from_data_frame(a )
  visIgraph(igraph_filter_distance(aig,from = "24068",1))%>%
    visEdges(smooth = T)

  av <- vdata(aig)

  for (i.av in av$name) {

    av.ig <- igraph_filter_distance(aig,from = i.av,1)
    av.ed <- edata(av.ig)

    av.edf  <- av.ed %>%
      dplyr::filter(from == i.av)
    av.edt <- av.ed%>%
      dplyr::filter(to == i.av)
    av.candi <- c(av.edf$from.adduct,av.edt$to.adduct)


  }


}


# Wed Oct 29 18:50:33 2025 ------------------------------
{
  xcms.xcms <- msdev.pave.480@xcmsData$PositiveMS1

  cn.peaks.high.cor <- cn.peaks%>%
    dplyr::filter(pave_cor >= 0.75)%>%
    dplyr::mutate(
      isotope = grepl("isotope",pave_junkremover),
      adduct = grepl("adduct",pave_junkremover),
      LowC = grepl("LowC",pave_junkremover),
      opposite_adduct  = grepl("opposite_adduct",pave_junkremover),
      dimer   = grepl("dimer",pave_junkremover),
      ringing    = grepl("ringing",pave_junkremover),
      others = grepl("ringing|dimer",pave_junkremover)
    )

  ggvenn::ggvenn(
    list(isotope = which(cn.peaks.high.cor$isotope),
         adduct = which(cn.peaks.high.cor$adduct),
         LowC = which(cn.peaks.high.cor$LowC),
         others = which(cn.peaks.high.cor$other)),
    stroke_color = "white",
    fill_color = ggsci::pal_npg()(5)
  )->p
  open_plot_win(p,8,4)


  ggplot(xcms.net.matched[sample(1:nrow(xcms.net.matched),nrow(xcms.net.matched)*0.1),],
         aes(mz.ppm,rt.diff))+
    geom_density_2d_filled()+
    geom_point(alpha = 0.01)

  p1 <- ggplot(xcms.net.matched)+
    geom_density(aes(x = mz.ppm))+
    theme_classic()

  p2 <- ggplot(xcms.net.matched)+
    geom_density(aes(x = rt.diff))+
    theme_classic()

  open_plot_win(p1+p2,6,2)

}

# Thu Oct 30 13:58:52 2025 WJY chromosome data------------------------------
{


  library(Biostrings)

  fa.files <- dir("d:/temp/FTP/",pattern = ".fa.gz",full.names = T)
  seqs <- readDNAStringSet(fa.files)

  base_counts <- letterFrequency(seqs,
                                 letters = c("A", "T", "C", "G", "N"))

  res <- data.frame(
    name = names(seqs),
    base_counts
  )


  openxlsx::write.xlsx(res,file = "d:/temp/chrom.ATCGN.count.xlsx")

}


# Thu Oct 30 14:35:09 2025 ------------------------------
{

  xcms.xcms <- msdev.pave.480@xcmsData$PositiveMS1

  mz.grouped <- groupMz(abs(xcms.net$mz.diff),return.type = "data.frame",ppm = 5)
  a <- mz.grouped%>%
    dplyr::group_by( mz.group)%>%
    dplyr::mutate(count = n())%>%
    dplyr::ungroup()

  b <- a%>%
    dplyr::distinct( mz.group,.keep_all = T)%>%
    dplyr::arrange(-count)%>%
    dplyr::mutate(ratio = count / nrow(mz.grouped))
  head(b)
  c <- dplyr::filter(b,ratio > 1e-3)
  match_mz(c$mz,adduct.diff$mass_diff,mz.ppm = 20)%>%is.na()%>%`!`%>%sum(.)
  match_mz(c$mz,pave_adduct$mass_diff,mz.ppm = 20)%>%is.na()%>%`!`%>%sum(.)

}

# Fri Oct 31 15:04:06 2025 ------------------------------
{

  object <- msdev.pave.480
  PAVE_get_atom_count()
  PAVE_junk_remover()


  ft.pave.high.cor <- cn.list%>%
    rbindlist()%>%
    dplyr::filter(pave_cor > 0.75)

  sum(duplicated(ft.pave.high.cor$feature_id))
  ft.dup <- ft.pave.high.cor %>%
    dplyr::filter(feature_id %in%ft.pave.high.cor$feature_id[duplicated(ft.pave.high.cor$feature_id)])


}

# Fri Oct 31 17:31:45 2025 ------------------------------
{
  system.time({
    xcms.net <- expand.grid(
      from = 1:nrow(xcms.fdf),
      to = 1:nrow(xcms.fdf)
    )
    xcms.net <- data.table::as.data.table(xcms.net)
    xcms.net <- xcms.net[from != to ][
      , rt.diff := abs(xcms.fdf$rtmed[to]-xcms.fdf$rtmed[from]) ][
        rt.diff < rt.tol,]
  })

  #   user  system elapsed
  #  21.57    2.51   21.30

 system.time({
   a <- lapply(FUN = function(from ){
     to <- 1:nrow(xcms.fdf)
     to.rt <- xcms.fdf$rtmed[to]
     from.rt <- xcms.fdf$rtmed[from]
     rt.diff <- abs(to.rt- from.rt)
     idx <- rt.diff < rt.tol
     data.table(from = from,
                to = to[idx],
                rt.diff = rt.diff[idx])
   }, 1:nrow(xcms.fdf))
 })


}


# Mon Nov  3 15:19:54 2025 ------------------------------
{

 a <- MSCC::elem_table%>%
    dplyr::filter(abundance > 0.01 &abundance <1)


 dt_x[,c("xmin","xmax") := .(
     x - ppm.base * ppm / 1e6,
   x + ppm.base * ppm / 1e6
 )]


 adduct.diff <- expand.grid(
   adduct.from = 1:nrow(adduct.table),
   adduct.to = 1:nrow(adduct.table)
 )%>%
   dplyr::filter(
     adduct.from != adduct.to,
     adduct.table$m_c[adduct.from] == adduct.table$m_c[adduct.to]
   )%>%
   dplyr::mutate(
     chemform_diff = chemform_calc(adduct.table$Formula_diff[adduct.to],
                                   adduct.table$Formula_diff[adduct.from],
                                   calc = "-",return = "chemform"),
     mass_diffc = chemform_mz(chemform_diff),
     mass_diff =adduct.table$Mass[adduct.to] - adduct.table$Mass[adduct.from],
     #charge = adduct.table$Charge[adduct.to],
     adduct.from =adduct.table$Adduct[adduct.from],
     adduct.to = adduct.table$Adduct[adduct.to],
     x = mass_diff-mass_diffc
   )


 chemform_adduct("C6H12O6","[M+Cl]-") - chemform_adduct("C6H12O6","[M-H2O-H]-")

}


# Tue Nov  4 13:38:08 2025 ------------------------------
{

  kegg.cp%>%
    as.data.frame()%>%
    dplyr::filter(Exact_mass < 1200)%>%
    dplyr::pull(Formula)%>%
    get_formula_ele_count("N")->kegg.n.count

  hist(kegg.n.count)
  quantile(kegg.n.count,c(0,0.9,0.99,0.999,0.9999))

}
