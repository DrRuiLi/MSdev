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

# Thu Nov  6 19:30:15 2025 ------------------------------
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
  xcms.net.matched <- xcms.net.matched%>%
    dplyr::filter(is_CN)
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


  cn.ig <- igraph::graph_from_data_frame( cn.net.hit  )



  idx <- 1:1e5
  a <- matrixSub(rt[idx],rt[idx])

}


{

  rt <- xcms.fdf$rtmed [1:1000]

  x <- data.table(point = rt,x = rt)
  range.x <- x[,.(id,start =  x - 5, end = x + 5)]
  setkey(x,point)
  setkey(range.x,start,end)
  a <- foverlaps(x, range.x, type="any", which=TRUE)



  object@projectInfo$MSdevFile

}
# Fri Nov  7 18:32:52 2025 ------------------------------
{

  p.r <- ggplot(cn.net.eval)+
    geom_histogram(aes(y = rt.diff,x = after_stat(ndensity), fill = cn.hit),
                   position = "dodge",#stat = "density",
                   binwidth = 1,col = "white")+
    stat_ecdf(aes(y = rt.diff, col = cn.hit),linewidth = 2)+
    scale_fill_manual(values =cols)+
    scale_color_manual(values =cols)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)#,limits = c(0,10)
                       )+
    labs(x = NULL, y = NULL)+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          legend.position = "none",
          axis.ticks = element_blank())

  ggplot(cn.net.eval)+
    stat_ecdf(aes(x = rt.diff,y = after_stat(density)))

  ggplot(cn.net.eval)+
    geom_bar(aes(x= 0,y=0, fill = cn.hit),stat = "identity")+
    scale_fill_manual(values =cols)+
    labs(x = NULL, y = NULL,fill = "CN labeled")+
    theme_void()+
    theme(legend.position = "inside")

  x.to.h <- ad.net%>%
    dplyr::filter( from %in% cn.net.hit$from,
                   to %in% cn.net.hit$from)
  xcms.ig%>%
    igraph_filter_distance(c("3162","3202"),1)%>%
    visNetwork::visIgraph()
}

# Mon Nov 10 13:25:26 2025 ------------------------------
{

  x <- cn.net.hit$mz.ppm
  ec <- ecdf(mz.ppm)
  cdf <- data.frame(x = knots(ec), y = ec(knots(ec)))

  # Find inflection point
  inflection::findiplist(cdf$x, cdf$y,index = 0)
  plot(ec,xlim = c(-5,5))
  abline(v = -1.0455205    )


  library(segmented)
  x <- cdf$x
  y <- cdf$y
  model <- lm(y ~ x)
  seg <- segmented(model, seg.Z = ~x, psi = 5)  # initial guess
  summary(seg)

  mu <- median(x)
  sigma <- mad(x)
  x[abs(x - mu) < 3 * sigma]
  mu-3*sigma
  x_clean <- x[abs(x - mu) < 3 * sigma]

  mu2 <- mean(x_clean)
  sigma2 <- sd(x_clean)
  ci95 <- c(mu2 - 3 * sigma2, mu2 + 3 * sigma2)
  ci95
}

# Mon Nov 10 19:20:14 2025 ------------------------------
{

  vda <- vdata(ig)%>%
    dplyr::mutate(id = name, .before = name)
  visNetwork::visNetwork(nodes = vda,edges = edata(ig))%>%
    visNetwork::visOptions(width = "100%",
                           height = "100%")



}

# Tue Nov 11 10:20:45 2025 ------------------------------
{
  ###
  node.count <- vda%>%
    dplyr::distinct(node.group,n)

  cn.seed.net <- xcms.net.matched%>%
    dplyr::filter(as.character(from) %in% cn.seed,
                  as.character(to) %in% cn.seed)%>%
    dplyr::mutate(from.cn = cn.seed.formula[as.character(from)],
                  to.cn = cn.seed.formula[as.character(to)],
                  chemfrom.diff = chemform_calc(from.cn ,to.cn,"-",return = "chemform"))


  x1 <- x[, .SD[order(type != "CN_label")[1]], by = chemform_diff]
  x1

  ### TEMP
  {


    xcms.net.candidate <- split(xcms.net.candidate,xcms.net.candidate$ion1)
    xcms.net.matched <- xcms.net[as.numeric(names(xcms.net.candidate)),]
    names(xcms.net.candidate) <- xcms.net.matched$eid


    xcms.net.integrated <- bplapply(
      xcms.net.candidate,
      function(x){

        x <- x[,temp := fcase(
          type == "isotope",chemform_diff,
          type == "CN_label" , pave_pattern,
          type == "adduct" , chemform_diff,
          default = ""
        )][,label :=  paste0(type,": ",temp)
        ][, .SD[order(type != "CN_label")[1]], by = chemform_diff]

        nrow(x)
      },BPPARAM = SerialParam(progressbar = T))


    xcms.ig <- igraph::graph_from_data_frame(xcms.net.matched)

    cn.ig <- igraph_filter_vertex(xcms.ig,  as.character(cn.net.hit$from))
    node.group <- igraph::components(cn.ig)$membership
    vda <-  vdata(cn.ig)%>%
      dplyr::mutate(node.group = node.group ,
                    pave_formula = cn.seed.formula[name],
                    label = pave_formula)%>%
      dplyr::group_by(node.group)%>%
      dplyr::mutate(n = n())%>%
      dplyr::ungroup()
    eda <- edata(cn.ig)
    vdata(cn.ig) <- vda

    igraph_filter_vertex(cn.ig, node.group== 276 )%>%
      vis_igraph()%>%
      visEdges(arrows = "to")

    igraph_filter_distance(xcms.ig,"928",10)%>%
      vis_igraph()%>%
      visEdges(arrows = "to")


    ### ig
    {

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

  }


  obj <- MSdev_load("d:/data/2025.10.10.PVAE/PAVE_With_Params/QEplus_ppm10_sn10.rdata")
  PAVE2(obj)

  obj <- MSdev_load("d:/data/2025.10.10.PVAE/data8600/MSdev_2025_10_16.Rdata")
  PAVE2(obj)

  obj <- MSdev_load("d:/data/2025.10.10.PVAE/PAVE_With_Params/TOF8600_ppm10_sn10.rdata")
  PAVE2(obj)

  obj <- MSdev_load("d:/data/2025.10.10.PVAE/PAVE_With_Params/TOF8600_ppm25_sn100.rdata")
  PAVE2(obj)

  obj <- MSdev_load("d:/data/2025.10.10.PVAE/PAVE_With_Params/TOF8600_DDA_ppm10_sn10.rdata")
  PAVE2(obj)

  obj <- MSdev_load("d:/data/2025.10.10.PVAE/PAVE_With_Params/TOF7600_ppm10_sn10.rdata")
  PAVE2(obj)

  obj <- MSdev_load("d:/data/2025.10.10.PVAE/PAVE_With_Params/OE480_480k_ppm10_sn10.rdata")
  PAVE2(obj)

  obj <- MSdev_load("d:/data/2025.10.10.PVAE/PAVE_With_Params/OE480_120k_ppm10_sn10.rdata")
  PAVE2(obj)

  obj <- MSdev_load("d:/data/2025.10.10.PVAE/PAVE_With_Params/OE480_60k_ppm10_sn10.rdata")
  PAVE2(obj)

}

# Wed Nov 19 16:24:10 2025 ------------------------------
{

  for (iv in i.cn.seed.vda$name) {

    for (jv in i.cn.seed.vda$name) {

      ijp <- igraph::all_simple_paths(
        i.cn.seed.ig,
        from = iv,to = jv,mode  = "all")
      if (length(ijp) >1) {

        i.path <- get_edges_from_epath(ig = i.cn.seed.ig, ijp[[1]],directed = F)

      }

    }

  }




}



# Thu Nov 20 15:37:25 2025 ------------------------------
{

  ### list all path i to j
  ### test if chemform cum exist

  vis_igraph(i.cn.seed.ig)

  for(iv in seq_len(length(V(i.cn.seed.ig))) ){

    ap <- igraph::shortest_paths(i.cn.seed.ig,
                                 from = iv,
                                 output = "both",
                                 mode = "all")

    for (jv in seq_along(ap$vpath)) {

      if (length( ap$vpath[[jv]])<=2) {
        next
      }

      j.dir <- get_path_direction(i.cn.seed.ig,
                                  ap$vpath[[jv]],
                                  ap$epath[[jv]])

      j.cd <- ap$epath[[jv]]$chemform_diff
      j.cd <- MSCC::chemform_multi(j.cd,j.dir,return = "chemform")

      j.cd.temp <- sapply(seq_along(j.cd),function(x){
        MSCC:::chemform_sum(j.cd[1:x])
      })
      j.cd.temp <- chemform_simplify(j.cd.temp)
      j.cd.temp <- chemform_remove_iso(j.cd.temp)
      j.cd.exist <- (j.cd.temp %in% ad.mass.diff$chemform_diff)



      if (all(j.cd.exist)) {

        message(iv," ",jv," ",sum(j.cd.exist))
        #message(paste0(i.ep$eid,collapse = ";"))
        #edge.remain[[i.ring]] <- i.ep$eid

      }

    }



    #if (i.cd.cum!="") break

  }

}



# Thu Nov 20 14:19:14 2025 ------------------------------
{

  ap <- igraph::all_shortest_paths(i.cn.seed.ig,
                       from = 1,
                       mode = "all")

  for(i.ring in seq_len(length(ap$vpath)) ){

    if (length( ap$vpath[[i.ring]])<=1) {
      next
    }
    i.dir <- get_path_direction(i.cn.seed.ig,
                                ap$vpath[[i.ring]],
                                ap$epath[[i.ring]])
    i.cd <- ap$epath[[i.ring]]$chemform_diff
    i.cd <- MSCC::chemform_multi(i.cd,i.dir,return = "chemform")
    i.cd.cum <- MSCC:::chemform_sum(i.cd)
    message(i.cd.cum)
    sapply(seq_along(i.cd),function(x){
      MSCC:::chemform_sum(i.cd[1:x])
    })
    #if (i.cd.cum!="") break

  }

}

# Thu Nov 20 15:12:09 2025 ------------------------------
{

  for (i.cn.seed in i.cn.seed.vda$name) {
    all.from <- c( incident(i.cn.seed.ig,i.cn.seed,mode = "out")$adduct.from,
                   incident(i.cn.seed.ig,i.cn.seed,mode = "in")$adduct.to)
    print(all.from)
  }




}

# Thu Nov 20 13:04:18 2025 ------------------------------
{

  ### find all ring
  ### chemform cum in ring
  ### chenform cum exist


  i.cn.seed.ig.ring <- igraph::simple_cycles(
    i.cn.seed.ig,mode = "all")

  ring.node.forms <- list()
  for(i.ring in seq_len(length(i.cn.seed.ig.ring$vertices)) ){

    i.ep <-  i.cn.seed.ig.ring$edges[[i.ring]]
    i.dir <- get_path_direction(i.cn.seed.ig,
                                i.cn.seed.ig.ring$vertices[[i.ring]],
                                i.ep)
    i.cd <- i.cn.seed.ig.ring$edges[[i.ring]]$chemform_diff
    i.cd <- MSCC::chemform_multi(i.cd,i.dir,return = "chemform")
    #i.cd.cum <- MSCC:::chemform_sum(i.cd)
    i.cd.temp <- sapply(seq_along(i.cd),function(x){
      MSCC:::chemform_sum(i.cd[1:x])
    })
    i.cd.temp <- chemform_simplify(i.cd.temp)
    i.cd.temp <- chemform_remove_iso(i.cd.temp)
    i.cd.exist <- (i.cd.temp%in% ad.mass.diff$chemform_diff)

    if (all(i.cd.exist)) {

      i.ring.ig  <- igraph_filter_edge(i.cn.seed.ig,
                                       which(E(i.cn.seed.ig)$eid %in%  i.ep$eid))
      i.ring.vform <- get_pave_ring_vertex_form(i.ring.ig)
      if (any(lengths(i.ring.vform) > 1)) next
      i.ring.vform <- unlist(i.ring.vform)
      #print(i.ring.vform)
      ring.node.forms[[i.ring]] <- i.ring.vform
    }
  }
  ring.node.form <- do.call(bind_rows, ring.node.forms)
  ring.node.form.group(ring.node.form)

  i.cn.seed.ig.loop <- igraph_filter_edge(i.cn.seed.ig,
                                          which(E(i.cn.seed.ig)$eid %in%
                                                  unique(unlist(edge.remain))))
  i.cn.seed.ig.exclude.ring <- igraph_remove_vertex(
    i.cn.seed.ig, names(V(i.cn.seed.ig.loop)))


}

# Thu Nov 27 11:33:51 2025 ------------------------------
p1 <- plot_xcms_TIC(res$OE480_120k_ppm10_sn10.rdata@xcmsData$NegativeMS1,title = "480_120k")
p2 <- plot_xcms_TIC(res$OE480_480k_ppm10_sn10.rdata@xcmsData$NegativeMS1,title = "480_480k")
p3 <- plot_xcms_TIC(res$QEplus_ppm10_sn10.rdata@xcmsData$NegativeMS1,title = "QE_plus")
p4 <- plot_xcms_TIC(res$Astral_ppm10_sn10.rdata@xcmsData$NegativeMS1,title = "Astral")

p <- p1/p2/p3/p4


p1 <- plot_xcms_TIC(res$OE480_120k_ppm10_sn10.rdata@xcmsData$PositiveMS1,title = "480_120k")
p2 <- plot_xcms_TIC(res$OE480_480k_ppm10_sn10.rdata@xcmsData$PositiveMS1,title = "480_480k")
p3 <- plot_xcms_TIC(res$QEplus_ppm10_sn10.rdata@xcmsData$PositiveMS1,title = "QE_plus")
p4 <- plot_xcms_TIC(res$Astral_ppm10_sn10.rdata@xcmsData$PositiveMS1,title = "Astral")

p <- p1/p2/p3/p4
open_plot_win(p,10,10)



# Thu Nov 27 13:46:04 2025 ------------------------------
{
  library(tidyverse)
  {
    hs <- read.csv("d:/Share/Hang Seng TECH Historical Data.csv")%>%
      dplyr::mutate(Date = as.Date(Date, format = "%m/%d/%Y")+1,
                    Price = parse_number(Price),
                    diff = c(0,diff(Price)),
                    change = as.numeric(sub("%", "", Change..)) ,
                    type = "hs")

    ns <- read.csv("d:/Share/NASDAQ Composite Historical Data.csv")%>%
      dplyr::mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
                    Price = parse_number(Price),
                    diff = c(0,diff(Price)),
                    change = as.numeric(sub("%", "", Change..)) ,
                    type = "ns")

    overlap.date <- intersect(hs$Date,ns$Date)%>%sort()
    hs <- hs[match(overlap.date,hs$Date),]
    ns <- ns[match(overlap.date,ns$Date),]

    plot.data <- rbind(hs,ns)


    cor.test(hs$diff,ns$diff)
    #plot(hs$change,ns$change)


    plot.data <-
      data.frame(date = overlap.date,hs = hs$change,ns = ns$change)


    p1 <-ggplot(plot.data,aes(x = ns , y = hs))+
      geom_point()+
      geom_smooth()+
      ggpubr::stat_cor()+
      labs(x = "NASDAQ (day 0)", y = "HS TECH (day -1)")+
      theme_bw()


  }
  {
    hs <- read.csv("d:/Share/Hang Seng TECH Historical Data.csv")%>%
      dplyr::mutate(Date = as.Date(Date, format = "%m/%d/%Y")+0,
                    Price = parse_number(Price),
                    diff = c(0,diff(Price)),
                    change = as.numeric(sub("%", "", Change..)) ,
                    type = "hs")

    ns <- read.csv("d:/Share/NASDAQ Composite Historical Data.csv")%>%
      dplyr::mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
                    Price = parse_number(Price),
                    diff = c(0,diff(Price)),
                    change = as.numeric(sub("%", "", Change..)) ,
                    type = "ns")

    overlap.date <- intersect(hs$Date,ns$Date)%>%sort()
    hs <- hs[match(overlap.date,hs$Date),]
    ns <- ns[match(overlap.date,ns$Date),]

    plot.data <- rbind(hs,ns)


    cor.test(hs$diff,ns$diff)
    #plot(hs$change,ns$change)


    plot.data <-
      data.frame(date = overlap.date,hs = hs$change,ns = ns$change)


    p2 <-ggplot(plot.data,aes(x = ns , y = hs))+
      geom_point()+
      geom_smooth()+
      ggpubr::stat_cor()+
      labs(x = "NASDAQ (day 0)", y = "HS TECH (day 0)")+
      theme_bw()


  }
  {
    hs <- read.csv("d:/Share/Hang Seng TECH Historical Data.csv")%>%
      dplyr::mutate(Date = as.Date(Date, format = "%m/%d/%Y")-1,
                    Price = parse_number(Price),
                    High = parse_number(High),
                    Low = parse_number(Low),
                    diff = c(0,diff(Price)),
                    change = as.numeric(sub("%", "", Change..)) ,
                    hr = (1-Price/(1+change/100)/High)*100,
                    lr = (1-Price/(1+change/100)/Low)*100,
                    type = "hs")%>%
      dplyr::arrange(Date)

    ns <- read.csv("d:/Share/NASDAQ Composite Historical Data.csv")%>%
      dplyr::mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
                    Price = parse_number(Price),
                    High = parse_number(High),
                    Low = parse_number(Low),
                    diff = c(0,diff(Price)),
                    change = as.numeric(sub("%", "", Change..)) ,
                    hr = (1-Price/(1+change/100)/High)*100,
                    lr = (1-Price/(1+change/100)/Low)*100,
                    type = "hs")%>%
      dplyr::arrange(Date)

    overlap.date <- intersect(hs$Date,ns$Date)%>%sort()
    hs <- hs[match(overlap.date,hs$Date),]
    ns <- ns[match(overlap.date,ns$Date),]

    plot.data <- rbind(hs,ns)


    cor.test(hs$diff,ns$diff)
    #plot(hs$change,ns$change)


    plot.data <-
      data.frame(date = overlap.date,hs = hs$change,ns = ns$change)


    p3 <-ggplot(plot.data,aes(x = ns , y = hs))+
      geom_point()+
      geom_smooth()+
      ggpubr::stat_cor()+
      labs(x = "NASDAQ (day 0)", y = "HS TECH (day 1)")+
      theme_bw()
    p3

    diff.thresh <- 2
    start.point <- 0
    strategy <- data.frame(date = overlap.date,
               hs = hs$change,
               ns = ns$change)%>%
      dplyr::mutate(
        x = case_when(
          ns > start.point ~ (hs > ns - diff.thresh),
          ns < -start.point ~ (hs < ns + diff.thresh),
          T ~ T
        ),
        y = hs > ns - diff.thresh,
        z = hs < ns + diff.thresh,
        s = case_when(
          ns>0&hs>0~"11",
          ns>0&hs<0~"10",
          ns<0&hs>0~"01",
          ns<0&hs<0~"00"
        )
      )

    x.exist = case_when(
      ns$change > start.point ~ (hs$lr  < ns$change - diff.thresh),
      ns$change < -start.point ~ (hs$hr > ns$change + diff.thresh),
      T~F
    )
    sum(strategy$x)/length(overlap.date)
    sum(x.exist)/length(overlap.date)


    sum(strategy$y)/length(overlap.date)
    sum(strategy$z)/length(overlap.date)
    table(strategy$s)

  }
  {
    hs <- read.csv("d:/Share/Hang Seng TECH Historical Data.csv")%>%
      dplyr::mutate(Date = as.Date(Date, format = "%m/%d/%Y")-2,
                    Price = parse_number(Price),
                    diff = c(0,diff(Price)),
                    change = as.numeric(sub("%", "", Change..)) ,
                    type = "hs")

    ns <- read.csv("d:/Share/NASDAQ Composite Historical Data.csv")%>%
      dplyr::mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
                    Price = parse_number(Price),
                    diff = c(0,diff(Price)),
                    change = as.numeric(sub("%", "", Change..)) ,
                    type = "ns")

    overlap.date <- intersect(hs$Date,ns$Date)%>%sort()
    hs <- hs[match(overlap.date,hs$Date),]
    ns <- ns[match(overlap.date,ns$Date),]

    plot.data <- rbind(hs,ns)


    cor.test(hs$diff,ns$diff)
    #plot(hs$change,ns$change)


    plot.data <-
      data.frame(date = overlap.date,hs = hs$change,ns = ns$change)


    p4 <-ggplot(plot.data,aes(x = ns , y = hs))+
      geom_point()+
      geom_smooth()+
      ggpubr::stat_cor()+
      labs(x = "NASDAQ (day 0)", y = "HS TECH (day 2)")+
      theme_bw()


  }

  p <- p1+p2+p3 + p4+plot_layout(nrow = 1)
  open_plot_win(p, 16,4)


  {
    hs <- read.csv("d:/Share/Hang Seng TECH Historical Data.csv")%>%
      dplyr::mutate(Date = as.Date(Date, format = "%m/%d/%Y")-1,
                    Price = parse_number(Price),
                    High = parse_number(High),
                    Low = parse_number(Low),
                    diff = c(0,diff(Price)),
                    change = as.numeric(sub("%", "", Change..)) ,
                    hr = (1-Price/(1+change/100)/High)*100,
                    lr = (1-Price/(1+change/100)/Low)*100,
                    type = "hs")%>%
      dplyr::arrange(Date)

    ns <- read.csv("d:/Share/NASDAQ Composite Historical Data.csv")%>%
      dplyr::mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
                    Price = parse_number(Price),
                    High = parse_number(High),
                    Low = parse_number(Low),
                    diff = c(0,diff(Price)),
                    change = as.numeric(sub("%", "", Change..)) ,
                    hr = (1-Price/(1+change/100)/High)*100,
                    lr = (1-Price/(1+change/100)/Low)*100,
                    type = "hs")%>%
      dplyr::arrange(Date)

    overlap.date <- intersect(hs$Date,ns$Date)%>%sort()
    hs <- hs[match(overlap.date,hs$Date),]
    ns <- ns[match(overlap.date,ns$Date),]

    plot.data <- rbind(hs,ns)


    cor.test(hs$diff,ns$diff)
    #plot(hs$change,ns$change)


    plot.data <-
      data.frame(date = overlap.date,hs = hs$change,ns = ns$change)


    p3 <-ggplot(plot.data,aes(x = ns , y = hs))+
      geom_point()+
      geom_smooth()+
      ggpubr::stat_cor()+
      labs(x = "NASDAQ (day 0)", y = "HS TECH (day 1)")+
      theme_bw()
    p3

    diff.thresh <- 2
    start.point <- 0
    strategy <- data.frame(date = overlap.date,
                           hs = hs$change,
                           ns = ns$change)%>%
      dplyr::mutate(
        x = case_when(
          ns > start.point ~ (hs > ns - diff.thresh),
          ns < -start.point ~ (hs < ns + diff.thresh),
          T ~ T
        ),
        y = hs > ns - diff.thresh,
        z = hs < ns + diff.thresh,
        s = case_when(
          ns>0&hs>0~"11",
          ns>0&hs<0~"10",
          ns<0&hs>0~"01",
          ns<0&hs<0~"00"
        )
      )

    x.exist = case_when(
      ns$change > start.point ~ (hs$lr  < ns$change - diff.thresh),
      ns$change < -start.point ~ (hs$hr > ns$change + diff.thresh),
      T~F
    )
    sum(strategy$x)/length(overlap.date)
    sum(x.exist)/length(overlap.date)


    sum(strategy$y)/length(overlap.date)
    sum(strategy$z)/length(overlap.date)
    table(strategy$s)

    start.point <- 0
    diff.thresh <- 2.9
    earn.point <- 1.3
    earn <- rep(NA,nrow(hs))

    for (i in 1:nrow(hs)) {
      i.nc <- ns$change[i]
      if (abs(i.nc) > start.point) {

        if ( hs$lr[i] < i.nc - diff.thresh ) {
          i.b <- mean(   (i.nc - diff.thresh)   ,
                         hs$lr[i]    )
          i.b <- (i.nc - diff.thresh)
          if (hs$hr[i] > i.nc - diff.thresh + earn.point ) {
            i.s <- ( i.nc - diff.thresh + earn.point )
          }else{
            i.s <- hs$change[i]
          }

          i.earn <- (i.s - i.b)
          earn[[i]] <- i.earn
          #message(i,";",i.b,";",i.s)
        }
      }
    }

    hist(earn)
    sum(earn ,na.rm = T)
    sum(!is.na(earn) ,na.rm = T)


    sum(between.range(hs$change,cbind(ns$change-2.9,0+2.9)))/938

  }

}

# Thu Nov 27 16:25:58 2025 ------------------------------
{
  library(ggplot2)
  library(ggpp)

  # Example data
  df <- data.frame(
    x = rnorm(50),
    y = rnorm(50)
  )

  # A small table to display
  tbl <- data.frame(
    Stat = c("Mean X", "Mean Y"),
    Value = c(mean(df$x), mean(df$y))
  )

  p <- ggplot() +
    #geom_point() +
    geom_table_npc(
      data = data.frame(x = NA,
                        y = NA,
                        tbl = I(list(tbl))),
      aes(npcx = 0.5, npcy = 0.5, label = tbl)
    ) +
    theme_minimal()

  open_plot_win(p,10,10)

}

# Tue Dec  2 10:24:24 2025 ------------------------------
{

  g.mmu <- buildGraphFromKEGGREST(
    organism = "hsa")

  buildDataFromGraph(
    keggdata.graph = g.mmu,
    internalDir = TRUE,
    matrices = c("hypergeom", "diffusion", "pagerank"),
    normality = c("diffusion", "pagerank") )

}

# Thu Dec  4 14:57:13 2025 ------------------------------
{
  total_iterations <- 100

  # 1. Initialize the progress bar object
  # :bar displays the bar, :percent shows progress, :elapsed shows time spent
  pbar <- progress_bar$new(
    format = "[:bar] :percent in :elapsed, current index: :current",
    total = 100,
    clear = FALSE
  )

  # 2. Start the for loop
  for (i in 1:total_iterations) {

    # Simulate work
    Sys.sleep(0.05)

    # 3. Increment the progress bar
    pbar$tick()
  }


  # 1. Define the parameters
  total_iterations <- 50
  pb <- txtProgressBar(min = 0, max = total_iterations, style = 5) # Style 3 is the most common text format

  # 2. Start the for loop
  for (i in 1:total_iterations) {

    # --- START OF YOUR LOOP BODY ---

    # Simulate some work
    Sys.sleep(0.1)

    # You can store results here if needed

    # --- END OF YOUR LOOP BODY ---

    # 3. Update the progress bar
    setTxtProgressBar(pb, i)
  }

  # 4. Close the progress bar when the loop is done
  close(pb)


}

# Mon Dec  8 19:30:24 2025 ------------------------------
{

  dur <- 60*6000
  act.time <- 25
  bh.cd <- 140
  gt.cd <- 280
  gb.cd <- 105
  bh.time <- seq(0,dur,bh.cd)
  gt.time <- seq(0,dur,gt.cd)
  gb.time <- seq(0,dur,gb.cd)

  all.tp <- unique(c(bh.time,gt.time,gb.time))
  all.times <- data.table(
    time = all.tp
  )%>%
    dplyr::mutate(
      bh =   time %in% bh.time,
      gt =   time %in% gt.time,
      gb =   time %in% gb.time,

      bh = ifelse(bh,2.5,1),
      gt = ifelse(gt,8 * 1.8,1),
      gb = ifelse(gb,2.2,1),
      tt = gt*bh*gb
      )%>%
    dplyr::arrange(time)

  print(
    (sum(all.times$tt * act.time) + dur - nrow(all.times) * act.time ) / dur
  )


}

# Wed Dec 10 19:39:16 2025 ------------------------------
{

  x.mz <- normalize_max_min(select.fdf$mzmed)*0.6+0.2
  mz.label.size = 6
  sp_anno_fun_list <- list(
    function(x, y, w, h) {
      grid.segments(x0 = unit(x.mz[1], "npc"), y0 = unit(0, "npc"),
                    x1 = unit(x.mz[1], "npc"), y1 = unit(0.8, "npc"))

      grid.text(label = str_digit(select.fdf$mzmed[1],4),
                gp = gpar(fontsize = mz.label.size),
                x = unit(x.mz[1], "npc"),  y = unit(0.9, "npc"))

      grid.segments(x0 = unit(0, "npc"), y0 = unit(0, "npc"),
                    x1 = unit(1, "npc"), y1 = unit(0, "npc"))
    },
    function(x, y, w, h) {
      grid.segments(x0 = unit(x.mz[2], "npc"), y0 = unit(0, "npc"),
                    x1 = unit(x.mz[2], "npc"), y1 = unit(0.8, "npc"))
      grid.text(label = str_digit(select.fdf$mzmed[2],4),
                gp = gpar(fontsize = mz.label.size),
                x = unit(x.mz[2], "npc"),  y = unit(0.9, "npc"))

      grid.segments(x0 = unit(0, "npc"), y0 = unit(0, "npc"),
                    x1 = unit(1, "npc"), y1 = unit(0, "npc"))
    },
    function(x, y, w, h) {
      grid.segments(x0 = unit(x.mz[3], "npc"), y0 = unit(0, "npc"),
                    x1 = unit(x.mz[3], "npc"), y1 = unit(0.8, "npc"))
      grid.text(label = str_digit(select.fdf$mzmed[3],4),
                gp = gpar(fontsize = mz.label.size),
                x = unit(x.mz[3], "npc"),  y = unit(0.9, "npc"))

      grid.segments(x0 = unit(0, "npc"), y0 = unit(0, "npc"),
                    x1 = unit(1, "npc"), y1 = unit(0, "npc"))
    },
    function(x, y, w, h) {
      grid.segments(x0 = unit(x.mz[4], "npc"), y0 = unit(0, "npc"),
                    x1 = unit(x.mz[4], "npc"), y1 = unit(0.8, "npc"))
      grid.text(label = str_digit(select.fdf$mzmed[4],4),
                gp = gpar(fontsize = mz.label.size),
                x = unit(x.mz[4], "npc"),  y = unit(0.9, "npc"))

      grid.segments(x0 = unit(0, "npc"), y0 = unit(0, "npc"),
                    x1 = unit(1, "npc"), y1 = unit(0, "npc"))
    }
  )
  names(sp_anno_fun_list) <- select.node$pave_pattern

  hm <- Heatmap(pave.mat,
                show_heatmap_legend = F,
                col = colramp(),
                cluster_columns = F,
                cluster_rows = F,,
                row_names_side = "left",
                rect_gp = gpar(col = "black"),
                right_annotation = HeatmapAnnotation(
                  mz = anno_customize(
                    x= select.node$pave_pattern,
                    graphics = sp_anno_fun_list
                    ),
                  which = "row",
                  border  = T,
                  annotation_name_rot = 0,
                  width =  unit(2,"inch")
                ),
                column_names_rot = -45)
  hm

}

# Sat Dec 13 22:21:22 2025 ------------------------------
{
  obj <- MSdev("d:/data/2025.12.13.lie/rawdata/")
  obj <- MSdev_msConvert(obj)
  #obj <- MSdev_checkSampleInfo(obj)
  sam.info <- obj@sampleInfo%>%
    dplyr::filter(file.exists(msData.files))

  ef <- c()
  for ( i.f in sam.info$msData.files) {
    try.res <- ""
    try.res <- try(xcms.xcms <- readSRMData(i.f))
    if (any(grepl("Error ",try.res))) {
      ef <- c(ef,i.f)
    }
  }

  xcms.xcms <- readSRMData(sam.info$msData.files)
  xcms.pks <- findChromPeaks(xcms.xcms,param = CentWaveParam())
  xcms.pks <- groupChromPeaks(xcms.pks,
                              param = PeakDensityParam(sampleGroups = sam.info$group))


}

# Sat Dec 13 23:14:11 2025 ------------------------------
{
  sample.info <- openxlsx::read.xlsx("d:/data/2025.12.13.lie/Lipid_SampleInfo.xlsx")
  data <-  openxlsx::read.xlsx("d:/data/2025.12.13.lie/Lipid_SampleInfo.xlsx",sheet = 2)
  sample.info <- sample.info%>%
    dplyr::filter(Sample.Name %in% data$Sample.Name)%>%
    dplyr::distinct(Sample.Name,.keep_all = T)
  data.matrix <-data %>%
    dplyr::filter(Sample.Name %in% sample.info$Sample.Name)%>%
    dplyr::distinct(Sample.Name,.keep_all = T)%>%
    dplyr::select(!contains("IS"))%>%
    column_to_rownames("Sample.Name")
  data.matrix <- data.matrix[sample.info$Sample.Name,]

  p.before <- plot_PCA(data.matrix,pca.group = sample.info$Sample.Type,
                       col = c("#EB5640","#3B5B8B"))+
    labs(title = "PCA before normalization")

  vsn_fit <- vsn2(as.matrix(data.matrix))
  mat_vsn <- predict(vsn_fit, as.matrix(data.matrix))
  p.vsn <- plot_PCA(mat_vsn,pca.group = sample.info$Sample.Type,
                    col = c("#EB5640","#3B5B8B"))+
    labs(title = "PCA after normalization")

  open_plot_win(p.before+p.vsn,6,3)
  export_graph2pdf(p.before+p.vsn,file_path = "d:/temp/figures.pdf",
                   width = 6,height = 3)


  ### rt
  {

    sample.info <- openxlsx::read.xlsx("d:/data/2025.12.13.lie/Lipid_SampleInfo.xlsx")
    data <-  openxlsx::read.xlsx("d:/data/2025.12.13.lie/Lipid_SampleInfo.xlsx",sheet = 3)
    sample.info <- sample.info%>%
      dplyr::filter(Sample.Name %in% data$Sample.Name)%>%
      dplyr::distinct(Sample.Name,.keep_all = T)%>%
      dplyr::arrange(Inject.Order)
    plot.data <-data %>%
      dplyr::filter(Sample.Name %in% sample.info$Sample.Name)%>%
      dplyr::distinct(Sample.Name,.keep_all = T)%>%
      dplyr::select(!contains("IS"),Sample.Name)%>%
      column_to_rownames("Sample.Name")%>%
      pivot_longer(everything(),names_to = "name",values_to = "value")%>%
      dplyr::group_by(name)%>%
      dplyr::mutate(value = value - mean(value,na.rm = T))


    p1 <- ggplot(plot.data)+
      geom_boxplot(aes(x = name,y = value ),outliers = F)+
      geom_hline(yintercept = c(-0.5,0.5) ,linetype = "dashed" , col = "grey" )+
      geom_jitter(aes(x = name,y = value ),pch = 21,
                  col = "black",fill = "red",alpha = 0.5)+
      labs(x = NULL, y="Retention Time error")+
      ylim(c(-1,1))+
      theme_classic()
    p1

    open_plot_win(p1,5,3)







  }


  ### intensity
  {

    sample.info <- openxlsx::read.xlsx("d:/data/2025.12.13.lie/Lipid_SampleInfo.xlsx")
    data <-  openxlsx::read.xlsx("d:/data/2025.12.13.lie/Lipid_SampleInfo.xlsx",sheet = 2)
    sample.info <- sample.info%>%
      dplyr::filter(Sample.Name %in% data$Sample.Name)%>%
      dplyr::distinct(Sample.Name,.keep_all = T)%>%
      dplyr::arrange(Inject.Order)
    plot.data <-data %>%
      dplyr::filter(Sample.Name %in% sample.info$Sample.Name)%>%
      dplyr::distinct(Sample.Name,.keep_all = T)%>%
      dplyr::select(!contains("IS"),Sample.Name)%>%
      column_to_rownames("Sample.Name")%>%
      pivot_longer(everything(),names_to = "name",values_to = "value")%>%
      dplyr::group_by(name)%>%
      dplyr::mutate(value = value/median(value,na.rm = T)*100)


    p2 <- ggplot(plot.data)+
      geom_boxplot(aes(x = name,y = value ),outliers = F)+
      geom_hline(yintercept = c(70,130) ,linetype = "dashed" , col = "grey" )+
      geom_jitter(aes(x = name,y = value ),pch = 21,
                  col = "black",fill = "red",alpha = 0.5)+
      labs(x = NULL, y="Variation in Relative Area (%)")+
      ylim(c(0,300))+
      theme_classic()
    p2

    #open_plot_win(p2,5,3)

    #open_plot_win(p1+p2,6,4)
    export_graph2pdf(p1+p2,file_path = "d:/temp/figures.pdf",
                     width = 6,height = 4,append = T)




  }



  ### int-injec,
  {


    data.matrix <-data %>%
      dplyr::filter(Sample.Name %in% sample.info$Sample.Name)%>%
      dplyr::distinct(Sample.Name,.keep_all = T)%>%
      column_to_rownames("Sample.Name")
    data.matrix <- data.matrix[sample.info$Sample.Name,]
    vsn_fit <- vsn2(as.matrix(data.matrix))
    mat_vsn <- predict(vsn_fit, as.matrix(data.matrix))%>%as.data.frame()


    plot.data <- sample.info %>%
      dplyr::mutate(value = data.matrix[Sample.Name,"IS_DSPCD70"])%>%
      dplyr::distinct(Sample.Name,.keep_all = T)%>%
      dplyr::mutate(
        injection.order = Inject.Order,
        group = Sample.Type,
        injection.order = as.numeric(injection.order)
      )

    x <- plot.data$value
    p1 <- ggplot(plot.data)+
      geom_point(aes(x = injection.order, y = value,fill= group),
                 pch = 21)+
      geom_hline(
        yintercept = c(mean(x) + sd(x) * c(-1,1)))+
      geom_hline(
        yintercept = c(mean(x) + sd(x) * c(-0.3,0.3)), linetype = "dashed")+
      scale_fill_manual(values = c("#EB5640","#3B5B8B"))+
      ylim(c(0,quantile(plot.data$value,0.98)))+
      labs(title = "Raw Peak",x = "Injection Order", y="Intensity",fill = "")+
      theme_classic()+
      theme(legend.position = "top",plot.title = element_text(hjust = 0.5))

    p1
    #open_plot_win(p,5,3)


    plot.data <- sample.info %>%
      dplyr::mutate(value = mat_vsn[Sample.Name,"IS_PED31"],
                    value=2^value)%>%
      dplyr::distinct(Sample.Name,.keep_all = T)%>%
      dplyr::mutate(
        injection.order = Inject.Order,
        group = Sample.Type,
        injection.order = as.numeric(injection.order)
      )

    x <- plot.data$value
    p2 <- ggplot(plot.data)+
      geom_point(aes(x = injection.order, y = value,fill= group),
                 pch = 21)+
      geom_hline(
        yintercept = c(mean(x) + sd(x) * c(-1,1)))+
      geom_hline(
        yintercept = c(mean(x) + sd(x) * c(-0.3,0.3)), linetype = "dashed")+
      ylim(c(0,quantile(plot.data$value,0.98)))+
      scale_fill_manual(values = c("#EB5640","#3B5B8B"))+
      labs(title = "Corrected Peak",x = "Injection Order", y="Intensity",fill = "")+
      theme_classic()+
      theme(legend.position = "top",plot.title = element_text(hjust = 0.5))

    p2
    open_plot_win(p2,5,3)

    open_plot_win(p1/p2,6,8)
    export_graph2pdf(p1/p2,file_path = "d:/temp/figures.pdf",
                     width = 6,height = 8,append = T)
  }


}
# Mon Dec 15 15:44:16 2025 ------------------------------
{

  MSdev_annotation(object,cpdb_path = cpdb_path)




  ###temp
  {
    system.time(
      matched.df <- match_mz_rt(mz1 = xcms.featuredef$mzmed,
                                mz2 = cp.adduct$chemform.adduct.mz,
                                mz.ppm = mz.ppm)
    )

    system.time(
      matched.df2 <- match_mz_foverlaps(mz1 = xcms.featuredef$mzmed,
                                mz2 = cp.adduct$chemform.adduct.mz,
                                ppm = mz.ppm)
    )
  }
}


# Mon Dec 15 20:46:52 2025 ------------------------------
{

  Sys.setenv(
    http_proxy  = "http://127.0.0.1:7897",
    https_proxy = "http://127.0.0.1:7897"
  )
  getSymbols("1810.HK", src = "yahoo",
             from = "2000-01-02",
             to   = TODAY())

  plot(`1810.HK`)

  plot(yearlyReturn(SPY))
  periodReturn(SPY)


}

# Mon Dec 15 21:58:24 2025 ------------------------------
{

  library(tidyquant)
  library(dplyr)
  library(TTR)

  getSymbols("SPY", from = "2005-01-01")
  prices <- tq_get("SPY", from = "2005-01-01")


  data <- prices %>%
    tq_transmute(
      select     = adjusted,
      mutate_fun = dailyReturn,
      col_rename = "ret"
    ) %>%
    left_join(prices, by = "date") %>%
    mutate(
      ma200  = SMA(adjusted, 200),
      signal = ifelse(adjusted > ma200, 1, 0),
      strategy_ret = signal * ret
    )




  library(tidyquant)
  library(dplyr)
  library(TTR)

  prices <- tq_get("SPY", from = "2005-01-01")

  data <- prices %>%
    mutate(
      ma200 = SMA(adjusted, n = 200),
      signal = ifelse(adjusted > ma200, 1, 0),
      ret = dailyReturn(adjusted),
      strategy_ret = signal * ret
    )

  # 累计收益
  cumprod(1 + na.omit(data$strategy_ret))


  library(tidyquant)
  library(dplyr)
  library(TTR)

  prices <- tq_get("SPY", from = "1998-01-01",to = "2025-12-11")

  tsmom <- prices %>%
    mutate(
      ret = log(adjusted / lag(adjusted)),
      ma200 = SMA(adjusted, 200),
      signal = ifelse(adjusted > ma200, 1, 0),
      signal = lag(signal, 1),
      strategy_ret = signal * ret
    ) %>%
    filter(!is.na(strategy_ret))

  # 累计收益
  tsmom <- tsmom %>%
    mutate(
      cum_strategy = exp(cumsum(strategy_ret)),
      cum_buyhold  = exp(cumsum(ret))
    )

  ggplot(tsmom)+
    geom_point(aes(x = date, y = cum_buyhold),col = "blue")+
    geom_point(aes(x = date, y = cum_strategy),col = "red")


  max_dd <- function(x) {
    cummax_x <- cummax(x)
    min((x - cummax_x) / cummax_x)
  }

  max_dd(tsmom$cum_strategy)
  max_dd(tsmom$cum_buyhold)

  sd(tsmom$strategy_ret, na.rm = TRUE) * sqrt(252)
  sd(tsmom$ret, na.rm = TRUE) * sqrt(252)


  mean(tsmom$strategy_ret, na.rm = TRUE) /
    sd(tsmom$strategy_ret, na.rm = TRUE) * sqrt(252)


  {

    symbols <- c("SPY", "IEF", "GLD","QQQ")
    symbols <- c("VOO","QQQ","1810.HK")
   # symbols <-"1810.HK"
    monthly_prices <- tq_get(
      symbols,
      from = "2015-01-01",
      to = TODAY(),
      periodicity = "daily"
    ) %>%
      select(symbol, date, adjusted) %>%
      rename(close = adjusted)%>%
      dplyr::filter(!is.na(close))

    tsmom <- monthly_prices %>%
      group_by(symbol) %>%
      mutate(
        ret    = log(close / lag(close)),
        ma10   = SMA(close, 30),
        signal = ifelse(close > ma10, 1, 0),
        signal = lag(signal, 1),
        strat_ret = signal * ret
      ) %>%
      filter(!is.na(strat_ret)) %>%
      ungroup()

    portfolio <- tsmom %>%
      group_by(date) %>%
      summarise(
        portfolio_ret = mean(strat_ret, na.rm = TRUE),
        buyhold_ret   = mean(ret, na.rm = TRUE)
      ) %>%
      mutate(
        cum_trend   = exp(cumsum(portfolio_ret)),
        cum_buyhold = exp(cumsum(buyhold_ret))
      )

    p <- portfolio %>%
      select(date, cum_trend, cum_buyhold) %>%
      pivot_longer(-date) %>%
      ggplot(aes(date, value, color = name)) +
      geom_line(linewidth = 1) +
      #scale_y_log10() +
      labs(
        title = "Multi-Asset Trend Following vs Buy & Hold",
        y = "Cumulative Return (log scale)",
        color = ""
      ) +
      theme_minimal()
    print(p)

    max_dd(portfolio$cum_trend)
    max_dd(portfolio$cum_buyhold)



    sd(portfolio$cum_trend, na.rm = TRUE) * sqrt(252)
    sd(portfolio$cum_buyhold, na.rm = TRUE) * sqrt(252)



    }

}
# Thu Dec 18 13:00:07 2025 ------------------------------
{


  cfmd <- get_CFM_data_from_smiles(
    smiles = "Nc1ccc(O)cc1",compound_id = "gshSSSS")

  shiny_vis_cfmd_trans(cfmd)

  {


  }









}
# Sun Dec 28 20:49:31 2025 ------------------------------
{

  obj <- MSdev("d:/data/2025.12.26.PAVE2/2021222_PAVE/data/")
  obj <- MSdev_msConvert(obj)
  obj <- MSdev_checkSampleInfo(obj)
  obj <- MSdev_set_param(
    obj,
    findChromPeaks =
      xcms::CentWaveParam(
        ppm = 10,
        prefilter = c(3,1000),
        peakwidth = c(10,30),
        snthresh = 10,
        fitgauss = T),
    groupChromPeaks =
      xcms::PeakDensityParam(
        sampleGroups = "A",
        minFraction = 0.6,
        binSize = 0.002,
        bw = 12,
        ppm = 10))
  obj <- MSdev_xcmsProcessing(obj)
  MSdev_save(obj)

}
# Tue Dec 30 10:09:25 2025 ------------------------------
{

  ppm = 10
  sn = 100
  ### 120K
  {
    obj <- MSdev("d:/temp/astral/data/")
    obj <- MSdev_msConvert(obj)
    obj <- MSdev_checkSampleInfo(obj)
    obj <- MSdev_set_param(
      obj,
      findChromPeaks =
        xcms::CentWaveParam(
          ppm = ppm,
          prefilter = c(3,1000),
          peakwidth = c(10,30),
          snthresh = sn,
          fitgauss = T),
      groupChromPeaks =
        xcms::PeakDensityParam(
          sampleGroups = "A",
          minFraction = 0.6,
          binSize = 0.002,
          bw = 12,
          ppm = ppm))
    obj <- MSdev_xcmsProcessing(obj)
    MSdev_save(obj)


    xcms <- obj@xcmsData$PositiveMS1
    pks <- table(chromPeaks(xcms)[,"sample"])
    names(pks) <- pData(xcms)[,1]
    print(pks)
    nrow(featureDefinitions(xcms))

  }


  ppm = 10
  sn = 100
  {

    OE480 <- MSdev("d:/temp/oe480/data/")
    OE480 <- MSdev_msConvert(OE480)
    OE480 <- MSdev_checkSampleInfo(OE480)
    OE480 <- MSdev_set_param(
      OE480,
      findChromPeaks =
        xcms::CentWaveParam(
          ppm = ppm,
          prefilter = c(3,1000),
          peakwidth = c(10,30),
          snthresh = sn,
          fitgauss = T),
      groupChromPeaks =
        xcms::PeakDensityParam(
          sampleGroups = "A",
          minFraction = 0.6,
          binSize = 0.002,
          bw = 12,
          ppm = ppm))
    OE480 <- MSdev_xcmsProcessing(OE480)
    xcms <- OE480@xcmsData$PositiveMS1
    pks <- table(chromPeaks(xcms)[,"sample"])
    names(pks) <- pData(xcms)[,1]
    print(pks)
    nrow(featureDefinitions(xcms))

  }


  obj <- load_demo()
  a <- MSdev_annotation(obj)

  obj <- MSdev_load("d:/data/2025.12.26.PAVE2/PAVE_With_Params/OE480_480k_ppm10_sn10.rdata")
  plot_xcms_adjustedRT(obj@xcmsData$PositiveMS1)
  plot_xcms_adjustedRT(obj@xcmsData$NegativeMS1)

}
