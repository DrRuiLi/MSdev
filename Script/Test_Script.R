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
  files <- dir("d:/data/2025.10.10.PVAE/data8600/rawdata/",full.names = T)
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
