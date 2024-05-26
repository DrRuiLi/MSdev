#' @title  expand_isoadduct_from_formula
#' @description Using formula, exact.mass and ion_mode, expand adduct by enviPat::adduct
#'
#' @param chem_formula a chemical formula, such as "C4H4Cl1N3O1".
#' some template see "MSCC::chem_formula_template"
#' @param ion_mode 1 OR 0
#'
#' @return df
#' @export
#'

expand_isoadduct_from_formula <- function(chem_formula,
                                          adduct.table = adduct.table,
                                          ion_mode = "positive") {


  #ion_mode <- "positive"
  #chem_formula <- chem_formula_template[1]
  adduct.rule <- MSCC::adduct.table %>%
    dplyr::filter(Ion_mode == ion_mode)


  x.formula <- chem_formula
  x.mass <- enviPat::check_chemform(isotopes,x.formula)$monoisotopic_mass
  x.isopat <- chemform_isotopes_pattern_enviPat(x.formula)%>%
    mutate(form = paste0(isotope_element , "M"))

  x_adduct <-data.frame()
  for (i in 1:nrow(x.isopat)) {

    m_formula <- x.isopat$formula[i]
    m_mass <- x.isopat$m.z[i]
    m_form <- x.isopat$form[i]
    m_adduct <- adduct.rule %>%
      rowwise() %>%
      mutate(
        formula =  enviPat::multiform(formula_in = m_formula , fact = Multi),
        formula = chemform_calculate_lc8(formula , Formula_diff),
        adduct = Adduct,
        isoadduct = sub(pattern = "M",replacement = m_form , x = adduct),
        charge = Charge,
        multi = Multi,
        ion_mode = Ion_mode,
        form = "adduct",
        exact.mz = m_mass * Multi / abs(Charge) + Mass
      ) %>%
      select("formula"  ,
             "adduct",
             "charge",
             "multi",
             "ion_mode",
             "form",
             "exact.mz")
    m_adduct
    x_adduct <- rbind(x_adduct , m_adduct)

  }



  x_adduct


  return(x_adduct)

}











#' @title match_isotopes_to_xcms_feature
#'
#' @param MS.network MS.network
#' @param xcms.features xcms.features
#' @param ppm.thresh ppm.thresh
#'
#' @return xcms
#' @export
#'

match_isotopes_to_xcms_feature <-
  function(isotopes.network , xcms.xcms, ppm.thresh = 20,rt.tol = 10,value="ratio_sub_nature")
  {


    xcms.features.def <- featureDefinitions(xcms.xcms)%>%as.data.frame()
    xcms.features.val <- featureValues(xcms.xcms , missing = "rowmin_half")
    xcms.features.intb.mean <- apply(xcms.features.val,
                                     2,mean)
    xcms.features.intb <-xcms.features.val[,which.max(xcms.features.intb.mean)]
    match.isotope <-function(x){


      isotopes.table.matched <- match_isotopes_to_featuredef(isotopes.network = x,
                                                             featuredef = xcms.features.def)

      isotopes.calced <- match_isotopes_to_featureval(isotopes.table.matched,
                                                      xcms.features.val,
                                                      value = value)





      return(isotopes.calced)
    }

    MS.network <- lapply(isotopes.network ,match.isotope )

    return(MS.network)

  }


match_isotopes_to_featuredef <- function(isotopes.network,
                                         featuredef,
                                         ppm.thresh = 20,
                                         rt.tol = 10){


  isotope.candidate <- isotopes.network
  isotope.mz <- isotope.candidate$m.z
  feature.mz <- featuredef$mzmed
  feature.id <-rownames(featuredef)
  isotope.matrix <- matrix(rep(isotope.mz,length(feature.mz)) , nrow = length(isotope.mz))
  feature.matrix <- matrix(rep(feature.mz,length(isotope.mz)) ,
                           nrow = length(isotope.mz),
                           byrow = T)
  sub.matrix <- isotope.matrix - feature.matrix
  tol.matrix <- feature.matrix * ppm.thresh*1e-6
  pass.matrix <- abs(sub.matrix) < tol.matrix

  matched.id <- which(pass.matrix,arr.ind = T)
  matched.id

  if (nrow(matched.id)==0) {
    return(NULL)
  }
  isotopes.matched <- data.frame( isotope.candidate[matched.id[,1],],
                                  feature.id = feature.id[matched.id[,2]],
                                  feature.mz = featuredef[matched.id[,2],"mzmed"],
                                  feature.rt = featuredef[matched.id[,2],"rtmed"])
  if (nrow(isotopes.matched) >1) {
    isotopes.matched <- isotopes.matched%>%
      dplyr::mutate(rt.cluster =   cutree(hclust(dist( feature.rt )),h = rt.tol))

  }else{
    isotopes.matched <- isotopes.matched%>%
      dplyr::mutate(rt.cluster =   1)
  }

  isotopes.matched <- isotopes.matched%>%
    dplyr::mutate(feature.mz.error = abs(feature.mz-m.z)/m.z*1e6)%>%
    dplyr::group_by(rt.cluster)%>%
    dplyr::filter(any( NA%in% isotope_element ))%>%
    dplyr::arrange(rt.cluster,-abundance)%>%
    dplyr::ungroup()%>%
    dplyr::group_by(rt.cluster,isotope_element)%>%
    dplyr::slice_min(feature.mz.error)%>%
    dplyr::ungroup()


  if (nrow(isotopes.matched)==0) {
    return(NULL)
  }

  isotopes.table.matched <- isotope.candidate %>%
    MSdev:::add_multi_column(unique(isotopes.matched$rt.cluster))%>%
    tidyr::pivot_longer(as.character(unique(isotopes.matched$rt.cluster)),names_to = "rt.cluster")%>%
    dplyr::select(-value)%>%
    dplyr::mutate(rt.cluster = as.numeric(rt.cluster))%>%
    dplyr::bind_rows(isotopes.matched)%>%
    dplyr::group_by(rt.cluster,formula)%>%
    dplyr::arrange(rt.cluster,feature.mz)%>%
    dplyr::slice_head(n=1)%>%
    dplyr::ungroup()%>%
    dplyr::arrange(rt.cluster,-abundance)


  return(isotopes.table.matched)









}


match_isotopes_to_featureval <- function(isotopes.table.matched,
                                         featureval,
                                         value = "ratio_sub_nature"){

  if (!value%in% c("intensity","ratio_to_base","ratio_sub_nature")) {
    stop("value should be one of \"intensity\",\"ratio_to_base\",\"ratio_sub_nature\"")
  }

  isotopes.rtcluster <-isotopes.table.matched
  isotopes.calced <- list()

  for (i in unique(isotopes.rtcluster$rt.cluster)) {

    isotopes.this.rtcluster <- isotopes.rtcluster%>%
      dplyr::filter(rt.cluster == i)%>%
      dplyr::group_by(feature.id)%>%
      dplyr::mutate(theory.abundance = sum(abundance) )%>%
      #dplyr::distinct(feature.id,.keep_all = T)%>%
      dplyr::ungroup()%>%
      dplyr::mutate(theory.abundance = case_when(is.na(feature.id)~abundance,
                                                 T~ theory.abundance))

    feature.base <- isotopes.this.rtcluster %>%
      dplyr::filter(is.na(isotope_element))%>%
      dplyr::pull(feature.id)

    matrix.int <- featureval[match(isotopes.this.rtcluster$feature.id,
                                          rownames(featureval)),colnames(featureval),drop=F]

    matrix.int.normal <- t(t(matrix.int)/(matrix.int[feature.base,]))*100

    matrix.theory <- rep(isotopes.this.rtcluster$theory.abundance,ncol(matrix.int))%>%
      matrix(ncol = ncol(matrix.int))%>%
      `rownames<-`(isotopes.this.rtcluster$feature.id)

    matrix.diff <- matrix.int.normal- matrix.theory

    matrix.output <- switch(value,
                            "intensity" = matrix.int,
                            "ratio_to_base" = matrix.int.normal,
                            "ratio_sub_nature"=matrix.diff
                            )
    colnames(matrix.output) <- colnames(matrix.output)

    isotopes.this.rtcluster.calced <- cbind(isotopes.this.rtcluster,matrix.output)%>%
      #dplyr::mutate(relative.intensity = apply(matrix.int.normal,1,mean))%>%
      dplyr::mutate(feature.isotopes = feature.base,.after = feature.id)%>%
      dplyr::arrange(-theory.abundance)%>%
      tibble::remove_rownames()


    isotopes.calced[[i]] <- isotopes.this.rtcluster.calced
  }


  isotopes.calced <- do.call("rbind",isotopes.calced)

  return(isotopes.calced)



}



