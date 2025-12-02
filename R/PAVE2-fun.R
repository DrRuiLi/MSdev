get_pave_ig_vertex_form <- function(ring.ig){


  v.form <- lapply(names(V(ring.ig)),function(x){

    #message(x)
    e <- igraph::incident(ring.ig,x,mode = "all")
    e.ends <- igraph::ends(ring.ig,e)
    ef <- e[e.ends[,1] == x]
    et <- e[e.ends[,2] == x]

    e.iso <- et$element%>%na.omit()%>%unique()
    if (isEmpty(e.iso)) e.iso <- "Natural"
    e.adduct <- c(ef$adduct.from,et$adduct.to)%>%na.omit()%>%unique()
    e.frag <- ef$fragment%>%na.omit()
    if (length(e.frag) > 0) e.frag <- MSCC::chemform_sum(e.frag)
    paste0(e.iso,";",e.adduct,";",e.frag)

  })
  names(v.form) <- names(V(ring.ig))
  return(v.form)
}

ring.node.form.group <- function(x){

  if (nrow(x)==1) return(setNames(1,1))
  xcomb <- combn(1:nrow(x),2)
  xcomb.eq <- rep(F,ncol(xcomb))
  for (i in 1:ncol(xcomb)) {
    f1 <- x[xcomb[1,i],]%>%str_split(pattern = ";")%>%lapply(function(y){
      y[y==""]<-NA
      y})
    f2 <- x[xcomb[2,i],]%>%str_split(pattern = ";")%>%lapply(function(y){
      y[y==""]<-NA
      y})
    xcomb.eq[i] <- sapply(seq_along(f1),function(z){
      all(f1[[z]] == f2[[z]],na.rm = T)
    })%>%all(na.rm = T)
    #xcomb.eq[i] <-all(x[xcomb[1,i],]==x[xcomb[2,i],],na.rm = T)

  }
  x.group <- get_igraph_membership(igraph::graph_from_data_frame(
    t(xcomb[,xcomb.eq]),vertices = 1:nrow(x)))
  return(x.group)
}





get_xcms_feature_connect <- function(xcms.xcms,rt.tol = 5){


  xcms.fdf <- featureDefinitions(xcms.xcms)

  {
    rt <- xcms.fdf$rtmed
    x <- data.table(rt = rt, start = rt,end = rt)
    y <- x[,.(id = seq_along(rt),rt,start = rt - rt.tol, end = rt + rt.tol)]
    data.table::setkey(y,start,end)
    rtm <- data.table::foverlaps(x, y, type="any", which=TRUE)
    rtm <- rtm[,.(xid,yid = y$id[yid])]
  }


  {
    xcms.net <- rtm[,.(from = xid, to = yid)]
    xcms.net <- xcms.net[from < to ][
      , rt.diff := (xcms.fdf$rtmed[to]-xcms.fdf$rtmed[from]) ][
        abs(rt.diff) < rt.tol,][
          , c("from.mz","to.mz") := .( xcms.fdf$mzmed[from], xcms.fdf$mzmed[to])][
            ,c("mz.diff","mz.mean") := .(to.mz-from.mz, (from.mz+to.mz)/2)]
  }


  return(xcms.net)



}



chemform_simplify <- function(chemform){

  ele.matrix <- MSCC::chemform_parse(chemform)
  ele.matrix <- ele.matrix[,order(colnames(ele.matrix))]
  MSCC:::chemform_from_ele_matrix(ele.matrix)

}


chemform_remove_iso <- function(chemform){

  ele.matrix <- MSCC::chemform_parse(chemform)
  iso.idx <- colnames(ele.matrix)[is.isotope(colnames(ele.matrix))]
  ele.matrix[,get_ele_uniso(iso.idx)] <- ele.matrix[,get_ele_uniso(iso.idx),drop = F]+
    ele.matrix[,iso.idx,drop = F]
  ele.matrix <- ele.matrix[,!is.isotope(colnames(ele.matrix)),drop = F]
  MSCC:::chemform_from_ele_matrix(ele.matrix)

}





get_adduct_mass_diff <- function(polarity = 0,direction = 1){


  pol <- ifelse(polarity==1,"positive","negative")

  adduct.table <- MSCC::adduct.table%>%
    dplyr::filter( Ion_mode == pol)%>%
    dplyr::mutate(m_c = Multi/Charge)


  adduct.diff <- expand.grid(
    adduct.from = 1:nrow(adduct.table),
    adduct.to = 1:nrow(adduct.table)
  )%>%
    dplyr::filter(
      adduct.table$m_c[adduct.from] == adduct.table$m_c[adduct.to]
    )%>%
    dplyr::mutate(
      chemform_diff = MSCC::chemform_calc(adduct.table$Formula_diff[adduct.to],
                                          adduct.table$Formula_diff[adduct.from],
                                          calc = "-",return = "chemform"),
      chemform_diff = chemform_simplify(chemform_diff),
      mass_diff = MSCC::chemform_mz(chemform_diff),
      #mass_diff =adduct.table$Mass[adduct.to] - adduct.table$Mass[adduct.from],
      #charge = adduct.table$Charge[adduct.to],
      adduct.from = adduct.table$Adduct[adduct.from],
      adduct.to = adduct.table$Adduct[adduct.to]
    )

  #which(upper.tri(diag(nrow(adduct.table)),diag = F),arr.ind = T)
  adduct.diff <- data.table::as.data.table(adduct.diff)

  return(adduct.diff)

}


get_iso_mass_diff <- function(){


  iso.ele <- c("[13]C", "[2]H", "[18]O","[15]N","[34]S",
               #"[41]K","[44]Ca","[10]B","[29]Si","[30]Si","[53]Cr", "[60]Ni","[62]Ni"
               "[37]Cl","[81]Br"
  )

  ele <- MSCC::elem_table%>%
    dplyr::mutate(ele.base = get_ele_uniso(element))%>%
    dplyr::group_by(ele.base)%>%
    dplyr::filter(any(element %in% iso.ele))%>%
    dplyr::arrange(mass)%>%
    dplyr::mutate(chemform_diff = paste0(element,"1",element[1],"-1"),
                  mass_diff = MSCC::chemform_mz(chemform_diff))%>%
    dplyr::ungroup()%>%
    dplyr::filter(mass_diff!=0)%>%
    dplyr::select("element","chemform_diff","mass_diff")

  data.table::as.data.table(ele)


}


get_fragment_mass_diff <- function(){

  data.table(chemform_diff = c(
    "C1O2","C1H2O1","H2O1","N1H3",### from netID
    "H1N1O1","C1H2","C2H4","C3H6" ### from data
  ))[,mass_diff := MSCC::chemform_mz(chemform_diff)][
    ,chemform_diff := chemform_simplify(chemform_diff)
  ][
    ,fragment := chemform_diff ]

}
