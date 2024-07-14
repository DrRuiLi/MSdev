get_MSIPCoreData <- function(sp.iso,
                             cfmd,
                             iso_count,
                             ppm = 10){

  sp.frag.data <- CFM_annotate_isotopologues(sp.iso,
                                             cfmd  = cfmd,
                                             ppm = ppm,
                                             iso_count = iso_count)
  sp.frag.data <- CFM_spectra_data_int_weight(sp.frag.data,iso_count)
  fg.map <- get_MSIPFragmentMap(sp.frag.data,cfmd,iso_count = iso_count)


  MSIPCoreData <- new("MSIPCoreData")
  MSIPCoreData@Spectra_data <- sp.frag.data
  MSIPCoreData@FG_map <- fg.map

  return(MSIPCoreData)


}


MSIPCore_correct_natural <- function(MSIPCoreData,
                                     cfmd,
                                     natural.ratio){

  if (is.null(MSIPCoreData@solve$MSIPIsoformMap))
    MSIPCoreData@solve$MSIPIsoformMap <- get_MSIPIsoformMap(MSIPCoreData)

  sp.frag.data <- CFM_spectra_data_remove_natural(
    sp.data = MSIPCoreData@Spectra_data,
    if.map = MSIPCoreData@solve$MSIPIsoformMap,
    natural.ratio = natural.ratio)

  fg.map <- get_MSIPFragmentMap(sp.frag.data,
                                cfmd,
                                iso_count = max(str_extract_num(colnames(MSIPCoreData@FG_map@fragment.ratio.matrix))))

  MSIPCoreData@Spectra_data <- sp.frag.data
  MSIPCoreData@FG_map <- fg.map
  return(MSIPCoreData)

}


MSIPCore_solve <- function(MSIPCoreData,
                           max_prob_map = F,
                           int_thresh = 1e3){

  MSIPFragmentMap_reduced <- MSIPFragmentMap_include_fragment(MSIPCoreData@FG_map)
  #MSIPFragmentMap_reduced <- MSIPFragmentMap_reduce_fragment(MSIPFragmentMap_reduced)
  #MSIPFragmentMap_reduced <- MSIPFragmentMap_filter_fragment(MSIPFragmentMap_reduced,int_thresh = int_thresh)

  if (isEmpty(MSIPFragmentMap_reduced))
    return(MSIPCoreData)
  if (max_prob_map) {

  }
  MSIPCoreData.temp <- MSIPCoreData
  MSIPCoreData.temp@FG_map <- MSIPFragmentMap_reduced
  MSIPIsoformMap <- get_MSIPIsoformMap(MSIPCoreData.temp)
  MSIPIsoformMap <- MSIPIsoformMap_set_split(MSIPIsoformMap,
                                       MSIPFragmentMap_reduced)
  #MSIPIsoformMap <- MSIPIsoformMap_set_solve_GLPK(MSIPIsoformMap)
  MSIPIsoformMap <- MSIPIsoformMap_set_solve_QP(MSIPIsoformMap)
  MSIPCoreData@solve$MSIPIsoformMap <- MSIPIsoformMap
  MSIPCoreData@solve$Atom_prob <- get_atom_prob_from_MSIPIsoformMap(MSIPIsoformMap)

  return(MSIPCoreData)
}

MSIPCore_drop <- function(MSIPCoreData){

  solve <- MSIPCoreData@solve
  solve <- solve["Atom_prob"]
  solve -> MSIPCoreData@solve
  return(MSIPCoreData)
}

get_MSIPFragmentMap <- function(sp.frag.data,
                                cfmd,
                                iso_count){

  ### frag group to label fraction
  {
    sp.frag.data <- sp.frag.data%>%
      dplyr::filter(sp.id== "combined_sp")
    if (nrow(sp.frag.data)==0) return(new("MSIPFragmentMap"))
    fg.idx <- split(1:nrow(sp.frag.data),sp.frag.data$fragment_group)
    frag.ratio.matrix <- matrix(
      nrow = length(fg.idx),ncol = iso_count+1,
      dimnames = list(names(fg.idx),paste0("M",0:iso_count)))
    frag.int.sum <- c()
    for (i.fg in seq_along(fg.idx)) {
      x.df <- sp.frag.data[fg.idx[[i.fg]],]
      x.int <- x.df%>%
        dplyr::select(-mz,-collisionEnergy)%>%
        tidyr::pivot_wider(names_from ="iso",
                           id_cols = "sp.id",
                           values_from = "intensity",
                           values_fn = sum)%>%
        tibble::column_to_rownames("sp.id")%>%
        dplyr::select(dplyr::starts_with("M"))%>%
        as.matrix()
      to.add <- setdiff(paste0("M",0:iso_count),colnames(x.int))
      x.int <- cbind(matrix(0,nrow(x.int),length(to.add),
                            dimnames = list(NULL,to.add)),x.int)
      x.int <- x.int[,paste0("M",0:iso_count),drop =F]
      x.int[is.na(x.int)] <- 0
      x.weight <- rowSums(x.int)
      x.int <- t(apply(x.int,1,function(z) z/sum(z)))
      x.int.weighted <- apply(x.int,2,weighted.mean,w = x.weight)
      frag.ratio.matrix[i.fg,] <- x.int.weighted
      frag.int.sum[i.fg] <- sum(x.weight)
    }
    names(frag.int.sum) <-names(fg.idx)
  }


  ### frag group to C atom prob
  {
    c_ele <- get_sdf_igraph_atom(get_cfm_data_sdf_igraph(cfmd),"C")
    frag.atom.matrix <- matrix(ncol = length(c_ele),
                            nrow = length(fg.idx),
                            dimnames = list(names(fg.idx),
                                            c_ele))
    for (i.fg in seq_along(fg.idx)) {

      this.frag.group <- names(fg.idx)[i.fg]
      this.frags <- cfmd@fragment_define[cfmd@fragment_define$fragment_group==this.frag.group,]
      this.frag.ratio <-frag.ratio.matrix[i.fg,]
      this.frag.atom <- get_cfm_data_fragment_group_atom_map(cfmd,this.frag.group)
      this.frag.c <- this.frag.atom[c_ele]
      #this.iso.expectation <- sum(str_extract_num(names(this.frag.ratio))*this.frag.ratio)
      #this.frag.c <- this.frag.c*this.iso.expectation/sum(this.frag.c)
      #this.frag.c <- this.frag.c[this.frag.c!=0]
      frag.atom.matrix[this.frag.group,names(this.frag.c)] <- this.frag.c
    }

  }


  ### return
  {
    #if (atom_prob) {
#
    #}else{
    #  frag.atom.matrix.max <- apply(frag.atom.matrix,1,function(x){
    #    x.sum <- sum(x)
    #    x.new <- rep(0,length(x))
    #    x.new[tail(order(x),x.sum)] <- 1
    #    x.new
    #  })%>%t
    #  colnames(frag.atom.matrix.max) <- colnames(frag.atom.matrix)
    #  frag.atom.matrix <- frag.atom.matrix.max
    #}

    fg.map <- new("MSIPFragmentMap")
    fg.map@fragment.atom.matrix <- frag.atom.matrix
    fg.map@fragment.ratio.matrix <- frag.ratio.matrix
    fg.map@fragment.intensity <- frag.int.sum
    fg.map@fragment.include <- make_vector(T, names(frag.int.sum))

    return(fg.map)

  }




}

get_MSIPIsoformMap <- function(MSIPCoreData){



  ### required info
  {

    frag.atom.matrix <- MSIPCoreData@FG_map@fragment.atom.matrix
    frag.ratio.matrix <-MSIPCoreData@FG_map@fragment.ratio.matrix
    frag.max.iso <- ncol(frag.ratio.matrix)-1
    MSIPIsoformMap <- new("MSIPIsoformMap")
    if (!nrow(frag.atom.matrix)) {
      return(MSIPIsoformMap)
    }

  }
  ### all possible iso form
  {
    if.combn <- choose(ncol(frag.atom.matrix),frag.max.iso)
    iso.form <- combn(colnames(frag.atom.matrix),frag.max.iso,simplify = F)
    names(iso.form) <- paste0("iso_form_",num2str(1:length(iso.form)))
  }
  ### iso form map to iso ratio
  {
    iso.form.maps <- lapply(seq_along(iso.form),
                            function(if.id){
                              lapply(rownames(frag.atom.matrix),function(fg.id){

                                get_iso_prob(frag.atom.matrix[fg.id,],
                                             iso.form[[ if.id ]])
                              })->mp
                              names(mp) <- rownames(frag.atom.matrix)
                              unlist(mp)
                            })
    iso.form.map <- t(do.call(bind_rows,iso.form.maps))
    iso.form.map <- iso.form.map[order(rownames(iso.form.map)),,drop = F]
    iso.form.map[is.na(iso.form.map)] <- 0
    rownames(iso.form.map) <- sub(pattern = ".",x = rownames(iso.form.map),
                                  replacement = "_",fixed = T)
    colnames(iso.form.map) <- names(iso.form)
  }

  ### iso form ratio
  {

    iso.form.ratio <- lapply(rownames(frag.ratio.matrix),function(fg.id){
      z <- frag.ratio.matrix[fg.id,]
      names(z) <- paste0(fg.id,"_",names(z))
      return(z)
    })
    iso.form.ratio <- unlist(iso.form.ratio)

  }

  ### iso form intensity
  {
    iso.form.intensity <- MSIPCoreData@FG_map@fragment.intensity[gsub(x = names(iso.form.ratio),
                                                                          pattern = "_M.*",replacement = "")]
    names(iso.form.intensity) <- names(iso.form.ratio)
  }



  MSIPIsoformMap@isoform.defination <-iso.form
  MSIPIsoformMap@isoform.map <-iso.form.map
  MSIPIsoformMap@isoform.ratio <-iso.form.ratio
  MSIPIsoformMap@isoform.intensity <-iso.form.intensity
  return(MSIPIsoformMap)

}



get_MSIPIsoformMap_from_atom <- function(atom,iso_count){

  x <- new("MSIPIsoformMap")

  exp <- paste0( "combn(",
          vector2str(atom),",",
          iso_count,",","simplify = F)"
  )
  exp <- str2expression(exp)
  x@isoform.defination <- eval(exp)

  return(x)

}



MSIPFragmentMap_reduce_fragment <- function(MSIPFragmentMap){


  MSIPFragmentMap <- MSIPFragmentMap_merge_duplicate(MSIPFragmentMap)
  MSIPFragmentMap <- MSIPFragmentMap_merge_complementary(MSIPFragmentMap)


  return(MSIPFragmentMap)
}

MSIPFragmentMap_include_fragment <- function(MSIPFragmentMap){

  frag.include <- MSIPFragmentMap@fragment.include
  frag.include <- frag.include[frag.include]
  frag.include <- names(frag.include)

  MSIPFragmentMap@fragment.atom.matrix <- MSIPFragmentMap@fragment.atom.matrix[frag.include,,drop = F]
  MSIPFragmentMap@fragment.ratio.matrix <- MSIPFragmentMap@fragment.ratio.matrix[frag.include,,drop = F]
  MSIPFragmentMap@fragment.intensity <- MSIPFragmentMap@fragment.intensity[frag.include]
  MSIPFragmentMap@fragment.include <-  MSIPFragmentMap@fragment.include[frag.include]
  return(MSIPFragmentMap)
}

MSIPFragmentMap_filter_fragment <- function(MSIPFragmentMap,
                                            int_thresh = 1E3){

  frag.int <- MSIPFragmentMap@fragment.intensity
  frag.include <- frag.int>int_thresh
  frag.include <- frag.include[frag.include]
  frag.include <- names(frag.include)

  MSIPFragmentMap@fragment.atom.matrix <- MSIPFragmentMap@fragment.atom.matrix[frag.include,,drop = F]
  MSIPFragmentMap@fragment.ratio.matrix <- MSIPFragmentMap@fragment.ratio.matrix[frag.include,,drop = F]
  MSIPFragmentMap@fragment.intensity <- MSIPFragmentMap@fragment.intensity[frag.include]
  MSIPFragmentMap@fragment.include <-  MSIPFragmentMap@fragment.include[frag.include]
  return(MSIPFragmentMap)
}

MSIPFragmentMap_merge_duplicate <- function(MSIPFragmentMap){


  frag.atom.matrix <- MSIPFragmentMap@fragment.atom.matrix
  frag.ratio.matrix <-  MSIPFragmentMap@fragment.ratio.matrix
  frag.int.sum <- MSIPFragmentMap@fragment.intensity
  if (isEmpty(MSIPFragmentMap))
    return(MSIPFragmentMap)
  ### merge duplicate
  {

    z <- frag.atom.matrix
    frag.atom.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.atom.matrix))
    frag.ratio.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.ratio.matrix))
    frag.int.sum1 <- c()
    z.split <- split(1:nrow(z),apply(z,1,paste0,collapse = ";"))
    for (i.z in seq_along(z.split)) {
      idx <- z.split[[i.z]]
      this.frag.c <- apply(frag.atom.matrix[idx,,drop =F],2,
                           weighted.mean,w = frag.int.sum[idx])
      this.frag.iso <- apply(frag.ratio.matrix[idx,,drop =F],2,
                             weighted.mean,w = frag.int.sum[idx])
      frag.atom.matrix1 <- rbind(frag.atom.matrix1,this.frag.c)
      frag.ratio.matrix1 <- rbind(frag.ratio.matrix1,this.frag.iso)
      frag.int.sum1 <- c(frag.int.sum1, max(frag.int.sum[idx]) )
      rn <- rownames(frag.ratio.matrix)[idx][which.max(frag.int.sum[idx])]
      rownames(frag.atom.matrix1)[i.z] <- rownames(frag.ratio.matrix1)[i.z]  <- rn
      names(frag.int.sum1)[i.z] <- rn
    }


  }

  frag.atom.matrix1 -> MSIPFragmentMap@fragment.atom.matrix
  frag.ratio.matrix1 -> MSIPFragmentMap@fragment.ratio.matrix
  frag.int.sum1 ->  MSIPFragmentMap@fragment.intensity
  MSIPFragmentMap@fragment.include <- MSIPFragmentMap@fragment.include[names(frag.int.sum1)]
  return(MSIPFragmentMap)
}


MSIPFragmentMap_merge_complementary <- function(MSIPFragmentMap){

  if (isEmpty(MSIPFragmentMap))
    return(MSIPFragmentMap)

  frag.atom.matrix <- MSIPFragmentMap@fragment.atom.matrix
  frag.ratio.matrix <-  MSIPFragmentMap@fragment.ratio.matrix
  frag.int.sum <- MSIPFragmentMap@fragment.intensity

  ### merge duplicate
  {

    ### complementary
    z <- frag.atom.matrix
    #z[z>0] <- 1
    z.comple <- apply(z,1,function(x){
      z1 <- t(t(z)+x)
      apply(z1,1,function(xx){all(xx==1)})

    })
    dim(z.comple) <- c(nrow(z),nrow(z))
    z.comple <- which(z.comple,arr.ind = T)
    z.split <- lapply(1:nrow(z),function(x){
      x.c <- z.comple[z.comple[,1] == x,2]
      x.c <- unname(x.c)
      sort( c(x,x.c))
    })
    z.split <- unique(z.split)
    frag.atom.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.atom.matrix))
    frag.ratio.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.ratio.matrix))
    frag.int.sum1 <- c()
    for (i.z in seq_along(z.split)) {
      idx <- z.split[[i.z]]
      this.frag.c <- apply(frag.atom.matrix[idx,,drop =F],2,
                           weighted.mean,w = frag.int.sum[idx])
      this.frag.iso <- apply(frag.ratio.matrix[idx,,drop =F],2,
                             weighted.mean,w = frag.int.sum[idx])
      frag.atom.matrix1 <- rbind(frag.atom.matrix1,this.frag.c)
      frag.ratio.matrix1 <- rbind(frag.ratio.matrix1,this.frag.iso)
      frag.int.sum1 <- c(frag.int.sum1, max(frag.int.sum[idx]) )
      rn <- rownames(frag.ratio.matrix)[idx][which.max(frag.int.sum[idx])]
      rownames(frag.atom.matrix1)[i.z] <- rownames(frag.ratio.matrix1)[i.z]  <- rn
      names(frag.int.sum1)[i.z] <- rn
    }


  }

  frag.atom.matrix1 -> MSIPFragmentMap@fragment.atom.matrix
  frag.ratio.matrix1 -> MSIPFragmentMap@fragment.ratio.matrix
  frag.int.sum1 ->  MSIPFragmentMap@fragment.intensity
  MSIPFragmentMap@fragment.include <- MSIPFragmentMap@fragment.include[names(frag.int.sum1)]
  return(MSIPFragmentMap)
}



#' heatmap_MSIPFragmentMap
#' @describeIn MSIPCore heatmap fragment map
#'
#' @param MSIPFragmentMap MSIPFragmentMap
#'
#' @return heatmap
#' @export
#' @import ComplexHeatmap grid
heatmap_MSIPFragmentMap <- function(MSIPFragmentMap){

  ### size
  {
    ppi <- 72
    mhpx <- 500
    mwpx <- 800
    nr <- nrow(MSIPFragmentMap@fragment.atom.matrix)
    nc <- ncol(MSIPFragmentMap@fragment.atom.matrix)+
      ncol(MSIPFragmentMap@fragment.ratio.matrix)
    cell.unit <- "inches"
    length.unit <- 0.3
    max.height <- mhpx/ppi
    max.width <- mwpx/ppi-1
    if (nr*length.unit>max.height)
      length.unit <- max.height/nc
    if (nc*length.unit>max.width)
      length.unit <- max.width/nc

  }

  frag.atom.matrix <- MSIPFragmentMap@fragment.atom.matrix
  frag.ratio.matrix <- MSIPFragmentMap@fragment.ratio.matrix
  cf <- circlize::colorRamp2(breaks = c(0,0.5,1),
                             c("white","#888888","#111111"))
  h1 <- ComplexHeatmap::Heatmap(frag.atom.matrix,
                                na_col  ="#999999",
                                # width = unit(ncol(frag.atom.matrix)*length.unit,cell.unit),
                                # height =unit(nrow(frag.atom.matrix)*length.unit,cell.unit),
                                name = "Atom map\nprobability",
                                col = cf,
                                cluster_columns = F,
                                rect_gp = grid::gpar(lwd=2,col = "black",type = "none"),
                                cell_fun = function(j, i, x, y, width, height,color, fill){
                                  grid.rect(x = x, y = y, width = width, height = height,
                                            gp = grid::gpar(col = "grey", fill = NA))
                                  grid.circle(x = x, y = y, r = min(width,height),
                                              gp = grid::gpar(fill = color, col = color))
                                },
                                row_names_side  = "left",
                                column_names_rot = -45,
                                column_names_centered = F,
                                column_names_gp = grid::gpar(fontsize = 8),
                                #rect_gp =  grid::gpar(lwd=2,col = "white"),
                                cluster_rows = F)
  h1
  h2 <- ComplexHeatmap::Heatmap(frag.ratio.matrix,
                                na_col  ="#999999",
                                # width = unit( ncol(frag.ratio.matrix)*length.unit,cell.unit),
                                # height =unit( nrow(frag.atom.matrix)*length.unit,cell.unit),
                                name = "Isotope labeled\nratio",
                                col = circlize::colorRamp2(breaks = c(0,0.5,1),
                                                           c("white","#F7844F","#B20C26")),
                                right_annotation  = rowAnnotation(
                                  intensity = anno_numeric(round(log10(MSIPFragmentMap@fragment.intensity),1),
                                                           bg_gp = gpar(fill = "#AFAFAF", col = "black")),
                                  width  = unit(0.8,"inch"),
                                  annotation_label = list(intensity = "Log10\nIntensity"),
                                  annotation_name_rot  = 0,
                                  annotation_name_side  = "top"),
                                cluster_columns = F,
                                row_names_side  = "left",
                                column_names_side = "top",
                                column_names_rot = 0.5,
                                column_names_centered = T,
                                rect_gp =  grid::gpar(lwd=2,col = "black"),
                                cluster_rows = F)
  draw(h1+h2,#padding =  unit(c(max.height-nr*length.unit, 0, 10,
       #                  max.width-nc*length.unit), "points"),
       padding = unit(c(0,0,0.2,1), "inches"),
       background  = "#00000000")
  #open_plot_win(h1+h2,10,5)
}

heatmap_MSIPIsoformMap <- function(MSIPIsoformMap){

  iso.form.set.map <- MSIPIsoformMap@isoform.map
  iso.form.prob <- MSIPIsoformMap@isoform.probability
  iso.form.ratio <- MSIPIsoformMap@isoform.ratio
  iso.form.ratio <- iso.form.ratio[rownames(iso.form.set.map)]
  top.anno <- HeatmapAnnotation(
      ifp = iso.form.prob,
      col = list(ifp = colramp(c(0,0.00000001,max(iso.form.prob)),
                               c("grey","white","#0095D4"))),
      annotation_label  = list(ifp = "iso form probability"),
      show_annotation_name = F,
      which = "c"
    )


  Heatmap(iso.form.set.map,

          border = T,
          border_gp = gpar(col = "#808080"),
          row_split = str_extract(rownames(iso.form.set.map),".*(?=_M)"),
          row_title = NULL,

          cluster_rows = F,
          show_column_names = F,
          show_column_dend = T,
          show_row_names = T,
          row_names_side  = "l",
          column_title = "Iso form set",
          column_title_gp = gpar(fontsize = 30),
          cluster_columns = T,
          top_annotation = top.anno,
          left_annotation = HeatmapAnnotation(
            Ratio = iso.form.ratio,
            col = list(Ratio = colramp(c(0,0.0000001,0.5,1),
                                       c("grey","white","#F7844F","#B20C26"))),
            which = "row",
            show_annotation_name = F,
            width  = unit(20,"inch")),
          show_heatmap_legend = F,
          col = colramp(colors = c("white","white","black")))
}


MSIPIsoformMap_set_split <- function(MSIPIsoformMap,
                                     MSIPFragmentMap_reduced){


  if (isEmpty(MSIPIsoformMap))
    return(MSIPIsoformMap)

  ### iso form ratio
  {
    frag.ratio.matrix <- MSIPFragmentMap_reduced@fragment.ratio.matrix
    iso.form.ratio <- lapply(rownames(frag.ratio.matrix),function(fg.id){
      z <- frag.ratio.matrix[fg.id,]
      names(z) <- paste0(fg.id,"_",names(z))
      return(z)
    })
    iso.form.ratio <- unlist(iso.form.ratio)
  }


  ### iso form intensity
  {
    iso.form.intensity <- MSIPFragmentMap_reduced@fragment.intensity[gsub(x = names(iso.form.ratio),
         pattern = "_M.*",replacement = "")]
    names(iso.form.intensity) <- names(iso.form.ratio)
  }

  ### iso form map
  {

    iso.form.map <-MSIPIsoformMap@isoform.map
    to.add <- matrix(0,
                     nrow = length(setdiff(names(iso.form.ratio),
                                           rownames(iso.form.map))),
                     ncol = ncol(iso.form.map),
                     dimnames = list(rowname = setdiff(names(iso.form.ratio),
                                                       rownames(iso.form.map))))
    iso.form.map <- rbind(iso.form.map,to.add)
    iso.form.map <- iso.form.map[names(iso.form.ratio), ,drop = F]

  }

  ### iso form set
  {
    iso.form.split <- apply(iso.form.map, 2, function(x){paste0(x,collapse = ";")})
    iso.form.split <-  split(seq_along(MSIPIsoformMap@isoform.defination),iso.form.split)
    names(iso.form.split) <- paste0("iso_form_set_",num2str(seq_along(iso.form.split)))
    iso.form.set.map <- sapply(iso.form.split ,
                               function(x){ iso.form.map[,x[1]] },
                               USE.NAMES=F)
  }


  MSIPIsoformMap@solve$isoform.set <- iso.form.split
  MSIPIsoformMap@solve$isoform.set.map <- iso.form.set.map
  MSIPIsoformMap@solve$isoform.set.ratio <- iso.form.ratio
  MSIPIsoformMap@solve$isoform.set.intensity<- iso.form.intensity

  return(MSIPIsoformMap)
}

MSIPIsoformMap_set_solve_GLPK <- function(MSIPIsoformMap){

  if (isEmpty(MSIPIsoformMap))
    return(MSIPIsoformMap)

### GLPK
  {
  mat <- MSIPIsoformMap@solve$isoform.set.map
  obj <- rep(1,ncol(mat))
  dir <- rep(">=",nrow(mat))
  rhs <- MSIPIsoformMap@solve$isoform.set.ratio

  lp.result <- Rglpk::Rglpk_solve_LP(obj = obj,mat = mat,
                                     canonicalize_status = T,
                                     bounds = list(lower = list(ind = 1:ncol(mat),
                                                                val = rep(0,ncol(mat))),
                                                   upper = list(ind = 1:ncol(mat),
                                                                val = rep(1,ncol(mat)))),
                                     types = rep("C",ncol(mat)),
                                     dir = dir,
                                     rhs = rhs,
                                     max = F)
  lp.result$status
  lp.result$solution

  }

  ### assign set to isoform
  {
    isoform.set.prob <- lp.result$solution
    isoform.set <- MSIPIsoformMap@solve$isoform.set
    isoform.prob <- mapply(x = isoform.set.prob ,
                           y = isoform.set,function(x,y){
      make_vector(x/length(y),num2str(y,10))
    },SIMPLIFY = F)
    isoform.prob <- unlist(isoform.prob)
    isoform.prob <- isoform.prob[order(names(isoform.prob))]
    isoform.prob <- unname(isoform.prob)
  }
  MSIPIsoformMap@solve$isoform.set.prob <- lp.result$solution
  MSIPIsoformMap@solve$Rglpk <- lp.result
  MSIPIsoformMap@isoform.probability <- isoform.prob
  return(MSIPIsoformMap)
}

MSIPIsoformMap_set_solve_QP <- function(MSIPIsoformMap){

  if (isEmpty(MSIPIsoformMap))
    return(MSIPIsoformMap)

  ###  Quadratic Programming
  {

    # Load necessary library

    # Example data (I_i: observed isotopomer intensities, f_ij: contribution matrix)
    I <- MSIPIsoformMap@solve$isoform.set.ratio # 3 fragments
    f <- MSIPIsoformMap@solve$isoform.set.map # 4 atoms

    # Define weights based on fragment reliability (higher intensity = higher weight)
    int <- MSIPIsoformMap@solve$isoform.set.intensity # Example weights
    weights <- log10(int)/10
    weights <- weights^4
    # Incorporate weights into the contribution matrix and observed intensities
    W <- diag(weights)
    f_weighted <- W %*% f
    I_weighted <- W %*% I

    # Number of atoms
    n_atoms <- ncol(f_weighted)

    # Construct the matrices for the quadratic programming problem
    Dmat <- t(f_weighted) %*% f_weighted

    # Regularization term to make Dmat positive definite
    epsilon <- 1e-6
    Dmat <- Dmat + epsilon * diag(n_atoms)

    dvec <- t(f_weighted) %*% I_weighted

    # Constraints: probabilities should be between 0 and 1
    Amat <- cbind(diag(n_atoms), -diag(n_atoms))
    bvec <- c(rep(0, n_atoms), rep(-1, n_atoms))

    # Solve the quadratic programming problem
    qp_result <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 0)

    # Estimated labeling probabilities
    p_estimated <- qp_result$solution
    round(p_estimated,2)


    ratio.pred <- f%*%p_estimated
    ratio.df <- data.frame(
      val =round(I,4)  ,
      pred = round(ratio.pred,4),
      int = log10(int)
    )
    }


  ### assign set to isoform
  {
    isoform.set.prob <- p_estimated
    isoform.set <- MSIPIsoformMap@solve$isoform.set
    isoform.prob <- mapply(x = isoform.set.prob ,
                           y = isoform.set,function(x,y){
                             make_vector(x/length(y),num2str(y,10))
                           },SIMPLIFY = F)
    isoform.prob <- unlist(isoform.prob)
    isoform.prob <- isoform.prob[order(names(isoform.prob))]
    isoform.prob <- unname(isoform.prob)
  }
  MSIPIsoformMap@solve$isoform.set.prob <- isoform.set.prob
  MSIPIsoformMap@solve$solve.QP <- qp_result
  MSIPIsoformMap@isoform.probability <- isoform.prob
  return(MSIPIsoformMap)
}


get_atom_prob_from_MSIPIsoformMap <- function(MSIPIsoformMap){

  if (isEmpty(MSIPIsoformMap))
    return(NULL)


  atom <- unique(unlist(MSIPIsoformMap@isoform.defination))
  atom.prob <- make_vector(0,atom)
  iso_count <- length(MSIPIsoformMap@isoform.defination[[1]])
  atom.prob.m <- sapply(seq_along(MSIPIsoformMap@isoform.probability),
                    function(idx){

    x <- MSIPIsoformMap@isoform.probability[idx]
    x.atom <- MSIPIsoformMap@isoform.defination[[idx]]
    x.prob <- atom.prob
    x.prob[x.atom] <-x/iso_count

    return(x.prob)
  })

  x <- apply(atom.prob.m,1,sum)
  if (all(x==0)) {
    return(x)
  }
  x/sum(x)*iso_count

}



MSIPCore_vis_isoform <- function(MSIPCoreData,cfmd,id,...){

  ig <- get_cfm_data_sdf_igraph(cfmd)
  atom.to.show <- MSIPCoreData@solve$MSIPIsoformMap@isoform.defination[[id]]
  message(names(MSIPCoreData@solve$MSIPIsoformMap@isoform.defination)[id])
  vis_sdf_igraph(ig,highlight = atom.to.show,...)

}


