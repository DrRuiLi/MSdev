get_MSIPCoreData <- function(sp.iso,
                             cfmd,
                             iso_count,
                             ppm = 10){

  sp.frag.data <- CFM_annotate_isotopologues(sp.iso,
                                             cfmd  = cfmd,
                                             ppm = ppm,
                                             iso_count = iso_count)
  sp.frag.data <- CFM_spectra_data_merge(sp.frag.data,iso_count)
  fg.map <- get_MSIPFragmentMap(sp.frag.data,
                                cfmd,
                                iso_count = iso_count)


  MSIPCoreData <- new("MSIPCoreData")
  MSIPCoreData@Spectra_data <- sp.frag.data
  MSIPCoreData@FG_map <- fg.map

  return(MSIPCoreData)


}


MSIPCore_correct_natural <- function(MSIPCoreData,
                                     cfmd,
                                     natural.ratio){

  if (is.null(MSIPCoreData@solve$MSIPIsotopomerMap))
    MSIPCoreData@solve$MSIPIsotopomerMap <- get_MSIPIsotopomerMap(MSIPCoreData)

  sp.frag.data <- CFM_spectra_data_remove_natural(
    sp.data = MSIPCoreData@Spectra_data,
    if.map = MSIPCoreData@solve$MSIPIsotopomerMap,
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
                           int_thresh = 10^3.8,
                           certainty_thresh = 0.8,
                           re_split_isotopomers = T){

  if (isEmpty(MSIPCoreData))
    return(MSIPCoreData)

  ### set all include
  MSIPCoreData@FG_map@fragment.include <- MSIPCoreData@FG_map@fragment.include|T
  #MSIPFragmentMap_reduced <- MSIPFragmentMap_reduce_fragment(MSIPFragmentMap_reduced)
  MSIPFragmentMap_temp <- MSIPFragmentMap_filter_intensity(MSIPCoreData@FG_map,int_thresh = int_thresh)
  MSIPFragmentMap_temp <- MSIPFragmentMap_filter_certainty(MSIPFragmentMap_temp,certainty_thresh = certainty_thresh)
  MSIPCoreData@FG_map <- MSIPFragmentMap_temp

  MSIPFragmentMap_reduced <- MSIPFragmentMap_include_fragment(MSIPFragmentMap_temp)%>%
    MSIPFragmentMap_add_constraint()

  if (isEmpty(MSIPFragmentMap_reduced))
    return(MSIPCoreData)
  if (max_prob_map) {

  }
  MSIPCoreData.temp <- MSIPCoreData
  MSIPCoreData.temp@FG_map <- MSIPFragmentMap_reduced
  if(re_split_isotopomers)
    MSIPCoreData.temp@solve <- list()
  MSIPIsotopomerMap <- get_MSIPIsotopomerMap(MSIPCoreData.temp)
  MSIPIsotopomerMap <- MSIPIsotopomerMap_set_split(MSIPIsotopomerMap,
                                       MSIPFragmentMap_reduced)
  #MSIPIsotopomerMap <- MSIPIsotopomerMap_set_solve_GLPK(MSIPIsotopomerMap)
  MSIPIsotopomerMap <- MSIPIsotopomerMap_set_solve_QP(MSIPIsotopomerMap)
  MSIPCoreData@solve$MSIPIsotopomerMap <- MSIPIsotopomerMap
  MSIPCoreData@solve$Atom_prob <- get_atom_prob_from_MSIPIsotopomerMap(MSIPIsotopomerMap)
  MSIPCoreData@solve$int_thresh <- int_thresh
  MSIPCoreData@solve$certainty_thresh <- certainty_thresh
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
                                iso_ele = "[13]C",
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
        tidyr::pivot_wider(names_from ="iso_count",
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


  ### frag group to atom prob
  {
    target_atoms <- get_sdf_igraph_atom(get_cfm_data_sdf_igraph(cfmd),get_ele_uniso(iso_ele))
    frag.atom.matrix <- matrix(ncol = length(target_atoms),
                            nrow = length(fg.idx),
                            dimnames = list(names(fg.idx),
                                            target_atoms))
    for (i.fg in seq_along(fg.idx)) {

      this.frag.group <- names(fg.idx)[i.fg]
      this.frags <- cfmd@fragment_define[cfmd@fragment_define$fragment_group==this.frag.group,]
      this.frag.ratio <-frag.ratio.matrix[i.fg,]
      this.frag.atom <- get_cfm_data_fragment_group_atom_map(cfmd,this.frag.group)
      this.frag.c <- this.frag.atom[target_atoms]
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

get_MSIPIsotopomerMap <- function(MSIPCoreData){




  ### filter
  {
    if (!is.null(MSIPCoreData@solve$MSIPIsotopomerMap)) {
      MSIPIsotopomerMap <- MSIPCoreData@solve$MSIPIsotopomerMap
      frag.ratio.matrix <- MSIPCoreData@FG_map@fragment.ratio.matrix
      isotopomer.ratio <- lapply(rownames(frag.ratio.matrix),function(fg.id){
        z <- frag.ratio.matrix[fg.id,]
        names(z) <- paste0(fg.id,"_",names(z))
        return(z)
      })
      isotopomer.ratio <- unlist(isotopomer.ratio)
      isotopomer.intensity <- MSIPCoreData@FG_map@fragment.intensity[gsub(x = names(isotopomer.ratio),
                                                                        pattern = "_M.*",replacement = "")]
      names(isotopomer.intensity) <- names(isotopomer.ratio)

      #MSIPIsotopomerMap@isotopomer.map <- MSIPIsotopomerMap@isotopomer.map[names(isotopomer.ratio),,drop = F]
      #MSIPIsotopomerMap@isotopomer.ratio <- MSIPIsotopomerMap@isotopomer.ratio[names(isotopomer.ratio)]
      #MSIPIsotopomerMap@isotopomer.intensity <- isotopomer.intensity
      return(MSIPIsotopomerMap)

    }

  }

  ### required info
  {

    frag.atom.matrix <- MSIPCoreData@FG_map@fragment.atom.matrix
    frag.ratio.matrix <-MSIPCoreData@FG_map@fragment.ratio.matrix
    frag.max.iso <- ncol(frag.ratio.matrix)-1
    MSIPIsotopomerMap <- new("MSIPIsotopomerMap")
    if (!nrow(frag.atom.matrix)) {
      return(MSIPIsotopomerMap)
    }

  }
  ### all possible iso form
  {
    if (frag.max.iso == 0 ) {
      message_with_time("No isotopomers")
      return(MSIPIsotopomerMap)
    }
    if.combn <- choose(ncol(frag.atom.matrix),frag.max.iso)
    if(if.combn > 1e8){
      message_with_time("Isotopomers too many: ",format(1236547889,sci = T,digits = 2),", solve cancel")
      return(MSIPIsotopomerMap)
    }
    isotopomer <- combn(colnames(frag.atom.matrix),frag.max.iso,simplify = F)
    names(isotopomer) <- paste0("isotopomer_",num2str(1:length(isotopomer)))
  }
  ### iso form map to iso ratio
  {
    isotopomer.maps <- bplapply(seq_along(isotopomer),
                            function(if.id){
                              #message_with_time(if.id)
                              lapply(rownames(frag.atom.matrix),function(fg.id){

                                get_iso_prob_chatgpt(frag.atom.matrix[fg.id,],
                                             isotopomer[[ if.id ]])
                              })->mp
                              names(mp) <- rownames(frag.atom.matrix)
                              unlist(mp)
                            },BPPARAM = SerialParam(progressbar = F))

    isotopomer.map <- t(do.call(bind_rows,isotopomer.maps))
    isotopomer.map <- isotopomer.map[order(rownames(isotopomer.map)),,drop = F]
    isotopomer.map[is.na(isotopomer.map)] <- 0
    rownames(isotopomer.map) <- sub(pattern = ".",x = rownames(isotopomer.map),
                                  replacement = "_",fixed = T)
    colnames(isotopomer.map) <- names(isotopomer)
  }

  ### iso form ratio
  {

    isotopomer.ratio <- lapply(rownames(frag.ratio.matrix),function(fg.id){
      z <- frag.ratio.matrix[fg.id,]
      names(z) <- paste0(fg.id,"_",names(z))
      return(z)
    })
    isotopomer.ratio <- unlist(isotopomer.ratio)

  }

  ### iso form intensity
  {
    isotopomer.intensity <- MSIPCoreData@FG_map@fragment.intensity[gsub(x = names(isotopomer.ratio),
                                                                          pattern = "_M.*",replacement = "")]
    names(isotopomer.intensity) <- names(isotopomer.ratio)
  }



  MSIPIsotopomerMap@isotopomer.defination <-isotopomer
  MSIPIsotopomerMap@isotopomer.map <-isotopomer.map
  MSIPIsotopomerMap@isotopomer.ratio <-isotopomer.ratio
  MSIPIsotopomerMap@isotopomer.intensity <-isotopomer.intensity
  return(MSIPIsotopomerMap)

}



get_MSIPIsotopomerMap_from_atom <- function(atom,iso_count){

  x <- new("MSIPIsotopomerMap")

  exp <- paste0( "combn(",
          vector2str(atom),",",
          iso_count,",","simplify = F)"
  )
  exp <- str2expression(exp)
  x@isotopomer.defination <- eval(exp)

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

MSIPFragmentMap_filter_intensity <- function(MSIPFragmentMap,
                                            int_thresh = 1E3){

  frag.int <- MSIPFragmentMap@fragment.intensity
  frag.include <- frag.int > int_thresh

  frag.include <- frag.include & MSIPFragmentMap@fragment.include
  MSIPFragmentMap@fragment.include <-  frag.include

  return(MSIPFragmentMap)
}

MSIPFragmentMap_filter_certainty <- function(MSIPFragmentMap,
                                            certainty_thresh = 0.8 ){

  frag.certainty <- get_MSIPFragmentMap_certainty(MSIPFragmentMap)
  frag.include <- frag.certainty >= certainty_thresh

  frag.include <- frag.include & MSIPFragmentMap@fragment.include

  MSIPFragmentMap@fragment.include <-  frag.include
  return(MSIPFragmentMap)
}
MSIPFragmentMap_add_constraint <- function(MSIPFragmentMap){

  m <- MSIPFragmentMap@fragment.atom.matrix
  MSIPFragmentMap@fragment.atom.matrix <-
    rbind(m,matrix(c(0,1),2,ncol(m),dimnames = list(c("FGNULL","FGFULL"),colnames(m))))

  m <- MSIPFragmentMap@fragment.ratio.matrix
  m1 <- matrix(0,2,ncol(m),dimnames = list(c("FGNULL","FGFULL"),colnames(m)))
  m1[1] <- m1[length(m1)] <- 1
  MSIPFragmentMap@fragment.ratio.matrix <- rbind(m,m1)

  MSIPFragmentMap@fragment.intensity %<>% c(.,FGNULL = 1e9,FGFULL=1e9)
  MSIPFragmentMap@fragment.include %<>% c(.,FGNULL = T,FGFULL=T)


  return(MSIPFragmentMap)

}

get_MSIPFragmentMap_certainty <- function(MSIPFragmentMap){
  frag.matrix <- MSIPFragmentMap@fragment.atom.matrix
  frag.certainty <- apply(frag.matrix,1,function(x){
    sum(x==1)/sum(x)
  })
  frag.certainty[rowSums(frag.matrix)==0] <- 0
  return(frag.certainty)
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
  frag.include <- MSIPFragmentMap@fragment.include
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
                                  grid.circle(x = x, y = y, r = min(width,height)/2,
                                              gp = grid::gpar(fill = color, col = color))
                                },
                                row_names_side  = "left",
                                column_names_rot = -45,
                                column_names_centered = F,
                                column_names_gp = grid::gpar(fontsize = 8),
                                row_names_gp = gpar(col = ifelse(frag.include,"black","grey")),
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

heatmap_MSIPIsotopomerMap <- function(MSIPIsotopomerMap){

  isotopomer.set.map <- MSIPIsotopomerMap@isotopomer.map
  isotopomer.prob <- MSIPIsotopomerMap@isotopomer.probability
  isotopomer.ratio <- MSIPIsotopomerMap@isotopomer.ratio
  isotopomer.ratio <- isotopomer.ratio[rownames(isotopomer.set.map)]
  top.anno <- HeatmapAnnotation(
      ifp = isotopomer.prob,
      col = list(ifp = colramp(c(0,0.00000001,max(isotopomer.prob)),
                               c("grey","white","#0095D4"))),
      annotation_label  = list(ifp = "iso form probability"),
      show_annotation_name = F,
      which = "c"
    )


  Heatmap(isotopomer.set.map,

          border = T,
          border_gp = gpar(col = "#808080"),
          row_split = str_extract(rownames(isotopomer.set.map),".*(?=_M)"),
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
            Ratio = isotopomer.ratio,
            col = list(Ratio = colramp(c(0,0.0000001,0.5,1),
                                       c("grey","white","#F7844F","#B20C26"))),
            which = "row",
            show_annotation_name = F,
            width  = unit(20,"inch")),
          show_heatmap_legend = F,
          col = colramp(colors = c("white","white","black")))
}


MSIPIsotopomerMap_set_split <- function(MSIPIsotopomerMap,
                                     MSIPFragmentMap_reduced){


  if (isEmpty(MSIPIsotopomerMap))
    return(MSIPIsotopomerMap)

  ### iso form ratio
  {
    frag.ratio.matrix <- MSIPFragmentMap_reduced@fragment.ratio.matrix
    isotopomer.ratio <- lapply(rownames(frag.ratio.matrix),function(fg.id){
      z <- frag.ratio.matrix[fg.id,]
      names(z) <- paste0(fg.id,"_",names(z))
      return(z)
    })
    isotopomer.ratio <- unlist(isotopomer.ratio)
  }


  ### iso form intensity
  {
    isotopomer.intensity <- MSIPFragmentMap_reduced@fragment.intensity[gsub(x = names(isotopomer.ratio),
         pattern = "_M.*",replacement = "")]
    names(isotopomer.intensity) <- names(isotopomer.ratio)
  }

  ### iso form map
  {

    isotopomer.map <-MSIPIsotopomerMap@isotopomer.map
    to.add <- matrix(0,
                     nrow = length(setdiff(names(isotopomer.ratio),
                                           rownames(isotopomer.map))),
                     ncol = ncol(isotopomer.map),
                     dimnames = list(rowname = setdiff(names(isotopomer.ratio),
                                                       rownames(isotopomer.map))))
    isotopomer.map <- rbind(isotopomer.map,to.add)
    isotopomer.map <- isotopomer.map[names(isotopomer.ratio), ,drop = F]

  }

  ### iso form set
  {
    isotopomer.split <- apply(isotopomer.map, 2, function(x){paste0(x,collapse = ";")})
    isotopomer.split <-  split(seq_along(MSIPIsotopomerMap@isotopomer.defination),isotopomer.split)
    names(isotopomer.split) <- paste0("isotopomer_set_",num2str(seq_along(isotopomer.split)))
    isotopomer.set.map <- sapply(isotopomer.split ,
                               function(x){ isotopomer.map[,x[1]] },
                               USE.NAMES=F)
  }


  MSIPIsotopomerMap@solve$isotopomer.set <- isotopomer.split
  MSIPIsotopomerMap@solve$isotopomer.set.map <- isotopomer.set.map
  MSIPIsotopomerMap@solve$isotopomer.set.ratio <- isotopomer.ratio
  MSIPIsotopomerMap@solve$isotopomer.set.intensity<- isotopomer.intensity

  return(MSIPIsotopomerMap)
}

MSIPIsotopomerMap_set_solve_GLPK <- function(MSIPIsotopomerMap){

  if (isEmpty(MSIPIsotopomerMap))
    return(MSIPIsotopomerMap)

### GLPK
  {
  mat <- MSIPIsotopomerMap@solve$isotopomer.set.map
  obj <- rep(1,ncol(mat))
  dir <- rep(">=",nrow(mat))
  rhs <- MSIPIsotopomerMap@solve$isotopomer.set.ratio

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

  ### assign set to isotopomer
  {
    isotopomer.set.prob <- lp.result$solution
    isotopomer.set <- MSIPIsotopomerMap@solve$isotopomer.set
    isotopomer.prob <- mapply(x = isotopomer.set.prob ,
                           y = isotopomer.set,function(x,y){
      make_vector(x/length(y),num2str(y,10))
    },SIMPLIFY = F)
    isotopomer.prob <- unlist(isotopomer.prob)
    isotopomer.prob <- isotopomer.prob[order(names(isotopomer.prob))]
    isotopomer.prob <- unname(isotopomer.prob)
  }
  MSIPIsotopomerMap@solve$isotopomer.set.prob <- lp.result$solution
  MSIPIsotopomerMap@solve$Rglpk <- lp.result
  MSIPIsotopomerMap@isotopomer.probability <- isotopomer.prob
  return(MSIPIsotopomerMap)
}

MSIPIsotopomerMap_set_solve_QP <- function(MSIPIsotopomerMap){

  if (isEmpty(MSIPIsotopomerMap))
    return(MSIPIsotopomerMap)

  ###  Quadratic Programming
  {

    # Load necessary library

    # Example data (I_i: observed isotopomer intensities, f_ij: contribution matrix)
    I <- MSIPIsotopomerMap@solve$isotopomer.set.ratio # 3 fragments
    f <- MSIPIsotopomerMap@solve$isotopomer.set.map # 4 atoms
    message(nrow(f),"-",ncol(f))
    # Define weights based on fragment reliability (higher intensity = higher weight)
    int <- MSIPIsotopomerMap@solve$isotopomer.set.intensity # Example weights
    weights <- .intensity_weight(int)
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


  ### assign set to isotopomer
  {
    isotopomer.set.prob <- p_estimated
    isotopomer.set <- MSIPIsotopomerMap@solve$isotopomer.set
    isotopomer.prob <- mapply(x = isotopomer.set.prob ,
                           y = isotopomer.set,function(x,y){
                             make_vector(x/length(y),num2str(y,10))
                           },SIMPLIFY = F)
    isotopomer.prob <- unlist(isotopomer.prob)
    isotopomer.prob <- isotopomer.prob[order(names(isotopomer.prob))]
    isotopomer.prob <- unname(isotopomer.prob)
  }
  MSIPIsotopomerMap@solve$isotopomer.set.prob <- isotopomer.set.prob
  MSIPIsotopomerMap@solve$solve.QP <- qp_result
  MSIPIsotopomerMap@isotopomer.probability <- isotopomer.prob
  return(MSIPIsotopomerMap)
}


.intensity_weight <- function(intensity){

  weight <- log10(intensity)
  weight <- 1/(1+exp(-weight*5+20))
  return(weight)
}


get_atom_prob_from_MSIPIsotopomerMap <- function(MSIPIsotopomerMap){

  if (isEmpty(MSIPIsotopomerMap))
    return(NULL)


  atom <- unique(unlist(MSIPIsotopomerMap@isotopomer.defination))
  atom.prob <- make_vector(0,atom)
  iso_count <- length(MSIPIsotopomerMap@isotopomer.defination[[1]])
  atom.prob.m <- sapply(seq_along(MSIPIsotopomerMap@isotopomer.probability),
                    function(idx){

    x <- MSIPIsotopomerMap@isotopomer.probability[idx]
    x.atom <- MSIPIsotopomerMap@isotopomer.defination[[idx]]
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



MSIPCore_vis_isotopomer <- function(MSIPCoreData,cfmd,id,...){

  ig <- get_cfm_data_sdf_igraph(cfmd)
  atom.to.show <- MSIPCoreData@solve$MSIPIsotopomerMap@isotopomer.defination[[id]]
  message(names(MSIPCoreData@solve$MSIPIsotopomerMap@isotopomer.defination)[id])
  vis_sdf_igraph(ig,highlight = atom.to.show,...)

}


plotly_MSIPCore_pred_nature_prob  <- function(MSIPCoreData){

  im <- MSIPCoreData@solve$MSIPIsotopomerMap

  df <- data.frame(
    natural.prob =lengths(im@solve$isotopomer.set)/length(im@isotopomer.defination),
    predict.prob = im@solve$isotopomer.set.prob
  )


  p <- plot_ly(df)%>%
    add_markers(x = ~natural.prob,y = ~predict.prob,size = I(100))%>%
    add_lines(x = c(0,1),y= c(0,1),color = I("grey"))%>%
    layout(
           plot_bgcolor = 'rgba(0,0,0,0)',
           paper_bgcolor = 'rgba(0,0,0,0)',
           margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0),
           showlegend = FALSE)%>%
    config(displayModeBar = FALSE)

  #open_visNet(p)
  p




}

MSIPCore_merge <- function(MSIPCoreData1,
                           MSIPCoreData2,
                           suffix1 = "Positive",
                           suffix2 = "Negative"){


  if (!(is.null(MSIPCoreData1)|isEmpty(MSIPCoreData1))) {
    MSIPCoreData1 <- MSIPCore_FG_suffix(MSIPCoreData1,suffix1)
  }
  if (!(is.null(MSIPCoreData2)|isEmpty(MSIPCoreData2))) {
    MSIPCoreData2 <- MSIPCore_FG_suffix(MSIPCoreData2,suffix2)
  }


  if (is.null(MSIPCoreData1)|isEmpty(MSIPCoreData1)) {
    if (is.null(MSIPCoreData2)|isEmpty(MSIPCoreData2)) {
      return(NULL)
    }else{
      return(MSIPCoreData2)
    }
  }else{
    if (is.null(MSIPCoreData2)|isEmpty(MSIPCoreData2)) {
      return(MSIPCoreData1)
    }
    MSIPCoreData <- MSIPCoreData1
  }


  ### Spectra_data
  {


    MSIPCoreData@Spectra_data <-
      bind_rows(MSIPCoreData1@Spectra_data,
            MSIPCoreData2@Spectra_data)

  }


  ### FG_map
  {
    MSIPCoreData@FG_map@fragment.intensity
    ### fragment.atom.matrix
    {

      MSIPCoreData@FG_map@fragment.atom.matrix <-
        rbind(MSIPCoreData1@FG_map@fragment.atom.matrix,
            MSIPCoreData2@FG_map@fragment.atom.matrix)
    }

    ### fragment.ratio.matrix
    {

      MSIPCoreData@FG_map@fragment.ratio.matrix <-
        rbind(MSIPCoreData1@FG_map@fragment.ratio.matrix,
              MSIPCoreData2@FG_map@fragment.ratio.matrix)
    }

    ### fragment.intensity
    {

      MSIPCoreData@FG_map@fragment.intensity <-
        c(MSIPCoreData1@FG_map@fragment.intensity,
              MSIPCoreData2@FG_map@fragment.intensity)
    }

    ### fragment.include
    {

      MSIPCoreData@FG_map@fragment.include <-
        c(MSIPCoreData1@FG_map@fragment.include,
          MSIPCoreData2@FG_map@fragment.include)
    }

  }

  MSIPCoreData@solve <- list()

  return(MSIPCoreData)

}

MSIPCore_FG_suffix <- function(MSIPCoreData,
                               suffix = "suffix"){
  ### Spectra_data
  {

    MSIPCoreData@Spectra_data$fragment_group <-
      sapply(MSIPCoreData@Spectra_data$fragment_group,function(x){
        if(is.na(x)){
          return(NA)
        }else
          paste0(x,"_",suffix)
      })%>%as.character()


  }


  ### FG_map
  {
    ### fragment.atom.matrix
    {
      rownames(MSIPCoreData@FG_map@fragment.atom.matrix) %<>%
        paste0(.,"_",suffix)

    }

    ### fragment.ratio.matrix
    {
      rownames(MSIPCoreData@FG_map@fragment.ratio.matrix) %<>%
        paste0(.,"_",suffix)

    }

    ### fragment.intensity
    {
      names(MSIPCoreData@FG_map@fragment.intensity) %<>%
        paste0(.,"_",suffix)
    }

    ### fragment.include
    {
      names(MSIPCoreData@FG_map@fragment.include) %<>%
        paste0(.,"_",suffix)

    }

  }

  MSIPCoreData@solve <- list()

  return(MSIPCoreData)


}


