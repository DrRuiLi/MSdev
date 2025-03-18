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
                                     natural.ratio){

  ### check status
  {
    if (is.null(MSIPCoreData@solve$MSIPIsotopomerMap))
      MSIPCoreData@solve$MSIPIsotopomerMap <- get_MSIPIsotopomerMap(MSIPCoreData)


    natural_corrected <- MSIPCoreData@solve$natural_corrected
    natural_corrected <- ifelse(is.null(natural_corrected),F,natural_corrected)
    if (natural_corrected) {
      return(MSIPCoreData)
    }
    if (natural.ratio >= 1) {
      return(MSIPCoreData)
    }
    if (isEmpty(MSIPCoreData)) {
      return(MSIPCoreData)
    }
  }


  ### correct
  {

    MSIPIsotopomerMap <- MSIPCoreData@solve$MSIPIsotopomerMap
    FSIS.prob <- MSIPIsotopomerMap@solve$isotopomer.set.prob
    FSIS.natural <- lengths(MSIPIsotopomerMap@solve$isotopomer.set)%>%
      `/`(sum(.))
    prob <- FSIS.prob-FSIS.natural*natural.ratio
    prob[prob<0] <-0
    prob <- prob/sum(prob)

    MSIPIsotopomerMap@solve$isotopomer.set.prob <- prob


    isotopomer.prob <- mapply(x = prob ,
                              y = MSIPIsotopomerMap@solve$isotopomer.set,
                              function(x,y){
                                make_vector(x/length(y),y)
                              },SIMPLIFY = F)%>%
      unname()%>%unlist()
    isotopomer.prob <- isotopomer.prob[order(names(isotopomer.prob))]
    MSIPIsotopomerMap@isotopomer.probability <- isotopomer.prob


    MSIPCoreData@solve$Atom_prob <- get_atom_prob_from_MSIPIsotopomerMap(MSIPIsotopomerMap)
    MSIPCoreData@solve$MSIPIsotopomerMap <- MSIPIsotopomerMap
    MSIPCoreData@solve$natural_corrected <- T
    MSIPCoreData@solve$natural_ratio  <- natural.ratio
  }


  return(MSIPCoreData)

}


MSIPCore_solve <- function(MSIPCoreData,
                           max_prob_map = F,
                           int_thresh = 10^3.8,
                           certainty_thresh = 0.6,
                           weight_fun = .intensity_weight,
                           re_split_isotopomers = T){

  if (isEmpty(MSIPCoreData))
    return(MSIPCoreData)

  ### set all include
  MSIPCoreData@FG_map@FG.data$include <- MSIPCoreData@FG_map@FG.data$include |T
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
  if(re_split_isotopomers)  MSIPCoreData.temp@solve <- list()

  MSIPIsotopomerMap <- get_MSIPIsotopomerMap(MSIPCoreData.temp)
  MSIPIsotopomerMap <- MSIPIsotopomerMap_set_split(MSIPIsotopomerMap,
                                       MSIPFragmentMap_reduced)
  #MSIPIsotopomerMap <- MSIPIsotopomerMap_set_solve_GLPK(MSIPIsotopomerMap)
  MSIPIsotopomerMap <- MSIPIsotopomerMap_set_solve_QP(MSIPIsotopomerMap,
                                                      weight_fun = weight_fun)
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

  ### FG - Ratio
  {
    sp.frag.data.merged <- sp.frag.data%>%
      dplyr::filter(merged)
    if (nrow(sp.frag.data.merged)==0) return(new("MSIPFragmentMap"))
    fg.idx <- split(1:nrow(sp.frag.data.merged),sp.frag.data.merged$fragment_group)
    frag.ratio.matrix <- matrix(
      nrow = length(fg.idx),ncol = iso_count+1,
      dimnames = list(names(fg.idx),paste0("M",0:iso_count)))
    for (i.fg in seq_along(fg.idx)) {


      x.df <- sp.frag.data.merged[fg.idx[[i.fg]],]
      x.ratio <- x.df%>%
        tidyr::pivot_wider(names_from ="iso_count",
                           id_cols = "sp.id",
                           values_from = "ratio",
                           values_fn = mean
        )%>%
        tibble::column_to_rownames("sp.id")%>%
        dplyr::select(dplyr::starts_with("M"))%>%
        as.matrix()
      x.ratio <- get_matrix_value_fill_with_NA(
        x.ratio,
        rownames_vec = row.names(x.ratio),
        colnames_vec =colnames(frag.ratio.matrix),drop = F)
      x.ratio[is.na(x.ratio)] <- 0
      frag.ratio.matrix[i.fg,] <- x.ratio
    }
    frag.int.sum <- sp.frag.data.merged%>%
      dplyr::distinct(fragment_group,.keep_all = T)%>%
      dplyr::pull(int_sum,name = fragment_group)
    frag.int.sum <- frag.int.sum[rownames(frag.ratio.matrix)]

  }


  ### FG - atom
  {

    target_atoms <- get_sdf_igraph_atom(get_cfm_data_sdf_igraph(cfmd),get_ele_uniso(iso_ele))
    frag.atom.matrix <- matrix(ncol = length(target_atoms),
                            nrow = length(fg.idx),
                            dimnames = list(names(fg.idx),
                                            target_atoms))
    for (i.fg in seq_along(fg.idx)) {

      this.frag.group <- names(fg.idx)[i.fg]
      #this.frags <- cfmd@fragment_define[cfmd@fragment_define$fragment_group==this.frag.group,]
      #this.frag.ratio <-frag.ratio.matrix[i.fg,]
      #this.frag.atom <- get_cfm_data_fragment_group_atom_map(cfmd,this.frag.group)
      #this.frag.c <- this.frag.atom[target_atoms]
      #this.iso.expectation <- sum(str_extract_num(names(this.frag.ratio))*this.frag.ratio)
      #this.frag.c <- this.frag.c*this.iso.expectation/sum(this.frag.c)
      #this.frag.c <- this.frag.c[this.frag.c!=0]
      this.frag.c <-cfmd@fragment_group_map[this.frag.group,]
      frag.atom.matrix[this.frag.group,names(this.frag.c)] <- this.frag.c
    }

  }

  ### FG data
  {
    fg.data <- sp.frag.data.merged%>%
      dplyr::distinct(fragment_group,.keep_all = T)%>%
      dplyr::mutate(include =T )%>%
      dplyr::select(-c("sp.id",
                       "mz",
                       "intensity",
                       "collisionEnergy",
                       "iso_count" ,
                       "ratio",
                       "merged"
                       ))%>%
      as.data.frame()
    rownames(fg.data) <- fg.data$fragment_group
    fg.data <- fg.data[rownames(frag.atom.matrix),]
  }



    fg.map <- new("MSIPFragmentMap")
    fg.map@FG.atom.matrix <- frag.atom.matrix
    fg.map@FG.ratio.matrix <- frag.ratio.matrix
    fg.map@FG.data <- fg.data

    return(fg.map)






}


get_MSIPIsotopomerMap <- function(MSIPCoreData){




  ### filter
  {
    #if (!is.null(MSIPCoreData@solve$MSIPIsotopomerMap)) {
    if(F){
      MSIPIsotopomerMap <- MSIPCoreData@solve$MSIPIsotopomerMap
      frag.ratio.matrix <- MSIPCoreData@FG_map@FG.ratio.matrix
      isotopomer.ratio <- lapply(rownames(frag.ratio.matrix),function(fg.id){
        z <- frag.ratio.matrix[fg.id,]
        names(z) <- paste0(fg.id,"_",names(z))
        return(z)
      })
      isotopomer.ratio <- unlist(isotopomer.ratio)
      isotopomer.intensity <- MSIPCoreData@FG_map@FG.data[gsub(x = names(isotopomer.ratio),
                                                               pattern = "_M.*",replacement = ""),
                                                          "int_sum"]
      names(isotopomer.intensity) <- names(isotopomer.ratio)

      #MSIPIsotopomerMap@isotopomer.map <- MSIPIsotopomerMap@isotopomer.map[names(isotopomer.ratio),,drop = F]
      #MSIPIsotopomerMap@isotopomer.ratio <- MSIPIsotopomerMap@isotopomer.ratio[names(isotopomer.ratio)]
      #MSIPIsotopomerMap@Labeled.FG.data$intensity <- isotopomer.intensity
      return(MSIPIsotopomerMap)

    }

  }

  ### required info
  {

    frag.atom.matrix <- MSIPCoreData@FG_map@FG.atom.matrix
    frag.ratio.matrix <-MSIPCoreData@FG_map@FG.ratio.matrix
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


  ### isotopomers - labeled FG
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
    isotopomer.map <- isotopomer.map[!apply(isotopomer.map,1,function(x){all(x==0)}),,drop = F]
    rownames(isotopomer.map) <- sub(pattern = ".",x = rownames(isotopomer.map),
                                  replacement = "_",fixed = T)
    colnames(isotopomer.map) <- names(isotopomer)
  }

  ### LFG data
  {

    isotopomer.ratio <- lapply(rownames(frag.ratio.matrix),function(fg.id){
      z <- frag.ratio.matrix[fg.id,]
      names(z) <- paste0(fg.id,"_",names(z))
      return(z)
    })
    isotopomer.ratio <- unlist(isotopomer.ratio)
    Labeled.FG.data <- MSIPCoreData@FG_map@FG.data[gsub(x = names(isotopomer.ratio),
                                                             pattern = "_M.*",replacement = ""),
                                                        ]
    Labeled.FG.data$labeled.FG <- names(isotopomer.ratio)
    Labeled.FG.data$ratio <- isotopomer.ratio
    rownames(Labeled.FG.data) <-  names(isotopomer.ratio)
  }




  MSIPIsotopomerMap@isotopomer.defination <-isotopomer
  MSIPIsotopomerMap@isotopomer.map <-isotopomer.map
  MSIPIsotopomerMap@Labeled.FG.data <-Labeled.FG.data
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

  frag.include <- MSIPFragmentMap@FG.data%>%
    dplyr::filter(include)%>%
    dplyr::pull(fragment_group)

  MSIPFragmentMap@FG.atom.matrix <- MSIPFragmentMap@FG.atom.matrix[frag.include,,drop = F]
  MSIPFragmentMap@FG.ratio.matrix <- MSIPFragmentMap@FG.ratio.matrix[frag.include,,drop = F]
  MSIPFragmentMap@FG.data <-  MSIPFragmentMap@FG.data[frag.include,]
  return(MSIPFragmentMap)
}

MSIPFragmentMap_filter_intensity <- function(MSIPFragmentMap,
                                            int_thresh = 1E3){

  frag.int <- MSIPFragmentMap@FG.data$int_sum
  frag.include <- frag.int > int_thresh

  frag.include <- frag.include & MSIPFragmentMap@FG.data$include
  MSIPFragmentMap@FG.data$include <-  frag.include

  return(MSIPFragmentMap)
}

MSIPFragmentMap_filter_certainty <- function(MSIPFragmentMap,
                                            certainty_thresh = 0.6 ){

  frag.certainty <- get_MSIPFragmentMap_certainty(MSIPFragmentMap)
  frag.include <- frag.certainty >= certainty_thresh

  frag.include <- frag.include & MSIPFragmentMap@FG.data$include

  MSIPFragmentMap@FG.data$include <-  frag.include
  return(MSIPFragmentMap)
}
MSIPFragmentMap_add_constraint <- function(MSIPFragmentMap){

  m <- MSIPFragmentMap@FG.atom.matrix
  MSIPFragmentMap@FG.atom.matrix <-
    rbind(m,matrix(c(0,1),2,ncol(m),dimnames = list(c("FGNULL","FGFULL"),colnames(m))))

  m <- MSIPFragmentMap@FG.ratio.matrix
  m1 <- matrix(0,2,ncol(m),dimnames = list(c("FGNULL","FGFULL"),colnames(m)))
  m1[1] <- m1[length(m1)] <- 1
  MSIPFragmentMap@FG.ratio.matrix <- rbind(m,m1)

  MSIPFragmentMap@FG.data  <- MSIPFragmentMap@FG.data %>%
    bind_rows(.,   data.frame(
      fragment_group   = c("FGNULL","FGFULL"),
      int_sum = 1e9,
      peaks_count = 9999,
      icc = 1,cos=1,
      include = T
      ))
  rownames(MSIPFragmentMap@FG.data) <- MSIPFragmentMap@FG.data$fragment_group

  return(MSIPFragmentMap)

}

get_MSIPFragmentMap_certainty <- function(MSIPFragmentMap){
  frag.matrix <- MSIPFragmentMap@FG.atom.matrix
  frag.certainty <- apply(frag.matrix,1,function(x){
    sum(x==1)/sum(x)
  })
  frag.certainty[rowSums(frag.matrix)==0] <- 0
  return(frag.certainty)
}

MSIPFragmentMap_merge_duplicate <- function(MSIPFragmentMap){


  frag.atom.matrix <- MSIPFragmentMap@FG.atom.matrix
  frag.ratio.matrix <-  MSIPFragmentMap@FG.ratio.matrix
  frag.int.sum <- MSIPFragmentMap@FG.data$int_sum
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

  frag.atom.matrix1 -> MSIPFragmentMap@FG.atom.matrix
  frag.ratio.matrix1 -> MSIPFragmentMap@FG.ratio.matrix
  frag.int.sum1 ->  MSIPFragmentMap@FG.data$int_sum
 # MSIPFragmentMap@FG.data <- MSIPFragmentMap@FG.data[]

  return(MSIPFragmentMap)
}


MSIPFragmentMap_merge_complementary <- function(MSIPFragmentMap){

  if (isEmpty(MSIPFragmentMap))
    return(MSIPFragmentMap)

  frag.atom.matrix <- MSIPFragmentMap@FG.atom.matrix
  frag.ratio.matrix <-  MSIPFragmentMap@FG.ratio.matrix
  frag.int.sum <- MSIPFragmentMap@FG.intensity

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

  frag.atom.matrix1 -> MSIPFragmentMap@FG.atom.matrix
  frag.ratio.matrix1 -> MSIPFragmentMap@FG.ratio.matrix
  frag.int.sum1 ->  MSIPFragmentMap@FG.intensity
  MSIPFragmentMap@FG.data$include <- MSIPFragmentMap@FG.data$include[names(frag.int.sum1)]
  return(MSIPFragmentMap)
}



#' heatmap_MSIPFragmentMap
#' @describeIn MSIPCore heatmap fragment map
#'
#' @param MSIPFragmentMap MSIPFragmentMap
#'
#' @return heatmap
#' @export
heatmap_MSIPFragmentMap <- function(MSIPFragmentMap,
                                    show_ratio = F){

  ### size
  {
    ppi <- 72
    mhpx <- 500
    mwpx <- 800
    nr <- nrow(MSIPFragmentMap@FG.atom.matrix)
    nc <- ncol(MSIPFragmentMap@FG.atom.matrix)+
      ncol(MSIPFragmentMap@FG.ratio.matrix)
    cell.unit <- "inches"
    length.unit <- 0.3
    max.height <- mhpx/ppi
    max.width <- mwpx/ppi-1
    if (nr*length.unit>max.height)
      length.unit <- max.height/nc
    if (nc*length.unit>max.width)
      length.unit <- max.width/nc

  }

  frag.atom.matrix <- MSIPFragmentMap@FG.atom.matrix
  frag.ratio.matrix <- MSIPFragmentMap@FG.ratio.matrix
  frag.include <- MSIPFragmentMap@FG.data$include
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
                                #column_names_gp = grid::gpar(fontsize = 8),
                                row_names_gp = gpar(col = ifelse(frag.include,"black","grey")),
                                #rect_gp =  grid::gpar(lwd=2,col = "white"),
                                cluster_rows = F)
  h1
  if(show_ratio){
    cellfun <- function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", frag.ratio.matrix[i, j]), x, y, gp = gpar(col = "black", fontsize = 10))
    }
  }else cellfun <- NULL

  lgint <- round(log10(MSIPFragmentMap@FG.data$int_sum),1)
  lgint <- c(lgint,0)
  h2 <- ComplexHeatmap::Heatmap(frag.ratio.matrix,
                                na_col  ="#999999",
                                cell_fun = cellfun,
                                # width = unit( ncol(frag.ratio.matrix)*length.unit,cell.unit),
                                # height =unit( nrow(frag.atom.matrix)*length.unit,cell.unit),
                                name = "Isotope labeled\nratio",
                                col = circlize::colorRamp2(breaks = c(0,0.5,1),
                                                           c("white","#F7844F","#B20C26")),
                                right_annotation  = rowAnnotation(
                                  intensity = ComplexHeatmap::anno_numeric(lgint,
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
                                     MSIPFragmentMap_reduced,
                                     approx_split = T,binwidth = 0.02
                                     ){


  if (isEmpty(MSIPIsotopomerMap))
    return(MSIPIsotopomerMap)



  ### Isotopomers to FSIS
  {
    isotopomer.map <- MSIPIsotopomerMap@isotopomer.map
    if (approx_split) {
      bins <- seq(0,1,binwidth)
      closest_value <- function(x, values) {
        values[which.min(abs(values - x))]
      }
      isotopomer.map <- apply(isotopomer.map,c(1,2),closest_value,values = bins)
    }


    isotopomer.split <- apply(isotopomer.map, 2, function(x){paste0(x,collapse = ";")})
    isotopomer.split <-  split(colnames(isotopomer.map),
                               isotopomer.split)
    names(isotopomer.split) <- paste0("FSIS_",num2str(seq_along(isotopomer.split)))
    isotopomer.set.map <- sapply(isotopomer.split ,
                               function(x){ isotopomer.map[,x[1]] },
                               USE.NAMES=F)
  }


  MSIPIsotopomerMap@solve$isotopomer.set <- isotopomer.split
  MSIPIsotopomerMap@solve$isotopomer.set.map <- isotopomer.set.map

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

MSIPIsotopomerMap_set_solve_QP <- function(MSIPIsotopomerMap,
                                           weight_fun = .intensity_weight){

  if (isEmpty(MSIPIsotopomerMap))
    return(MSIPIsotopomerMap)

  ###  Quadratic Programming
  {

    # Load necessary library

    # Example data (I_i: observed isotopomer intensities, f_ij: contribution matrix)
    f <- MSIPIsotopomerMap@solve$isotopomer.set.map # 4 map
    I <- MSIPIsotopomerMap@Labeled.FG.data[rownames(f),"ratio"] # 3 FISIS
    message(nrow(f),"-",ncol(f))
    # Define weights based on fragment reliability (higher intensity = higher weight)
    Labeled.FG.data <- MSIPIsotopomerMap@Labeled.FG.data[rownames(f),]  # get weights
    cer <- apply(f,1,function(x){
      sum(x==1)/sum(x)
    })
    cer[is.nan(cer)] <- 1
    weights <- weight_fun(Labeled.FG.data$int_sum*I)*(Labeled.FG.data$cos-0.5)
    #weights <- make_vector(1,I)
    weight.df <- data.frame(
      frag.total = Labeled.FG.data$int_sum,
      peak.int = Labeled.FG.data$int_sum*I,
      ratio = I,
      w = weights
    )
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
      int = log10(Labeled.FG.data$int_sum)
    )
    }


  ### assign set to isotopomer
  {

    p_estimated[p_estimated<0] <-0
    isotopomer.set.prob <- p_estimated/sum(p_estimated)
    names(isotopomer.set.prob) <- colnames(f)
    isotopomer.set <- MSIPIsotopomerMap@solve$isotopomer.set
    isotopomer.prob <- mapply(x = isotopomer.set.prob ,
                           y = isotopomer.set[names(isotopomer.set.prob)],
                           function(x,y){
                             make_vector(x/length(y),y)
                           },SIMPLIFY = F)%>%
      unname()%>%unlist()
    isotopomer.prob <- isotopomer.prob[order(names(isotopomer.prob))]
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

.intensity_weight_astral <- function(intensity){
  weight <- log10(intensity)
  #x <- ifelse(weight>4,weight-4,0)
  weight <- 1/(1+exp(-weight*10+25))
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
  r2 <- summary(lm(natural.prob~predict.prob,data = df))$r.squared
  rmse <-  sqrt(mean((df$natural.prob - df$predict.prob) ^ 2))
  p <- plot_ly(df)%>%
    add_markers(x = ~natural.prob,y = ~predict.prob,size = I(100))%>%
    add_lines(x = c(0,1),y= c(0,1),color = I("grey"))%>%
    add_text(x = 0.3,y=0.9,text = paste0("R2 = ",str_digit(r2,4),"\n",
                                         "RMSE = ",str_digit(rmse,4)
                                         ))%>%
    layout(
           plot_bgcolor = 'rgba(0,0,0,0)',
           paper_bgcolor = 'rgba(0,0,0,0)',
           margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0),
           showlegend = FALSE)%>%
    plotly::config(displayModeBar = FALSE)

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

    ### FG.atom.matrix
    {

      MSIPCoreData@FG_map@FG.atom.matrix <-
        rbind(MSIPCoreData1@FG_map@FG.atom.matrix,
            MSIPCoreData2@FG_map@FG.atom.matrix)
    }

    ### FG.ratio.matrix
    {

      MSIPCoreData@FG_map@FG.ratio.matrix <-
        rbind(MSIPCoreData1@FG_map@FG.ratio.matrix,
              MSIPCoreData2@FG_map@FG.ratio.matrix)
    }

    ### FG.intensity
    {

      MSIPCoreData@FG_map@FG.data <-
        rbind(MSIPCoreData1@FG_map@FG.data,
              MSIPCoreData2@FG_map@FG.data)
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
    ### FG.atom.matrix
    {
      rownames(MSIPCoreData@FG_map@FG.atom.matrix) %<>%
        paste0(.,"_",suffix)

    }

    ### FG.ratio.matrix
    {
      rownames(MSIPCoreData@FG_map@FG.ratio.matrix) %<>%
        paste0(.,"_",suffix)

    }

    ### FG.data
    {
      rownames(MSIPCoreData@FG_map@FG.data) %<>%
        paste0(.,"_",suffix)
      MSIPCoreData@FG_map@FG.data$fragment_group%<>%
        paste0(.,"_",suffix)
    }



  }

  MSIPCoreData@solve <- list()

  return(MSIPCoreData)


}


get_MSIPCore_isotopomer_data <- function(MSIPCoreData){

  FSIS.def <- MSIPCoreData@solve$MSIPIsotopomerMap@solve$isotopomer.set
  FSIS.prob <- MSIPCoreData@solve$MSIPIsotopomerMap@solve$isotopomer.set.prob
  isotopomer.def <-  MSIPCoreData@solve$MSIPIsotopomerMap@isotopomer.defination

  data.frame(
    FSIS = names(FSIS.def)
  )
  isotopomer.iso_count <- unique(lengths(isotopomer.def))
  isotopomer.def.formate <- sapply(isotopomer.def,function(x){
    paste0("M+",isotopomer.iso_count,":",paste0(x,collapse = ","))
  })
  lapply(FSIS.def,function(x) data.frame(isotopomer = x))%>%
    data.table::rbindlist(idcol = "FSIS")%>%
    dplyr::mutate(prob = FSIS.prob[FSIS],
                  name = isotopomer.def.formate[isotopomer],
                  iso_count = isotopomer.iso_count
                  )

}



plot_MSIPCore_solve_weight_fun <- function(weight_fun = .intensity_weight){

  plot.data <- data.frame(
    x = seq(0,7,0.01)
  )%>%
    dplyr::mutate(int = 10^x,
                  w = weight_fun(int))
  ggplot(plot.data)+
    geom_point(aes(x = int, y = w ))+
    scale_x_log10()+
    labs(x = "Intensity", y = "Weight")+
    theme_bw()
}


vis_MSIPcore_isotopoer_set <- function(msip.core,
                                       mol.ig,
                                       is.idx = 1){

  if (is.null(is.idx)) {
    return(NULL)
  }
  if (class(is.idx)=="character") {
    is.idx <- match(is.idx,names(msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set))
  }
  is.prob <- msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set.prob[[is.idx]]
  isotopomer.idx <- msip.core@solve$MSIPIsotopomerMap@solve$isotopomer.set[[is.idx]]
  isotopomer.list <- msip.core@solve$MSIPIsotopomerMap@isotopomer.defination[isotopomer.idx]
  iso_count <- unique( lengths(isotopomer.list))
  atom.prob <- table(unlist(isotopomer.list))/length(isotopomer.list)*is.prob
  mol.ig%>%
    sdf_igraph_add_background_color(atom.prob)%>%
    sdf_igraph_add_border_color(make_vector(1,names(atom.prob)),
                                color.ramp =
                                  colramp(colors = c("grey","white","#000000")))%>%
    vis_sdf_igraph()

}



vis_MSIPcore_isotopomer <- function(msip.core,
                                       mol.ig,
                                       is.idx = 1){

  if (is.null(is.idx)) {
    return(NULL)
  }
  if (class(is.idx)=="character") {
    is.idx <- match(is.idx,names(msip.core@solve$MSIPIsotopomerMap@isotopomer.defination))
  }
  is.prob <- msip.core@solve$MSIPIsotopomerMap@isotopomer.probability[[is.idx]]
  isotopomer.idx <-is.idx
  isotopomer.list <- msip.core@solve$MSIPIsotopomerMap@isotopomer.defination[isotopomer.idx]
  iso_count <- unique( lengths(isotopomer.list))
  atom.prob <- table(unlist(isotopomer.list))/length(isotopomer.list)*is.prob
  mol.ig%>%
    sdf_igraph_add_background_color(atom.prob)%>%
    sdf_igraph_add_border_color(make_vector(1,names(atom.prob)),
                                color.ramp =
                                  colramp(colors = c("grey","white","#1485EE")))%>%
    vis_sdf_igraph()

}
