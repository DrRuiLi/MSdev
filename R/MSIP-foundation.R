get_frag_group_map <- function(sp.frag.data,cfmd,iso.count,atom_prob = F){

  ### frag group to label fraction
  {
    sp.frag.data <- sp.frag.data%>%
      dplyr::filter(sp.id== "combined_sp")
    if (nrow(sp.frag.data)==0) return(NA)
    fg.idx <- split(1:nrow(sp.frag.data),sp.frag.data$fragment_group)
    frag.iso.matrix <- matrix(
      nrow = length(fg.idx),ncol = iso.count+1,
      dimnames = list(names(fg.idx),paste0("M",0:iso.count)))
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
      to.add <- setdiff(paste0("M",0:iso.count),colnames(x.int))
      x.int <- cbind(matrix(0,nrow(x.int),length(to.add),
                            dimnames = list(NULL,to.add)),x.int)
      x.int <- x.int[,paste0("M",0:iso.count),drop =F]
      x.int[is.na(x.int)] <- 0
      x.weight <- rowSums(x.int)
      x.int <- t(apply(x.int,1,function(z) z/sum(z)))
      x.int.weighted <- apply(x.int,2,weighted.mean,w = x.weight)
      frag.iso.matrix[i.fg,] <- x.int.weighted
      frag.int.sum[i.fg] <- sum(x.weight)
    }
    names(frag.int.sum) <-names(fg.idx)
  }


  ### frag group to C atom prob
  {
    c_ele <- get_atom_from_igraph(get_cfm_data_sdf_igraph(cfmd),"C")
    frag.c.matrix <- matrix(ncol = length(c_ele),
                            nrow = length(fg.idx),
                            dimnames = list(names(fg.idx),
                                            c_ele))
    for (i.fg in seq_along(fg.idx)) {

      this.frag.group <- names(fg.idx)[i.fg]
      this.frags <- cfmd@fragment_define[cfmd@fragment_define$fragment_group==this.frag.group,]
      this.frag.ratio <-frag.iso.matrix[i.fg,]
      this.frag.atom <- get_cfm_data_fg_atom_map(cfmd,this.frag.group)
      this.frag.c <- this.frag.atom[c_ele]
      #this.iso.expectation <- sum(str_extract_num(names(this.frag.ratio))*this.frag.ratio)
      #this.frag.c <- this.frag.c*this.iso.expectation/sum(this.frag.c)
      #this.frag.c <- this.frag.c[this.frag.c!=0]
      frag.c.matrix[this.frag.group,names(this.frag.c)] <- this.frag.c
    }

  }


  ### return
  {
    if (atom_prob) {

    }else{
      frag.c.matrix.max <- apply(frag.c.matrix,1,function(x){
        x.sum <- sum(x)
        x.new <- rep(0,length(x))
        x.new[tail(order(x),x.sum)] <- 1
        x.new
      })%>%t
      colnames(frag.c.matrix.max) <- colnames(frag.c.matrix)
      frag.c.matrix <- frag.c.matrix.max
    }

  }

  x <- list(frag.c.matrix=frag.c.matrix,
            frag.iso.matrix=frag.iso.matrix,
            frag.int = frag.int.sum)
  return(x)

}


get_iso_form_map_old <- function(fg.map,atom_prob = T ){

  #load("temp.rda")
  #iso.count <- 5
  frag.c.matrix <- fg.map$frag.c.matrix
  frag.iso.matrix <-fg.map$frag.iso.matrix
  frag.max.iso <- ncol(frag.iso.matrix)-1

  ### all possible iso form
  iso.form <- combn(colnames(frag.c.matrix),frag.max.iso,simplify = F)
  #iso.form.iso.count.matrix <- matrix(NA,nrow(frag.c.matrix),length(iso.form))
  #iso.form.prob.matrix <- iso.form.iso.count.matrix


  ### iso form map to iso ratio
  if(!atom_prob){
    iso.form.map <- matrix(nrow = 0,ncol = length(iso.form))
    iso.form.ratio <- c()
    for (i.fg in 1:nrow(frag.c.matrix)) {

      this.c <- frag.c.matrix[i.fg,]
      this.iso <- frag.iso.matrix[i.fg,]
      this.possible.c <- names(this.c[this.c!=0])
      this.iso.count <- sapply(iso.form,function(x){sum(x %in% this.possible.c)})
      #iso.form.iso.count.matrix[i.fg,] <- this.iso.count
      #iso.form.prob.matrix[i.fg,] <- this.iso[paste0("M",this.iso.count)]

      this.iso.form.map <- lapply(0:(ncol(frag.iso.matrix)-1),function(x){
        z <-rep(0,length(this.iso.count))
        z[this.iso.count==x] <- 1
        z
      })%>%do.call("rbind",.)
      rownames(this.iso.form.map) <-paste0(rownames(frag.c.matrix)[i.fg],
                                           "_M",
                                           0:(ncol(frag.iso.matrix)-1))
      iso.form.map <- rbind(iso.form.map,this.iso.form.map)
      names(this.iso) <-paste0(rownames(frag.c.matrix)[i.fg],
                               "_M",
                               0:(ncol(frag.iso.matrix)-1))
      iso.form.ratio <-c(iso.form.ratio,this.iso)
    }
  }


  ### iso form map to iso ratio with atom prob
  if(atom_prob){

    iso.form.maps <- lapply(seq_along(iso.form),
                            function(if.id){
                              lapply(rownames(frag.c.matrix),function(fg.id){

                                get_iso_prob(frag.c.matrix[fg.id,],
                                             iso.form[[ if.id ]])
                              })->mp
                              names(mp) <- rownames(frag.c.matrix)
                              unlist(mp)
                            })
    iso.form.map <- t(do.call(bind_rows,iso.form.maps))
    iso.form.map <- iso.form.map[order(rownames(iso.form.map)),]
    iso.form.map[is.na(iso.form.map)] <- 0
    rownames(iso.form.map) <- sub(pattern = ".",x = rownames(iso.form.map),
                                  replacement = "_",fixed = T)
    iso.form.ratio <- lapply(rownames(frag.iso.matrix),function(fg.id){
      z <- frag.iso.matrix[fg.id,]
      names(z) <- paste0(fg.id,"_",names(z))
      return(z)
    })
    iso.form.ratio <- unlist(iso.form.ratio)
    iso.form.ratio <- iso.form.ratio[rownames(iso.form.map)]
  }


  ### iso form set
  {
    iso.form.split <- apply(iso.form.map, 2, function(x){paste0(x,collapse = ";")})
    iso.form.split <-  split(seq_along(iso.form),iso.form.split)
    names(iso.form.split) <- paste0("iso_form_set_",num2str(seq_along(iso.form.split)))
    iso.form.set.map <- sapply(iso.form.split ,
                                function(x){ iso.form.map[,x[1]] },
                                USE.NAMES=F)
    iso.form.set <- lapply( iso.form.split,
                            function(x){ iso.form[x] })
  }



  x <- list(iso.form.set.map = iso.form.set.map,
            iso.form.set=iso.form.set,
            iso.form.ratio=iso.form.ratio,
            iso.form.int=rep(fg.map$frag.int,each = ncol(frag.iso.matrix))
  )
  return(x)




}

get_iso_form_map <- function(fg.map,max_combn = 500000){


  ###
  {
    if (all(is.na(fg.map))) return(NA)
  }

  ### required info
  {

    frag.c.matrix <- fg.map$frag.c.matrix
    frag.iso.matrix <-fg.map$frag.iso.matrix
    frag.max.iso <- ncol(frag.iso.matrix)-1

  }
  ### all possible iso form
  {
    if.combn <- choose(ncol(frag.c.matrix),frag.max.iso)
    if (if.combn > max_combn) return(NA)
    iso.form <- combn(colnames(frag.c.matrix),frag.max.iso,simplify = F)
    names(iso.form) <- paste0("iso_form_",num2str(1:length(iso.form)))
  }
  ### iso form map to iso ratio
  {
    iso.form.maps <- lapply(seq_along(iso.form),
                            function(if.id){
                              lapply(rownames(frag.c.matrix),function(fg.id){

                                get_iso_prob(frag.c.matrix[fg.id,],
                                             iso.form[[ if.id ]])
                              })->mp
                              names(mp) <- rownames(frag.c.matrix)
                              unlist(mp)
                            })
    iso.form.map <- t(do.call(bind_rows,iso.form.maps))
    iso.form.map <- iso.form.map[order(rownames(iso.form.map)),,drop = F]
    iso.form.map[is.na(iso.form.map)] <- 0
    rownames(iso.form.map) <- sub(pattern = ".",x = rownames(iso.form.map),
                                  replacement = "_",fixed = T)
    colnames(iso.form.map) <- names(iso.form)
  }


  x <- list(iso.form = iso.form,
            iso.form.map = iso.form.map)
  return(x)


}

get_iso_form_set_map <- function(if.map,fg.map){



  ###
  {
    if (all(is.na(if.map)|all(is.na(fg.map)))) return(NA)

  }


  ### iso form ratio
  {
    frag.iso.matrix <- fg.map$frag.iso.matrix
    iso.form.ratio <- lapply(rownames(frag.iso.matrix),function(fg.id){
      z <- frag.iso.matrix[fg.id,]
      names(z) <- paste0(fg.id,"_",names(z))
      return(z)
    })
    iso.form.ratio <- unlist(iso.form.ratio)

    iso.form.map <-if.map$iso.form.map
    to.add <- matrix(0,
                     nrow = length(setdiff(names(iso.form.ratio),rownames(iso.form.map))),
                     ncol = ncol(iso.form.map),
                     dimnames = list(rowname = setdiff(names(iso.form.ratio),rownames(iso.form.map))))
    iso.form.map <- rbind(iso.form.map,to.add)
    iso.form.map <- iso.form.map[names(iso.form.ratio), ,drop = F]
    #iso.form.ratio <- iso.form.ratio[rownames(iso.form.set.map)]
    #iso.form.int <- rep(fg.map$frag.int,each = ncol(frag.iso.matrix))
  }

  ### iso form set
  {
    iso.form.split <- apply(iso.form.map, 2, function(x){paste0(x,collapse = ";")})
    iso.form.split <-  split(seq_along(if.map$iso.form),iso.form.split)
    names(iso.form.split) <- paste0("iso_form_set_",num2str(seq_along(iso.form.split)))
    iso.form.set.map <- sapply(iso.form.split ,
                               function(x){ iso.form.map[,x[1]] },
                               USE.NAMES=F)
  }

  if.map$iso.form.set <- iso.form.split
  if.map$iso.form.set.map <- iso.form.set.map
  if.map$iso.form.ratio <- iso.form.ratio
  #if.map$iso.form.int <- iso.form.int

  return(if.map)
}

get_iso_prob <- function(fc, ifc){

  ifc.prob <- fc[ifc]
  #ifc.prob[2:4] <- c(0.2,0.4,0.6)
  cv <- ifc.prob[which(ifc.prob>0&ifc.prob<1) ]
  cs <- ifc.prob[which(ifc.prob==0|ifc.prob==1) ]
  if (!length(cv)) {
    iso.p <- 1
    names(iso.p) <- paste0("M",sum(cs))
    return(iso.p)
  }

  cvpm <- matrix(rep(cv,2^length(cv)),ncol = length(cv),byrow = T)
  cvm <- expand.grid( rep(list(0:1),length(cv)))
  cvpm <- 1-cvpm-cvm
  cvpm[cvpm<0] <- -cvpm[cvpm<0]
  cv.p <- apply(cvpm,1,prod)
  cv.n <- apply(cvm,1,sum)
  iso.p <- sapply(split(cv.p,cv.n),sum)
  names(iso.p)<-paste0("M",as.numeric(names(iso.p))+sum(cs))
  iso.p
}


merge_frag_group_map <- function(fg.map){



  if (all(is.na(fg.map))) return(NA)
  if (nrow(fg.map$frag.c.matrix)<=1)  return(fg.map)

  ### merge duplicate
  {
    frag.c.matrix <- fg.map$frag.c.matrix
    frag.iso.matrix <- fg.map$frag.iso.matrix
    frag.int.sum <- fg.map$frag.int
    z <- frag.c.matrix
    #z[z>0] <- 1
    ### duplicated
    frag.c.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.c.matrix))
    frag.iso.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.iso.matrix))
    frag.int.sum1 <- c()
    z.split <- split(1:nrow(z),apply(z,1,paste0,collapse = ";"))
    for (i.z in seq_along(z.split)) {
      idx <- z.split[[i.z]]
      this.frag.c <- apply(frag.c.matrix[idx,,drop =F],2,
                           mean,weight = frag.int.sum[idx])
      this.frag.iso <- apply(frag.iso.matrix[idx,,drop =F],2,
                             mean,weight = frag.int.sum[idx])
      frag.c.matrix1 <- rbind(frag.c.matrix1,this.frag.c)
      frag.iso.matrix1 <- rbind(frag.iso.matrix1,this.frag.iso)
      frag.int.sum1 <- c(frag.int.sum1, max(frag.int.sum[idx]) )
      rn <- rownames(frag.iso.matrix)[idx][which.max(frag.int.sum[idx])]
      rownames(frag.c.matrix1)[i.z] <- rownames(frag.iso.matrix1)[i.z]  <- rn
    }
    frag.c.matrix <- frag.c.matrix1
    frag.iso.matrix <- frag.iso.matrix1
    frag.int.sum <- frag.int.sum1
    names(frag.int.sum) <- rownames(frag.c.matrix)

  }


  ### complementary
  {
    ### complementary
    z <- frag.c.matrix
    #z[z>0] <- 1
    z.comple <- apply(z,1,function(x){
      z1 <- t(t(z)+x)
      apply(z1,1,function(xx){all(xx==1)})

    })
    z.comple <- which(z.comple,arr.ind = T)
    if (any(z.comple)) {
      z.split <- lapply(1:nrow(z),function(x){
        x.c <- z.comple[z.comple[,1] == x,2]
        sort( c(x,x.c))
      })
      z.split <- unique(z.split)
      frag.c.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.c.matrix))
      frag.iso.matrix1 <- matrix(nrow = 0,ncol=ncol(frag.iso.matrix))
      frag.int.sum1 <- c()
      for (i.z in seq_along(z.split)) {
        idx <- z.split[[i.z]]
        this.frag.c <- frag.c.matrix[idx[1],,drop=F]
        this.frag.iso <- frag.iso.matrix[idx,,drop =F]
        this.frag.iso[setdiff(nrow(this.frag.iso),1),] <- this.frag.iso[setdiff(nrow(this.frag.iso),1),ncol(this.frag.iso):1]
        this.frag.iso <- apply(this.frag.iso,2,
                               weighted.mean,w = frag.int.sum[idx])
        frag.c.matrix1 <- rbind(frag.c.matrix1,this.frag.c)
        frag.iso.matrix1 <- rbind(frag.iso.matrix1,this.frag.iso)
        frag.int.sum1 <- c(frag.int.sum1, max(frag.int.sum[idx]) )
        rn <- rownames(frag.iso.matrix)[idx][which.max(frag.int.sum[idx])]
        rownames(frag.c.matrix1)[i.z] <- rownames(frag.iso.matrix1)[i.z]  <- rn
      }
      frag.c.matrix <- frag.c.matrix1
      frag.iso.matrix <- frag.iso.matrix1
      frag.int.sum <- frag.int.sum1
      names(frag.int.sum) <- rownames(frag.c.matrix)
    }

  }


  ### return
  {
    frag.c.matrix <- frag.c.matrix[order(rownames(frag.c.matrix)),,drop = F]
    frag.iso.matrix <- frag.iso.matrix[order(rownames(frag.iso.matrix)),,drop = F]
    frag.int.sum <- frag.int.sum[order(names(frag.int.sum))]

    fg.map$frag.c.matrix <- frag.c.matrix
    fg.map$frag.iso.matrix <- frag.iso.matrix
    fg.map$frag.int <- frag.int.sum
    return(fg.map)
  }


}

get_iso_form_prob_GLPK <- function(iso.form.map){

  {
    if (all(is.na(iso.form.map))) return(NA)
  }

  mat <- iso.form.map$iso.form.set.map
  obj <- rep(1,ncol(mat))
  # <- sample(c(0,1),ncol(mat),replace = T)
  # <- rnorm(ncol(mat))
  #obj[1] <- 2

  dir <- rep("==",nrow(mat))
  rhs <- iso.form.map$iso.form.ratio

  lp.result <- Rglpk::Rglpk_solve_LP(obj = obj,mat = mat,
                      canonicalize_status = T,
                      bounds = list(lower = list(ind = 1:ncol(mat),
                                                 val = rep(0,ncol(mat))),
                                    upper = list(ind = 1:ncol(mat),
                                                 val = rep(1,ncol(mat)))),
                      types = rep("C",ncol(mat)),
                      dir = dir,rhs = rhs,max = F)
  lp.result$status
  iso.form.map$iso.form.prob <- lp.result$solution
  iso.form.map$Rglpk <- lp.result
  return(iso.form.map)
}



.get_isotopologues_label_fraction <- function(sp.iso,
                                              cfmd,
                                              ppm = 10,
                                              iso.count,
                                              natural.ratio){

  sp.iso <- Spectra_filter_noise(sp.iso)
  #sp.iso <- combineSpectra_groupby_ce(sp.back,
  #                                    minProp = 0.5,plot = T,
  #                                    ppm = 20)
  sp.frag.data <- CFM_annotate_isotopologues(sp.iso,
                                             cfmd  = cfmd,
                                             ppm = ppm,
                                             iso.count = iso.count)
  sp.frag.data <- CFM_spectra_data_int_weight(sp.frag.data,iso.count)
  fg.map <- get_frag_group_map(sp.frag.data,cfmd,iso.count = iso.count)
  if.map <- get_iso_form_map(fg.map)
  ### adjust nature
  sp.frag.data <- CFM_spectra_data_remove_natural(sp.frag.data,natural.ratio,if.map)
  fg.map <- get_frag_group_map(sp.frag.data,cfmd,iso.count = iso.count)
  ### merge fragment
  fgn.map <- merge_frag_group_map(fg.map)

  ### calc by if-set
  if.map <- get_iso_form_set_map(if.map ,fgn.map)
  if.map <- get_iso_form_prob_GLPK(if.map)
  c.prob <- get_iso_from_C_prob(if.map, cfmd,iso.count)

  #heatmap.fg.map(fg.map)
  #heatmap.ifs.map(if.map)
  #sum(iso.form.map$iso.form.prob)




  ### vis
  # {
  #   this.dir <- paste0("d:/temp/nad_M",i.iso)
  #   dir.create(this.dir,showWarnings = F)
  #   hm <-  heatmap.fg.map(fg.map)
  #   export::graph2png(hm,
  #                     file = paste0(this.dir,"/Frag.map.png"),
  #                     width = 10,
  #                     height = nrow(fg.map$frag.c.matrix)*0.8)
  #   hm <- heatmap.ifs.map(iso.form.map)
  #   export::graph2png(draw(hm),
  #                     file = paste0(this.dir,"/Iso.form.map.png"),
  #                     width = 10,
  #                     height = nrow(fg.map$frag.c.matrix)*1.5)
  #   p <- vis_sdf_ig_prob(cfmd@fragment_igraph[[1]],c.prob,show.label = T)
  #   saveWidget(p,file = paste0(this.dir,"/Atom.prob.html"))
  # }

  ### SAVE
  {
    return(list(
      sp.data = sp.frag.data,
      fg.map = fg.map,
      fgn.map = fgn.map,
      if.map = if.map,
      c.prob = c.prob
    ))


  }


}

.re_calculate_isotopologues_label_fraction <- function(msip.core.data){

  ### select include
  {
    fg.map <- msip.core.data$fg.map
    #fg.map$frag.include <- fg.map$frag.int>1e3
    if (all(is.null(fg.map$frag.include))) {

      frag.include <- rep(T,nrow(fg.map$frag.c.matrix))
      names(frag.include) <- rownames(fg.map$frag.c.matrix)
    }else{
      frag.include <- fg.map$frag.include
    }
    fg.map$frag.c.matrix <- fg.map$frag.c.matrix [frag.include,]
    fg.map$frag.iso.matrix <- fg.map$frag.iso.matrix [frag.include,]
    fg.map$frag.int <- fg.map$frag.int [frag.include]
    fg.map$frag.include <- fg.map$frag.include [frag.include]
  }

  ### merge fragment
  fgn.map <- merge_frag_group_map(fg.map)

  ### calc by if-set
  if.map <- get_iso_form_set_map(msip.core.data$if.map ,fgn.map)
  if.map <- get_iso_form_prob_GLPK(if.map)
  c.prob <- get_iso_from_C_prob(if.map)

  ### return
  {
    msip.core.data$fgn.map <- fgn.map
    msip.core.data$if.map <- if.map
    msip.core.data$c.prob <- c.prob
    return(msip.core.data)
  }

}


get_iso_from_C_prob <- function(iso.form.map){

  if (all(is.na(iso.form.map))) return(NA)

  c <- unique(unlist(iso.form.map$iso.form))
  c.prob <- rep(0,length(c))
  names(c.prob) <- c
  c.prob.m<- sapply(seq_along(iso.form.map$iso.form.prob),function(idx){

    x <- c.prob
    if (iso.form.map$iso.form.prob[idx] == 0) {

    }else{
      c.count <- table(unlist(iso.form.map$iso.form[iso.form.map$iso.form.set[[idx]]]))
      x.prob <- iso.form.map$iso.form.prob[idx]*c.count/sum(c.count)
      x[names(x.prob)] <- x.prob
    }
    return(x)
  })

  x <- apply(c.prob.m,1,sum)
  if (all(x==0)) {
    return(x)
  }
  x/sum(x)*length(iso.form.map$iso.form[[1]])
}



heatmap.ifs.map <- function(iso.form.map){

  iso.form.set.map <- iso.form.map$iso.form.set.map
  iso.form.prob <- iso.form.map$iso.form.prob
  if ("iso.form.prob" %in% names(iso.form.map)) {
    top.anno <- HeatmapAnnotation(
      ifp = iso.form.prob,
      col = list(ifp = colramp(c(0,0.00000001,max(iso.form.prob)),
                               c("grey","white","#0095D4"))),
      annotation_label  = list(ifp = "iso form probability"),
      show_annotation_name = F,
      which = "c"
    )
  }

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
            Ratio = iso.form.map$iso.form.ratio,
            col = list(Ratio = colramp(c(0,0.0000001,0.5,1),
                                       c("grey","white","#F7844F","#B20C26"))),
            which = "row",
            show_annotation_name = F,
            width  = unit(20,"inch")),
          show_heatmap_legend = F,
          col = colramp(colors = c("white","white","black")))
}


get_ele_uniso <- function(iso_ele = "[13]C"){

  sub(pattern = "\\[.+\\]",
      x = iso_ele,
      replacement = "")

}

trans_iso_ele <- function(iso_ele = "[13]C"){

  if (grepl(pattern = "\\[",x = iso_ele)) {
    x <- paste0(str_extract(iso_ele,"[:alpha:]+"),
           str_extract(iso_ele,"[:digit:]+"))
  }else{
    x <- paste0("[",str_extract(iso_ele,"[:digit:]+"),"]",
                str_extract(iso_ele,"[:alpha:]+") )
  }
  return(x)

}


get_formula_ele_count <- function(formula,ele = "C"){

  .f<-function(formula){
    if (is.na(formula)) {
      return(NA)
    }
    atom.count <- MSCC:::chemform_parse(formula)
    if (!ele%in% colnames(atom.count)) {
      max.atom <- 0
    }else
      max.atom <- atom.count[,ele]
  }
  max.atom <- unname(sapply(formula,.f))
  return(max.atom)

}
