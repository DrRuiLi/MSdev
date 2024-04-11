get_frag_group_map <- function(sp.frag.data,iso.count){

  ### frag group to label fraction
  {
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
      x.int.weighted <- apply(x.int,2,weighted.mean,m = x.weight)
      frag.iso.matrix[i.fg,] <- x.int.weighted
      frag.int.sum[i.fg] <- sum(x.weight)
    }
    names(frag.int.sum) <-names(fg.idx)
  }


  ### frag group to C atom prob
  {
    frag.c.matrix <- matrix(ncol = length(c_ele),
                            nrow = length(fg.idx),
                            dimnames = list(names(fg.idx),
                                            names(c_ele)))
    for (i.fg in seq_along(fg.idx)) {

      this.frag.group <- names(fg.idx)[i.fg]
      this.frags <- cfmd@fragment_define[cfmd@fragment_define$fragment_group==this.frag.group,]
      this.frag.ratio <-frag.iso.matrix[i.fg,]
      this.frag.atom <- get_cfm_data_fg_atom_map(cfmd,this.frag.group)
      this.frag.c <- this.frag.atom[names(c_ele)]
      #this.iso.expectation <- sum(str_extract_num(names(this.frag.ratio))*this.frag.ratio)
      #this.frag.c <- this.frag.c*this.iso.expectation/sum(this.frag.c)
      #this.frag.c <- this.frag.c[this.frag.c!=0]
      frag.c.matrix[this.frag.group,names(this.frag.c)] <- this.frag.c
    }

  }

  x <- list(frag.c.matrix=frag.c.matrix,
            frag.iso.matrix=frag.iso.matrix,
            frag.int = frag.int.sum)
  return(x)

}


get_iso_form_map <- function(fg.map,atom_prob = T,vis = F){

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

      this.iso.form.map <- sapply(0:(ncol(frag.iso.matrix)-1),function(x){
        z <-rep(0,length(this.iso.count))
        z[this.iso.count==x] <- 1
        z
      })%>%t
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

  }


  ### iso form set
  {
    iso.form.split <- apply(iso.form.map, 2, function(x){paste0(x,collapse = ";")})
    iso.form.split <-  split(seq_along(iso.form),iso.form.split)
    names(iso.form.split) <- paste0("iso_form_set_",num2str(seq_along(iso.form.split)))
    iso.form.set.map <- sapply(iso.form.split ,
                                function(x){ iso.form.map[,x[1]] },
                                USE.NAMES=F)
    iso.form.set <- sapply( iso.form.split,
                            function(x){ iso.form[x] })
  }

  ### vis
  if(vis){
    #dim(iso.form.set.map)
    Heatmap(iso.form.set.map,

            border = T,
            border_gp = gpar(col = "#808080"),
            row_split = rep(1:nrow(frag.c.matrix),each = ncol(frag.iso.matrix)),
            row_title = NULL,


            cluster_rows = F,
            show_column_names = F,
            show_column_dend = T,
            show_row_names = T,
            row_names_side  = "l",
            column_title = "Iso form set",
            column_title_gp = gpar(fontsize = 30),
            cluster_columns = T,
            left_annotation = HeatmapAnnotation(
              Ratio = (iso.form.ratio),
              which = "row",
              show_annotation_name = F,
              width  = unit(20,"inch"),
              col = list(Ratio = colramp())),
            show_heatmap_legend = F,
            col = colramp(colors = c("white","white","black")))->p
    #p
    open_plot_win(p,height = nrow(iso.form.set.map)*0.2+1,
                  width = nrow(iso.form.set.map)*0.2*1.5)

  }


  x <- list(iso.form.set.map = iso.form.set.map,
            iso.form.set=iso.form.set,
            iso.form.ratio=iso.form.ratio,
            iso.form.int=rep(fg.map$frag.int,each = ncol(frag.iso.matrix))
  )
  return(x)




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
    z.split <- sapply(1:nrow(z),function(x){
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
    fg.map$frag.c.matrix <- frag.c.matrix1
    fg.map$frag.iso.matrix <- frag.iso.matrix1
    fg.map$frag.int <- frag.int.sum1

  }

  return(fg.map)
}


