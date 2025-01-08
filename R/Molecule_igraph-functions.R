get_Molecule_igraph_from_sdf <- function(sdf = sdfsample) {


  if (class(sdf) == "SDF")
    Molecule_igraphs <- sdf_to_Molecule_igraph(sdf)


  if (class(sdf) == "SDFset") {
    Molecule_igraphs <- list()
    sdf.valid <- validSDF(sdf)
    for (i in 1:length(sdf)) {
      Molecule_igraphs[[i]] <- sdf_to_Molecule_igraph(sdf[[i]])
    }
    names(Molecule_igraphs) <- cid(sdf)
  }

  return(Molecule_igraphs)

}


sdf_to_Molecule_igraph <- function(sdf) {


  atom.block <- atomblock(sdf) %>%
    as.data.frame() %>%
    rownames_to_column("id") %>%
    dplyr::mutate(name = id,
                  bonds(sdf),
                  element = atom,
                  x = C1,
                  y = C2,
                  .before = id)

  bond.block  <- bondblock(sdf) %>%
    as.data.frame() %>%
    dplyr::mutate(from = atom.block$id[C1],
                  to = atom.block$id[C2],
                  bond_type = C3,
                  .before = C1)

  sdf.igraph <- igraph::graph_from_data_frame(bond.block, vertices = atom.block)

  isotopomer.df <-  make_vector(atom.block$element,atom.block$id)%>%
    as.list()%>%
    as.data.frame()%>%
    dplyr::mutate(isotopomer = "natural",
                  isotopologue = "M0",
                  label = "natural",
                  abundance = 1,
                  .before = everything()
    )
  rownames(isotopomer.df) <- isotopomer.df$isotopomer
  Molecule_igraph <- new("Molecule_igraph",
                         sdf = sdf, igraph = sdf.igraph,isotopomer = isotopomer.df)
  Molecule_igraph@molecule_info$smiles <- unname(as.character(get_sdf_smiles(sdf)))
  return(Molecule_igraph)
}


get_Molecule_igraph_from_smiles <- function(smiles) {
  sdf <- get_smiles_sdf(smiles)
  if (length(sdf)==1) sdf <- sdf[[1]]
  get_Molecule_igraph_from_sdf(sdf)
}


#'  Add isotopomer
#' @describeIn Molecule_igraph Add isotopomer
#'
#' @param Molecule_igraph `Molecule_igraph`
#' @param isotopomer isotopomer name
#' @param iso_vec isotope of atom, such as c("C_1" = "[13]C")
#' @param abundance number
#'
#' @returns `Molecule_igraph`
#' @export
#'
Molecule_igraph_add_isotopomer <- function(
    Molecule_igraph , isotopomer = NULL,iso_vec = NULL,abundance=NA){

  if (identical(iso_vec,"all_C")) {
    iso_vec <- make_vector("[13]C",atom(glu.mi,"C"))
  }

  iso_label <- split(names(iso_vec),iso_vec)
  iso_label <- sapply(seq_along(iso_label),function(x){
    paste0("(",paste0(iso_label[[x]],collapse = ","),")",names(iso_label)[x])
  })%>%paste0(collapse = ";")
  ele_vec <-make_vector(vdata(Molecule_igraph)$element,
                        atom(Molecule_igraph))
  ele_vec[names(iso_vec)] <- iso_vec
  isotopologue <- sum(is.isotope(ele_vec),na.rm = T)
  isotopologue <- paste0("M",isotopologue)
  isotopomer.df <- Molecule_igraph@isotopomer
  if (is.null(isotopomer)){
    i <- 1
    while(paste0(isotopologue,"_",i) %in% isotopomer.df$isotopomer){
      i <- i+1
    }
    isotopomer <- paste0(isotopologue,"_",i)
  }
  to.add <- data.frame(isotopomer = isotopomer,
                       isotopologue = isotopologue,
                       label = iso_label ,
                       abundance = abundance
                       )
  rownames(to.add) <- to.add$isotopomer
  isotopomer.df[isotopomer,names(to.add)] <- to.add
  isotopomer.df[isotopomer,atom(Molecule_igraph)] <- ele_vec

  isotopomer.df -> Molecule_igraph@isotopomer
  return(Molecule_igraph)


}


get_Molecule_igraph_MS1 <- function(Molecule_igraph,polarity=1,adduct = NULL){

  if (is.null(adduct)  ) {
    adduct <- ifelse(polarity==0,"[M-H]-","[M+H]+")
  }


  ### ms spectra
  {

    mz.m0 <- MSCC::chemform_adduct(formula(Molecule_igraph),adduct = adduct)
    isotopomers.eles <- as.matrix(Molecule_igraph@isotopomer[,atom(Molecule_igraph)])
    ele.iso.diff <- make_vector(element_table$Mass_Dif,element_table$symbol)
    isotopomers.mz.diff <- ele.iso.diff[isotopomers.eles]
    dim(isotopomers.mz.diff) <- dim(isotopomers.eles)
    dimnames(isotopomers.mz.diff) <- dimnames(isotopomers.eles)
    mz.diff <- apply(isotopomers.mz.diff,1,sum)

    ms1.data <- Molecule_igraph@isotopomer%>%
      dplyr::mutate(mz = mz.m0+mz.diff)%>%
      dplyr::group_by(isotopologue)%>%
      dplyr::mutate(intensity = sum(abundance) )%>%
      dplyr::ungroup()%>%
      dplyr::select(any_of(c("isotopomer","isotopologue","label","abundance","mz","intensity")))

  nor.to <- ms1.data$intensity[which.min(ms1.data$mz)]
  ms1.data$intensity <-  ms1.data$intensity/nor.to*100
  }

  return(ms1.data)



}

get_Molecule_igraph_MS2 <- function(Molecule_igraph,cfmd){


  isotopomers <- Molecule_igraph@isotopomer
  FG.map <- cfmd@fragment_group_map

  #### prob
  {
    FG.map[FG.map<0.5] <- 0
    FG.map[FG.map>=0.5] <- 1
  }

  isotopomers.eles <- as.matrix(Molecule_igraph@isotopomer[,atom(Molecule_igraph)])
  ele.iso.diff <- make_vector(element_table$Mass_Dif,element_table$symbol)
  isotopomers.mz.diff <- ele.iso.diff[isotopomers.eles]
  dim(isotopomers.mz.diff) <- dim(isotopomers.eles)
  dimnames(isotopomers.mz.diff) <- dimnames(isotopomers.eles)

  isotopomers.mz.diff <- isotopomers.mz.diff[,colnames(FG.map)]
  fg.isotopomers.mz.diff <- isotopomers.mz.diff %*% t(FG.map)


  ms2.data <- cbind( cfmd@fragment_group,t(fg.isotopomers.mz.diff))%>%
    dplyr::filter(fragment_mz > 0)%>%
    tidyr::pivot_longer(rownames(isotopomers.eles),names_to = "isotopomer",values_to = "mzdiff")%>%
    dplyr::mutate( isotopomers[match(isotopomer,isotopomers$isotopomer),
                               c("isotopologue","label","abundance")],
                   mz = fragment_mz +mzdiff)%>%
    dplyr::group_by(mz)%>%
    dplyr::mutate(intensity = sum(abundance))



  return(ms2.data)


}


Molecule_igraph_get_C_order <- function(Molecule_igraph){

  dism <- igraph::distances(Molecule_igraph@igraph)
  ele <- element(Molecule_igraph)
  ele.dis.R <- apply(dism,1,function(x){
    sum(x[ele!="C"])
  })
  ele.dis.C <- apply(dism,1,function(x){
    sum(x[ele=="C"])
  })
  ele.dis.R[ele=="C"]
  ele.dis.C[ele=="C"]

}


vis_Molecule_igraph <- function(Molecule_igraph,show_id = F){

  Molecule_igraph.formated <- Molecule_igraph%>%
    Molecule_igraph_vis_format()%>%
    sdf_igraph_show_id(show_id)



  vis_igraph(Molecule_igraph.formated) %>%
    visNetwork::visPhysics(enabled = F) %>%
    visNetwork::visOptions(width = "100%", height = "100%")

}


Molecule_igraph_vis_format <- function(Molecule_igraph){

  vdata(Molecule_igraph) <-  vdata(Molecule_igraph)%>%
   # dplyr::group_by(atom)%>%
    dplyr::mutate(
      label = paste0(" ",atom," "),
      x =  ( x-mean(x))*100,
      y = ( y-mean(y)) *100,
      font.size = 30,
      borderWidth = 20,
      font.vadjust = 5,
      font.align = "center",
      color.border = "#AAAAAA",
      color.background = "#FFFFFF",
      borderWidth = 5,
      shape = "circle",
      physics = F)%>%
    dplyr::ungroup()


  edata(Molecule_igraph)  <- edata(Molecule_igraph)%>%
    dplyr::mutate(smooth = FALSE,
                  width = 10*bond_type)

  return(Molecule_igraph)

}



get_Molecule_atom_transfer_by_map <- function(mol.ig.from,
                                              mol.ig.to){



  sdf.parent <- mol.ig.from@sdf
  sdf.product <- mol.ig.to@sdf
  ig.parent <- mol.ig.from@igraph
  ig.product <- mol.ig.to@igraph

  ### mcs map
  {
    mcs <- fmcsR::fmcs(sdf.parent,
                       sdf.product,bu = 10)
    mcs.map <- get_mcs_atom_map(mcs)
    #mcs.map <- mcs.map.filter.duplicate(mcs.map,target_ele = get_ele_uniso(iso_ele))
  }


  ### complement mcs map
  {
    atom.map.matrix <- matrix(nrow = length(mcs.map),
                              ncol = length(atom(sdf.product)),
                              dimnames = list(seq_along(mcs.map),
                                              atom(sdf.product)))
    ring.diff <- length(rings(sdf.parent))- length(rings(sdf.product))
    bond.score <- rep(0,length(mcs.map))
    for (j in seq_along(mcs.map)) {
      this.map <- mcs.map[[j]]
      this.mapv <-this.map$mc1.atom
      names(this.mapv) <- this.map$mc2.atom
      this.mapv <- this.mapv[atom(sdf.product)]
      names(this.mapv) <- atom(sdf.product)


      ### ring re-assign
      {
        ring.solved <- F
        if (ring.diff&sum(is.na(this.mapv))) {
          ring.atom <- unname(unlist(rings(sdf.parent)))
          ring.atom.to.assign <- ring.atom[!ring.atom%in% this.mapv]
          ring.atom.to.assign <- unique(ring.atom.to.assign)
          adj <- sapply(ring.atom.to.assign,function(x){
            #x <- ring.atom.to.assign
            x.adj <- names(V(ig.parent))[distances(ig.parent,x)==1]
            x.adj <- x.adj[x.adj%in%this.mapv]
            y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
            y.candi <-apply(distances(ig.product,y.adj),1,function(z){
              zz <- names(z)[which(z==1)]
              zz[!zz%in% names(na.omit(this.mapv))&
                   str_extract(zz,"[:alpha:]*")==str_extract(x,"[:alpha:]*")]
            })%>%unlist()%>%unique()
            #message("ring candi: ",length(y.candi))
            unname(y.candi[1])
          })

          adj <- na.omit(unlist(adj))
          this.mapv[adj] <- names(adj)
          ring.solved <- ifelse(length(adj),T,F)
        }

      }

      ### ring nearest
      {
        if (ring.diff&sum(is.na(this.mapv))&ring.solved){
          #
          ring.nearest.to.assign <-apply(distances(ig.parent,ring.atom),
                                         1,function(z){
                                           zz <- names(z)[which(z==1)]
                                           zz[!zz%in% (na.omit(this.mapv))]
                                         })%>%unlist()%>%unique()

          adj <- sapply(ring.nearest.to.assign,function(x){
            #x <- ring.nearest.to.assign
            x.adj <- names(V(ig.parent))[distances(ig.parent,x)==1]
            x.adj <- x.adj[x.adj%in%this.mapv & x.adj%in%ring.atom]
            y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
            y.candi <-apply(distances(ig.product,y.adj),1,function(z){
              zz <- names(z)[which(z==1)]
              zz[!zz%in% names(na.omit(this.mapv))&
                   str_extract(zz,"[:alpha:]*")==str_extract(x,"[:alpha:]*")]
            })%>%unlist()%>%unique()
            #message("ring nearest candi: ",length(y.candi))
            unname(y.candi[1])
          })

          adj <- na.omit(unlist(adj))
          this.mapv[adj] <- names(adj)


        }

      }

      ### non-match nearest
      {

        if (sum(is.na(this.mapv))) {
          non.match.to.assign <- setdiff(atom(sdf.parent),this.mapv)
          adj <- sapply(non.match.to.assign,function(x){
            #x <- non.match.to.assign
            x.adj <- names(V(ig.parent))[distances(ig.parent,x)==1]
            x.adj <- x.adj[x.adj%in%this.mapv]
            y.adj <- names(this.mapv)[match(x.adj,this.mapv)]
            y.candi <-apply(distances(ig.product,y.adj),1,function(z){
              zz <- names(z)[which(z==1)]
              zz[!zz%in% names(na.omit(this.mapv))&
                   str_extract(zz,"[:alpha:]*")==str_extract(x,"[:alpha:]*")]
            })%>%unlist()%>%unique()
            #message("non match nearest candi: ",length(y.candi))
            unname(y.candi[1])
          })

          adj <- na.omit(unlist(adj))
          this.mapv[adj] <- names(adj)

        }

      }

      ### bond diff
      {
        temp.map <- na.omit(this.mapv)
        temp.map.t <- make_vector(names(temp.map),temp.map)
        ig.sub <- igraph_filter_distance(ig.parent,from = temp.map,dis = 1)

        m1 <- igraph::as_adjacency_matrix(ig.sub,
                                          attr = "bond_type",sparse = F)
        m1 <- m1+t(m1)
        m2 <- igraph::as_adjacency_matrix(ig.product,
                                          attr = "bond_type",sparse = F)[names(temp.map),names(temp.map)]
        m2 <- m2 + t(m2)
        m2 <- get_matrix_value_fill_with_NA(
          m2,temp.map.t[rownames(m1)],temp.map.t[colnames(m1)])
        m2[is.na(m2)] <- 0
        bond.score[j] <- 1- sum((m1-m2)!=0)/sum(m1!=0)


      }
      #p<-vis_sdf_igraph_compare(ig.parent,ig.product,temp.map,names(temp.map),show.label = T)

      atom.map.matrix[j,] <- this.mapv


    }
  }


  ### class
  {

    atom.map.sum <- data.frame(
      from.smiles = mol.ig.from@molecule_info$smiles,
      to.smiles = mol.ig.to@molecule_info$smiles,
      transfer = paste0("Atom_map",num2str(1:nrow(atom.map.matrix)))
    )

    rownames(atom.map.matrix) <- atom.map.sum$transfer


    new("Molecule_atom_transfer",
        transfer_def = atom.map.sum,
        transfer_matrix = atom.map.matrix)



  }

}




