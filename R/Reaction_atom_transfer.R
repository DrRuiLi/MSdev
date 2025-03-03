Metabolic_flux_atom_transfer_by_reaction <- function(rat,
                                                     mol.igs,
                                                     target_ele = "C"
){

  ### map
  {
    atom_transfer.data <- rat@atom_transfer%>%
      dplyr::mutate(from.id = paste0(from.molecule.id,"_",from.atom.id),
                    to.id = paste0(to.molecule.id,"_",to.atom.id)
                    )


  }

  ### from combination
  {

    cp.mol <- atom_transfer.data%>%
      dplyr::distinct(from.compound.id,from.molecule.id)%>%
      dplyr::pull(from.compound.id,from.molecule.id)
    from.id <- names(rat@reaction_info$from.smiles)
    from.id <- from.id[from.id%in% names(cp.mol)]
    available.isotopomers <- lapply(from.id,
           function(m.id){
             cp.id <- cp.mol[m.id]
             cp.isotopomer <- mol.igs[[cp.id]]@isotopomer$isotopomer
             cp.isotopomer
           }
           )
    names(available.isotopomers) <- from.id

    isotopomer.combn <- do.call(expand.grid,available.isotopomers)

  }

  ### atom transfer matrix
  {
    to.atom.matrix <- matrix(nrow = nrow(isotopomer.combn),
                             ncol = nrow(atom_transfer.data),
                             dimnames = list( paste0("combn",1:nrow(isotopomer.combn)),
                                             atom_transfer.data$to.id))
    path.record <- c()
    for (i in 1:nrow(isotopomer.combn)) {

      all.atoms <- c()
      paths <- list()
      for (j in colnames(isotopomer.combn)) {
        cp.id <- cp.mol[j]
        mol.ig <- mol.igs[[cp.id]]

        mol.atoms <- mol.ig@isotopomer[isotopomer.combn[i,j],atom(mol.ig,target_ele),drop = F]%>%
          unlist()
        names(mol.atoms) <- paste0(j,"_",names(mol.atoms))
        all.atoms <- c(all.atoms,mol.atoms)

        path <-  mol.ig@isotopomer[isotopomer.combn[i,j],"path"]
        path <- ifelse(is.null(path),"",path)
        paths[j] <- path

      }
      path.record[i] <- do.call(path_merge,paths)
      atom_transfer.data.this.combn <- atom_transfer.data%>%
        dplyr::mutate(atom.type = all.atoms[ from.id])
      to.atom.matrix[i, ] <- atom_transfer.data.this.combn$atom.type
    }



  }

  ### add labeled to.cp to mol.igs
  {
    #to.atom.matrix <- to.atom.matrix[!duplicated(to.atom.matrix),,drop = F]
    combn.labeled <- apply(to.atom.matrix,1,function(x){any(is.isotope(x))})

    to.cp.mol <- atom_transfer.data%>%
      dplyr::distinct(to.compound.id,to.molecule.id)%>%
      dplyr::pull(to.compound.id,to.molecule.id)
    to.mol.igs <- mol.igs[unique(to.cp.mol)]
    for (i in which(combn.labeled)) {

      atom_transfer.data.this.combn <- atom_transfer.data%>%
        dplyr::mutate(atom.type = to.atom.matrix[ i,])
      this.path <- path_merge(path.record[i],rat@reaction_info$reaction_id)

      for (i.m in unique(atom_transfer.data.this.combn$to.molecule.id)) {
        cp.id <- to.cp.mol[i.m]
        this.mig <- to.mol.igs[[cp.id]]

        iso_vec <- atom_transfer.data.this.combn %>%
          dplyr::filter(to.molecule.id== i.m)%>%
          dplyr::pull(atom.type,to.atom.id)
        if (any(is.isotope(iso_vec))) {
          iso.name <- get_isotopomer_name(iso_vec)
          if (any(this.mig@isotopomer$label==iso.name #&this.mig@isotopomer$path ==this.path
                  ,na.rm = T)) {

          }else{
            this.mig <- Molecule_igraph_add_isotopomer(
              this.mig,iso_vec  = iso_vec,path = this.path
            )
          }

        }

        this.mig -> to.mol.igs[[cp.id]]

      }

    }

  }

  ### merge molecule to compound
  {
    return(to.mol.igs)
  }

}


Reaction_atom_transfer_reverse <- function(rat){

  rat.r <- rat


  ### reaction_info
  {
    rat.r@reaction_info$from.smiles <- rat@reaction_info$to.smiles
    rat.r@reaction_info$to.smiles <- rat@reaction_info$from.smiles
    rat.r@reaction_info$reversed <- T
  }


  ### atom_transfer
  {
    x <- colnames(rat.r@atom_transfer)
    x<- sub(pattern = "^from",replacement = "tempstr",x = x)
    x<- sub(pattern = "^to",replacement = "from",x = x)
    x<- sub(pattern = "^tempstr",replacement = "to",x = x)

    x -> colnames(rat.r@atom_transfer)
  }

  return(rat.r)

}



path_split <- function(p){

  p <- ifelse(is.null(p),"",p)
  str_split(p,";")[[1]]

}

path_merge <- function(...){

  p <- str_split(c(...),";")%>%
    unlist()%>%
    setdiff("")
  p <- sort(p)
  paste0(p,collapse = ";")

}
