#' MSIPAtomMap Class
#'
#' @title MSIPAtomMap Class for Atom Mapping in Mass Spectrometry
#' @description An S4 class that stores fragment structure information including
#' igraph representations, SDF objects, atom-to-atom mappings, and fragment groupings.
#' This class contains all the data needed for atom tracing analysis.
#'
#' @slot peak_assignment Data frame containing peak assignment data with fragment mappings
#' @slot fragment_define Data frame containing fragment definitions with SMILES and formulas
#' @slot fragment_transition Data frame containing fragment transition relationships
#' @slot fragment_igraph List of Molecule_igraph objects representing molecular graphs for each fragment
#' @slot fragment_sdf SDFset containing SDF representations of fragments
#' @slot fragment_atom_map List of matrices containing atom-to-atom mapping probabilities
#' @slot fragment_group Data frame containing fragment group definitions and statistics
#' @slot fragment_group_map Matrix mapping fragment groups to atoms
#'
#' @export
setClass("MSIPAtomMap",
         slots = list(
           peak_assignment = "data.frame",
           fragment_define = "data.frame",
           fragment_transition = "data.frame",
           fragment_igraph = "list",
           fragment_sdf = "SDFset",
           fragment_atom_map = "list",
           fragment_group = "data.frame",
           fragment_group_map = "matrix"
         ),
         prototype = list(
           peak_assignment = data.frame(),
           fragment_define = data.frame(),
           fragment_transition = data.frame(),
           fragment_igraph = list(),
           fragment_sdf = new("SDFset"),
           fragment_atom_map = list(),
           fragment_group = data.frame(),
           fragment_group_map = matrix()
         ))


#' Show method for MSIPAtomMap
#' @param object An MSIPAtomMap object
#' @export
setMethod("show", signature = "MSIPAtomMap", definition =
            function(object){
              n_frags <- length(object@fragment_igraph)
              n_groups <- if(nrow(object@fragment_group) > 0) nrow(object@fragment_group) else 0
              polarity_str <- if(nrow(object@fragment_define) > 0 && "polarity" %in% colnames(object@fragment_define)) {
                pol <- object@fragment_define$polarity[1]
                if(!is.na(pol)) {
                  if(pol == 0) "negative" else "positive"
                } else {
                  "unknown"
                }
              } else {
                "unknown"
              }
              msg <- paste0("MSIPAtomMap (", polarity_str, ") with ", n_frags, " fragments")
              if(n_groups > 0) msg <- paste0(msg, " in ", n_groups, " groups")
              message(msg)
            })


#' Create MSIPAtomMap from CFM_data
#' @title Create MSIPAtomMap from CFM_data
#' @description Initializes an MSIPAtomMap object from a CFM_data object.
#'
#' @param cfm_data A CFM_data object
#'
#' @return A new MSIPAtomMap object with data from CFM_data
#' @export
MSIPAtomMap_from_CFM_data <- function(cfm_data){
  object <- new("MSIPAtomMap")
  object@peak_assignment <- cfm_data@peak_assignment
  object@fragment_define <- cfm_data@fragment_define
  object@fragment_transition <- cfm_data@fragment_transition
  return(object)
}


#' Add Seed Fragment to MSIPAtomMap
#' @title Add Seed Fragment to MSIPAtomMap
#' @description Adds a seed fragment (Fragment00) with the original molecule SMILES
#' to the MSIPAtomMap object. This fragment represents the intact parent molecule
#' before any fragmentation occurs.
#'
#' @param object An MSIPAtomMap object
#' @param smiles SMILES string of the original molecule
#'
#' @return The updated MSIPAtomMap object with seed fragment added
#' @export
MSIPAtomMap_add_seed <- function(object, smiles){

  # Get the polarity suffix
  polarity_suffix <- ""
  if (nrow(object@fragment_define) > 0 && "polarity" %in% colnames(object@fragment_define)) {
    polarity <- object@fragment_define$polarity[1]
    if (!is.na(polarity)) {
      polarity_suffix <- ifelse(polarity == 0, "_0", "_1")
    }
  }

  # Create seed fragment row (Fragment00 with polarity suffix)
  seed_fragment <- data.frame(
    fragment_id = paste0("Fragment00", polarity_suffix),
    fragment_mz = 0,
    smiles = smiles,
    fragment_group = NA,
    formula = get_smile_formula(smiles)
  )

  # Add polarity column if it exists in fragment_define
  if ("polarity" %in% colnames(object@fragment_define)) {
    seed_fragment$polarity <- object@fragment_define$polarity[1]
  }

  # Set row names
  rownames(seed_fragment) <- seed_fragment$fragment_id

  # Add seed fragment to fragment_define (prepend)
  object@fragment_define <- dplyr::bind_rows(seed_fragment, object@fragment_define)

  # Get first fragment ID for transition
  first_frag_id <- object@fragment_define$fragment_id[2]  # Second row is now the original first fragment

  # Create transition from seed to first fragment
  seed_transition <- data.frame(
    from = paste0("Fragment00", polarity_suffix),
    to = first_frag_id
  )

  # Add polarity column if it exists in fragment_transition
  if ("polarity" %in% colnames(object@fragment_transition)) {
    seed_transition$polarity <- object@fragment_transition$polarity[1]
  }

  # Add seed transition to fragment_transition (prepend)
  object@fragment_transition <- dplyr::bind_rows(seed_transition, object@fragment_transition)

  return(object)
}


#' Get Igraph Objects for MSIPAtomMap Fragments
#' @title Get Igraph Objects for MSIPAtomMap Fragments
#' @description Converts fragment SMILES structures to igraph objects representing
#' molecular graphs. This function processes fragment SMILES strings to SDF format and then
#' creates igraph objects for structural analysis.
#' @describeIn MSIPAtomMap get fragment igraph objects
#'
#' @param object An MSIPAtomMap object
#'
#' @return The input MSIPAtomMap object with updated fragment_igraph, fragment_sdf, and fragment_define slots
#' @export
MSIPAtomMap_get_igraph <- function(object){

  ### fragment def
  {
    fragment.data <- object@fragment_define
    fragment.sdf <- suppressWarnings(get_smiles_sdf(fragment.data$smiles))
    fragment.data$formula <- get_sdf_formula(fragment.sdf)
  }

  ### filter valid smiles and sdf
  {
    fragment.data <- fragment.data %>%
      dplyr::filter(!is.na(formula))
    object@fragment_define <- fragment.data
    object@fragment_transition <-
      object@fragment_transition %>%
      dplyr::filter(from %in% object@fragment_define$fragment_id,
                    to %in% object@fragment_define$fragment_id)
  }

  ### igraph
  {
    fragment.sdf <- suppressWarnings(get_smiles_sdf(fragment.data$smiles))
    cid(fragment.sdf) <- fragment.data$fragment_id
    fragment.igraph <- get_sdf_igraph(fragment.sdf)
    names(fragment.igraph) <- fragment.data$fragment_id
  }

  ### save to object
  {
    object@fragment_igraph <- fragment.igraph
    object@fragment_sdf <- fragment.sdf
  }

  return(object)
}


#' Get Atom Tracing Maps for MSIPAtomMap
#' @title Get Atom Tracing Maps for MSIPAtomMap
#' @description Computes atom-to-atom mapping between parent and product fragments.
#' This function traces how atoms in the original molecule are distributed across fragments
#' through fragmentation transitions, which is essential for isotope tracing analysis.
#' @describeIn MSIPAtomMap compute atom maps
#'
#' @note MSIPAtomMap_get_atom_map can be computationally intensive, particularly in the
#' trans map step which uses the fmcs (Maximum Common Substructure) algorithm.
#'
#' @param object An MSIPAtomMap object containing fragment_igraph and fragment_sdf
#' @param iso_ele Isotope element specification for tracing, e.g., "\[13\]C" (default: "\[13\]C")
#' @param BPPARAM A BiocParallel backend for parallel processing (default: BiocParallel::SerialParam())
#'
#' @return The input MSIPAtomMap object with updated fragment_atom_map slot
#' @export
MSIPAtomMap_get_atom_map <- function(object,
                                     iso_ele = "[13]C",
                                     BPPARAM = SerialParam()){

  ### trans.map
  {
    message_with_time("trans map")
    fragment.trans <- object@fragment_transition
    trans.maps <- bplapply(1:nrow(fragment.trans),
                           function(trans_id, object, fragment.trans, iso_ele){
                             fragment.igraph <- object@fragment_igraph
                             fragment.sdf <- object@fragment_sdf
                             ig.parent <- fragment.igraph[[fragment.trans$from[trans_id]]]
                             ig.product <- fragment.igraph[[fragment.trans$to[trans_id]]]
                             sdf.parent <- fragment.sdf[[fragment.trans$from[trans_id]]]
                             sdf.product <- fragment.sdf[[fragment.trans$to[trans_id]]]
                             maps <- get_atom_map(sdf.parent, sdf.product,
                                                  ig.parent, ig.product, iso_ele = iso_ele)
                             return(maps)
                           },
                           object = object,
                           fragment.trans = fragment.trans,
                           iso_ele = iso_ele,
                           BPPARAM = BPPARAM)
    names(trans.maps) <- paste0(fragment.trans$from, fragment.trans$to)
    trans.maps.stat <- MSIPAtomMap_check_trans_map(object, trans.maps = trans.maps, iso_ele = iso_ele)
    fragment.trans[, colnames(trans.maps.stat)] <- trans.maps.stat
    object@fragment_transition <- fragment.trans
  }

  ### filter invalid trans
  {
    object_bak <- object
    object <- MSIPAtomMap_remove_trans(object)
  }

  ### fragment map
  {
    message_with_time("fragment map")
    fragment.atom.map <- list()
    ### weight for path selection
    fragment.trans <- object@fragment_transition %>%
      dplyr::mutate(id = paste0(from, to),
                    weight = loss.distance)
    fragment.data <- object@fragment_define
    fragment.data$ratio <- NA
    fragment.data$bond.score <- NA
    fragment.data$cumsum.loss.distance <- NA
    fragment.igraph <- object@fragment_igraph
    if (nrow(fragment.trans)) {
      ig.trans <- igraph::graph_from_data_frame(fragment.trans)
      for (i in 1:nrow(fragment.data)) {

        this.frag <- fragment.data$fragment_id[i]
        fragment.atom.map[[i]] <- NA
        if(i == 1){
          ele <- get_sdf_igraph_atom(fragment.igraph[[1]])
          maps <- diag(nrow = length(ele))
          rownames(maps) <- colnames(maps) <- ele
          fragment.atom.map[[i]] <- maps
          fragment.data$ratio[i] <- 1
          next
        }

        ### find path
        {
          this.paths <- igraph::all_shortest_paths(ig.trans,
                                                   from = 1,
                                                   to = this.frag,
                                                   mode = "out")$vpath
          this.epaths <- lapply(this.paths, function(path){
            epath <- paste0(names(path), names(path)[-1])
            epath[1:length(epath) - 1]
          })
          this.path.ratio <- sapply(this.epaths, function(this.epath){
            sapply(this.epath, function(epath){
              idx.path <- match(epath, fragment.trans$id)
              r <- fragment.trans$ratio[idx.path]
              return(r)
            }) %>% prod()
          })
          this.path.bond.score <- sapply(this.epaths, function(this.epath){
            sapply(this.epath, function(epath){
              idx.path <- match(epath, fragment.trans$id)
              r <- fragment.trans$bond.score[idx.path]
              return(r)
            }) %>% prod()
          })
          this.path.cumsum.loss.distance <- sapply(this.epaths, function(this.epath){
            sapply(this.epath, function(epath){
              idx.path <- match(epath, fragment.trans$id)
              r <- fragment.trans$loss.distance[idx.path]
              return(r)
            }) %>% sum()
          })
        }

        ### prod maps
        {
          maps.list <- lapply(this.epaths, function(this.epath){
            if (!length(this.epath) == 0) {
              maps <- trans.maps[this.epath]
              while(length(maps) > 1){
                maps[[2]] <- maps[[1]] %*% maps[[2]]
                maps[[1]] <- NULL
              }
              return(maps[[1]])
            }
          })
          map <- do.call(sum_matrix, maps.list) / length(maps.list)
        }

        ### return
        {
          fragment.atom.map[[i]] <- map
          fragment.data$ratio[i] <- mean(this.path.ratio)
          fragment.data$bond.score[i] <- mean(this.path.bond.score)
          fragment.data$cumsum.loss.distance[i] <- mean(this.path.cumsum.loss.distance)
        }
      }
      names(fragment.atom.map) <- fragment.data$fragment_id
    }
  }

  object@fragment_atom_map <- fragment.atom.map
  object@fragment_define <- fragment.data

  return(object)
}


#' Remove Invalid Transitions from MSIPAtomMap
#' @title Remove Invalid Transitions from MSIPAtomMap
#' @description Removes transitions marked as invalid and filters out unreachable fragments.
#'
#' @param object An MSIPAtomMap object
#'
#' @return The updated MSIPAtomMap object
#' @export
MSIPAtomMap_remove_trans <- function(object){

  .f <- function(object){
    fragment.trans <- object@fragment_transition
    fragment.trans <- fragment.trans[fragment.trans$volid, ]

    trans.ig <- igraph::graph_from_data_frame(fragment.trans)
    dis.to.fragment1 <- distances(trans.ig, mode = "out",
                                  v = object@fragment_define$fragment_id[1])
    reachable <- colnames(dis.to.fragment1)[!is.infinite(dis.to.fragment1)]
    to.remove <- setdiff(object@fragment_define$fragment_id, reachable)
    ### remove
    {
      object@peak_assignment <- object@peak_assignment %>%
        dplyr::filter(!fragment_id %in% to.remove)

      object@fragment_define <- object@fragment_define %>%
        dplyr::filter(!fragment_id %in% to.remove)

      object@fragment_transition <- fragment.trans %>%
        dplyr::filter(!(from %in% to.remove | to %in% to.remove))

      object@fragment_igraph <- object@fragment_igraph[
        !names(object@fragment_igraph) %in% to.remove]

      object@fragment_sdf <- object@fragment_sdf[
        !cid(object@fragment_sdf) %in% to.remove]

      object@fragment_atom_map <- object@fragment_atom_map[
        !names(object@fragment_atom_map) %in% to.remove]
    }

    return(object)
  }

  ### iteration
  i <- 1
  object <- .f(object)
  while(any(is.infinite(distances(get_MSIPAtomMap_trans_igraph(object), 1, mode = "out"))) & i <= 5){
    object <- .f(object)
  }
  if (i == 5)
    warning("MSIPAtomMap_remove_trans abnormal")

  return(object)
}


#' Check Transition Atom Maps
#' @title Check Transition Atom Maps
#' @description Validates and computes statistics for atom mapping transitions.
#'
#' @param object An MSIPAtomMap object
#' @param iso_ele Isotope element specification (default: "\[13\]C")
#' @param trans.maps List of transition mapping matrices (optional)
#'
#' @return Data frame with transition mapping statistics
#' @export
MSIPAtomMap_check_trans_map <- function(object,
                                        iso_ele = "[13]C",
                                        trans.maps = NULL){

  if (is.null(trans.maps)) {
    trans.maps <- object@fragment_atom_map
  }
  bond.score <- sapply(trans.maps, function(x) attributes(x)$bond.score)
  atom.ele <- get_ele_uniso(iso_ele)
  trans.maps <- lapply(trans.maps, function(x){
    x[grepl(atom.ele, rownames(x)), grepl(atom.ele, colnames(x)), drop = F]
  })
  n_parent <- sapply(trans.maps, function(x){
    nrow(x)
  })
  n_atoms <- sapply(trans.maps, function(x){
    ncol(x)
  })
  n_atoms_compose_map <- sapply(trans.maps, function(x){
    sum(rowSums(x) == 1)
  })
  n_atoms_certain_map <- sapply(trans.maps, function(x){
    max.prob <- apply(x, 2, max)
    sum(max.prob == 1)
  })
  n_atoms_noncertain_map <- sapply(trans.maps, function(x){
    max.prob <- apply(x, 2, max)
    sum(max.prob < 1 & max.prob > 0)
  })
  n_atoms_non_map <- sapply(trans.maps, function(x){
    max.prob <- apply(x, 2, max)
    sum(max.prob == 0)
  })

  map.stat <- data.frame(n_atoms,
                         n_parent,
                         n_atoms_compose_map,
                         n_atoms_certain_map,
                         n_atoms_noncertain_map,
                         n_atoms_non_map,
                         bond.score) %>%
    dplyr::mutate(volid = n_atoms_non_map == 0,
                  ratio = n_atoms_compose_map / n_atoms,
                  atom.loss = n_atoms - n_atoms_compose_map,
                  bond.loss = 1 - bond.score,
                  loss.distance = atom.loss + bond.loss)

  return(map.stat)
}


#' Get Fragment Group Map for MSIPAtomMap
#' @title Get Fragment Group Map for MSIPAtomMap
#' @description Groups fragments by m/z and computes atom mapping for each group.
#'
#' @param object An MSIPAtomMap object
#' @param iso_ele Isotope element specification (default: "\[13\]C")
#' @param ppm Mass tolerance in ppm for grouping (default: 5)
#'
#' @return The updated MSIPAtomMap object
#' @export
MSIPAtomMap_get_FG_map <- function(object, iso_ele = "[13]C", ppm = 5){

  ### Get polarity suffix
  {
    polarity_suffix <- ""
    if (nrow(object@fragment_define) > 0 && "polarity" %in% colnames(object@fragment_define)) {
      polarity <- object@fragment_define$polarity[1]
      if (!is.na(polarity)) {
        polarity_suffix <- ifelse(polarity == 0, "_0", "_1")
      }
    }
  }

  ### Fragment group
  {
    fg <- groupMz(object@fragment_define$fragment_mz, ppm)
    object@fragment_define$fragment_group <- paste0("FG", num2str(fg), polarity_suffix)
    object@peak_assignment$fragment_group <-
      object@fragment_define$fragment_group[match(
        object@peak_assignment$fragment_id,
        object@fragment_define$fragment_id)]
    fg.count <- table(object@fragment_define$fragment_group)
    fg.df <- data.frame(
      fragment_group = paste0("FG", num2str(sort(unique(fg))), polarity_suffix)
    ) %>%
      dplyr::mutate(
        fragment_count = as.numeric(fg.count[fragment_group]),
        fragment_mz = object@fragment_define$fragment_mz[match(fragment_group,
                                                                object@fragment_define$fragment_group)]
      )
  }

  ### FG map
  {
    target_atoms <- get_sdf_igraph_atom(get_MSIPAtomMap_sdf_igraph(object), get_ele_uniso(iso_ele))
    frag.atom.matrix <- matrix(ncol = length(target_atoms),
                               nrow = nrow(fg.df),
                               dimnames = list(fg.df$fragment_group,
                                               target_atoms))
    for (i.fg in seq_len(nrow(fg.df))) {

      this.frag.group <- fg.df$fragment_group[i.fg]
      this.frags <- object@fragment_define[object@fragment_define$fragment_group == this.frag.group, ]
      this.frag.atom <- get_MSIPAtomMap_fragment_group_atom_map(object, this.frag.group)
      this.frag.c <- this.frag.atom[target_atoms]
      frag.atom.matrix[this.frag.group, names(this.frag.c)] <- this.frag.c
    }
  }

  ### stat
  {
    frag.certainty <- apply(frag.atom.matrix, 1, function(x){
      sum(x == 1) / sum(x)
    })
    frag.certainty[rowSums(frag.atom.matrix) == 0] <- 0
    fg.df$certainty <- frag.certainty
  }

  object@fragment_group <- fg.df
  object@fragment_group_map <- frag.atom.matrix

  return(object)
}


#' Get Fragment Group Atom Map
#' @title Get Fragment Group Atom Map
#' @description Computes atom probability mapping for a specific fragment group.
#'
#' @param object An MSIPAtomMap object
#' @param frag.group Fragment group identifier
#'
#' @return Named numeric vector of atom probabilities
#' @export
get_MSIPAtomMap_fragment_group_atom_map <- function(object, frag.group){

  frag.idx <- which(object@fragment_define$fragment_group == frag.group)

  if (1 %in% frag.idx) {
    ele <- get_sdf_igraph_atom(get_MSIPAtomMap_sdf_igraph(object))
    frag.atoms.prob <- rep(1, length(ele))
    names(frag.atoms.prob) <- ele
  } else {
    frag.def <- object@fragment_define[frag.idx, ]
    frag.maps <- object@fragment_atom_map[frag.idx]
    frag.maps <- frag.maps[!sapply(frag.maps, is.null)]
    if (!length(frag.maps)) {
      ele <- get_sdf_igraph_atom(get_MSIPAtomMap_sdf_igraph(object))
      frag.atoms.prob <- rep(1, length(ele))
      names(frag.atoms.prob) <- ele
    } else {
      frag.atoms.prob <- sapply(frag.maps, rowSums) %>%
        rowMeans()
    }
  }

  return(frag.atoms.prob)
}


#' Get SDF Igraph from MSIPAtomMap
#' @title Get SDF Igraph from MSIPAtomMap
#' @description Retrieves the igraph object for a specific fragment.
#'
#' @param object An MSIPAtomMap object
#' @param fragment_id Fragment index or name (default: 1)
#'
#' @return An igraph object for the specified fragment
#' @export
get_MSIPAtomMap_sdf_igraph <- function(object, fragment_id = 1){
  object@fragment_igraph[[fragment_id]]
}


#' Get Transition Igraph for MSIPAtomMap Visualization
#' @title Get Transition Igraph for MSIPAtomMap
#' @description Creates an igraph object representing fragment transitions with
#' visualization attributes for use with visNetwork.
#'
#' @param object An MSIPAtomMap object
#'
#' @return An igraph object with node and edge attributes for visualization
#' @export
get_MSIPAtomMap_trans_igraph <- function(object){

  node.df <- object@fragment_define %>%
    dplyr::mutate(id = fragment_id,
                  label = id,
                  no = 1:n(),
                  color.border = case_when(
                    no == 1 ~ "rgba(100, 100, 100, 0.8)",
                    T ~ "rgba(100, 100, 100, 0.8)"
                  ),
                  color.background = case_when(
                    T ~ "rgba(255, 255, 255, 0.8)"
                  ),
                  font.size = 30,
                  borderWidth = 5,
                  size = case_when(no == 1 ~ 100,
                                   T ~ 50))
  edge.df <- object@fragment_transition %>%
    dplyr::mutate(match.str = paste0(from, to),
                  arrows.to.scaleFactor = 2,
                  atom.loss = n_atoms - n_atoms_compose_map,
                  bond.loss = 1 - bond.score,
                  color = "rgba(100, 100, 100, 0.2)",
                  width = 10,
                  smooth = F,
                  loss.distance = atom.loss + bond.loss,
                  length = normalize_max_min(loss.distance) * 1500,
                  length = 700,
                  no = 1:n())

  frag.trans.graph <- igraph::graph_from_data_frame(edge.df,
                                                     vertices = node.df)
  return(frag.trans.graph)
}


#' Get MSIPAtomMap from CFM_data
#' @title Get MSIPAtomMap from CFM_data
#' @description Creates MSIPAtomMap object from CFM_data by adding seed fragment,
#' igraph, atom maps, and fragment groups.
#'
#' @param cfm_data A CFM_data object
#' @param smiles SMILES string of the original molecule (required for adding seed fragment)
#' @param iso_ele Isotope element specification (default: "\[13\]C")
#'
#' @return An MSIPAtomMap object containing fragment data and atom mappings
#' @export
get_MSIPAtomMap_from_cfmd <- function(cfm_data,
                                      smiles = NULL,
                                      iso_ele = "[13]C"){

  message_with_time("Creating MSIPAtomMap from CFM_data")

  # Create MSIPAtomMap from CFM_data
  msipAtomMap <- MSIPAtomMap_from_CFM_data(cfm_data)

  # Add seed fragment if smiles is provided
  if (!is.null(smiles)) {
    msipAtomMap <- MSIPAtomMap_add_seed(msipAtomMap, smiles)
  }

  msipAtomMap <- MSIPAtomMap_get_igraph(msipAtomMap)

  message_with_time("MSIPAtomMap_get_atom_map")
  msipAtomMap <- MSIPAtomMap_get_atom_map(msipAtomMap, iso_ele = iso_ele)
  msipAtomMap <- MSIPAtomMap_get_FG_map(msipAtomMap, iso_ele = iso_ele)
  message_with_time("Done")

  return(msipAtomMap)
}


#' Get MSIPAtomMap from SMILES
#' @title Get MSIPAtomMap from SMILES
#' @description Wrapper function that creates MSIPAtomMap object from SMILES by running
#' CFM prediction/annotation and then constructing the atom map. Includes caching.
#'
#' @param smiles SMILES string of the molecule
#' @param compound_id Identifier for the compound (default: "temp_id")
#' @param ppm Mass tolerance in ppm (default: 5)
#' @param adduct Adduct type (default: "\[M+H\]+")
#' @param iso_ele Isotope element specification (default: "\[13\]C")
#' @param check_temp Whether to check/use temporary files (default: TRUE)
#' @param temp_dir Directory for temporary files
#' @param ... Additional arguments
#'
#' @return An MSIPAtomMap object containing fragment data and atom mappings
#' @export
get_MSIPAtomMap_from_smiles <- function(smiles = "NCC(O)=O",
                                        compound_id = "temp_id",
                                        ppm = 5,
                                        adduct = "[M+H]+",
                                        iso_ele = "[13]C",
                                        check_temp = T,
                                        temp_dir = get_dir_expand_from_onedrive("/Code/R/data/MSDB/CompoundDB/CFM_predicted_kegg.compdb_cfmd"),
                                        ...){
  # Check cache
  if(check_temp){
    if (!dir.exists(temp_dir)) dir.create(temp_dir,recursive = T,showWarnings = F)
    temp_file <- paste0(temp_dir,"/",compound_id,"_",adduct,"_msipAtomMap.rds")
    if(file.exists(temp_file)){
      message_with_time("loading MSIPAtomMap from temp:",temp_file)
      msipAtomMap <- readRDS(temp_file)
      return(msipAtomMap)
    }
  }

  start.time <- Sys.time()

  # Step 1: Get CFM_data
  cfm_data <- get_CFM_data_from_smiles(smiles = smiles,
                                        compound_id = compound_id,
                                        ppm = ppm,
                                        adduct = adduct)

  # Step 2: Create MSIPAtomMap from CFM_data
  msipAtomMap <- get_MSIPAtomMap_from_cfmd(cfm_data = cfm_data,
                                            smiles = smiles,
                                            iso_ele = iso_ele)

  map.time <- (Sys.time()-start.time)%>%as.numeric(units = "mins")

  # Save to cache
  if(check_temp){
    dir.create(temp_dir,showWarnings = F,recursive = T)
    temp_file <- paste0(temp_dir,"/",compound_id,"_",adduct,"_msipAtomMap.rds")
    saveRDS(msipAtomMap,file = temp_file)

    ### log
    log.info <- c()
    log.info["smiles"] <- smiles
    log.info["map.time"] <- map.time
    log.info["atom.count"] <- length(get_sdf_igraph_atom(get_MSIPAtomMap_sdf_igraph(msipAtomMap)))
    cat(paste0(paste0(log.info,collapse = ","),"\n"),
        file = paste0(temp_dir,"/atm_msip.log"),append = T)
  }

  return(msipAtomMap)
}
