get_MSIPFragmentMap_from_MSIPAtomMap <- function(msipAtomMap, iso_count_max = 0, target_ele){


  fg.map <- new("MSIPFragmentMap")

  cfmd.atom <- atom(get_Molecule_igraph_from_MSIPAtomMap(msipAtomMap),target_ele)

  fg.map@FG.atom.matrix <-  matrix(nrow = 0,ncol = length(cfmd.atom))
  colnames(fg.map@FG.atom.matrix) <- cfmd.atom

  fg.map@FG.ratio.matrix <-  matrix(nrow = 0,ncol = iso_count_max+1)
  colnames(fg.map@FG.ratio.matrix) <- format_isotopologue(0:iso_count_max,"M")


  return(fg.map)

}

get_MSIPFragmentMap_from_cfmd <- function(msipAtomMap, iso_count_max = 0, target_ele){
  .Deprecated("get_MSIPFragmentMap_from_MSIPAtomMap")
  get_MSIPFragmentMap_from_MSIPAtomMap(
    msipAtomMap = msipAtomMap,
    iso_count_max = iso_count_max,
    target_ele = target_ele
  )
}
