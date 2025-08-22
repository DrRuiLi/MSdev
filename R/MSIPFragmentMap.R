get_MSIPFragmentMap_from_cfmd <- function(cfmd,iso_count_max = 0,target_ele){


  fg.map <- new("MSIPFragmentMap")

  cfmd.atom <- atom(get_Molecule_igraph_from_cfmd(cfmd),target_ele)

  fg.map@FG.atom.matrix <-  matrix(nrow = 0,ncol = length(cfmd.atom))
  colnames(fg.map@FG.atom.matrix) <- cfmd.atom

  fg.map@FG.ratio.matrix <-  matrix(nrow = 0,ncol = iso_count_max+1)
  colnames(fg.map@FG.ratio.matrix) <- format_isotopologue(0:iso_count_max,"M")


  return(fg.map)

}
