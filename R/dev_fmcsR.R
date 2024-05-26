get_mcs_atom_map <- function(mcs){

  mcs.count <- length(mcs@mcs1[[2]])
  mcs1.atom <- rownames(atomblock(mcs@mcs1$query)[[1]])
  mcs2.atom <- rownames(atomblock(mcs@mcs2$target)[[1]])

  atom.map <- list()
  for (i in 1:mcs.count) {
    this.map <- data.frame(
      mc1.idx = mcs@mcs1$mcs1[[i]],
      mc2.idx = mcs@mcs2$mcs2[[i]]
    )
    this.map$mc1.atom <- mcs1.atom[this.map$mc1.idx]
    this.map$mc2.atom <- mcs2.atom[this.map$mc2.idx]
    atom.map[[i]] <- this.map
  }

  return(atom.map)

}
