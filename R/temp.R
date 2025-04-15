

atomcount_temp <- function (x, addH = FALSE, ...)
{
  if (addH == TRUE) {
    return(table(c(gsub("_.*", "", rownames(x@atomblock)),
                   rep("H", bonds(x, type = "addNH")))))
  }
  else {
    return(table(gsub("_.*", "", rownames(x@atomblock))))
  }
}


bonds_temp <- function (x, type = "addNH")
{
  if (!any(c("SDF", "SDFset") %in% class(x)))
    stop("x needs to be of class SDF or SDFset")
  if (!any(c("bonds", "charge", "addNH") %in% type))
    stop("type can only be assigned: bonds, charge or addNH")
  .bonds <- function(x, type = type) {
    atomMA <- atomblock(x)
    atoms <- gsub("_.*", "", rownames(atomMA))
    bondMA <- bondblock(x)
    Nbonds1 <- cbind(atoms = c(bondMA[, 1], bondMA[, 2]),
                     bonds = c(bondMA[, 3], bondMA[, "C3"]))
    Nbonds1 <- tapply(Nbonds1[, "bonds"], Nbonds1[, "atoms"],
                      sum)
    Nbonds <- rep(0, length(atomMA[, 1]))
    names(Nbonds) <- seq(along = atomMA[, 1])
    Nbonds[names(Nbonds1)] <- Nbonds1
    val <- c(`1` = 1, `17` = 1, `2` = 2, `16` = 2, `13` = 3,
             `15` = 3, `14` = 4)
    group <- as.numeric(atomprop$Group)
    names(group) <- as.character(atomprop$Symbol)
    Nbondrule <- val[as.character(group[atoms])]
    Nbondrule[is.na(Nbondrule)] <- 0
    Nbondrule[Nbondrule < Nbonds] <- Nbonds[Nbondrule < Nbonds]
    charge <- c(`0` = 0, `1` = 3, `2` = 2, `3` = 1, `4` = 0,
                `5` = -1, `6` = -2, `7` = -3)
    charge <- charge[as.character(atomMA[, 5])]
    Nbonds <- data.frame(atom = atoms, Nbondcount = Nbonds,
                         Nbondrule = Nbondrule, charge = charge)
    if (type == "bonds") {
      return(Nbonds)
    }
    if (type == "charge") {
      chargeindex <- Nbonds[, "charge"] != 0
      if (sum(chargeindex) == 0) {
        return(NULL)
      }
      else {
        chargeDF <- Nbonds[chargeindex, ]
        charge <- chargeDF[, "charge"]
        names(charge) <- chargeDF[, "atom"]
        return(charge)
      }
    }
    if (type == "addNH") {
      Nbonds[Nbonds[, "Nbondcount"] >= Nbonds[, "Nbondrule"],
             "charge"] <- 0
      Nbonds[Nbonds[, "Nbondcount"] == 0, c("Nbondrule",
                                            "charge")] <- 0
      NH <- sum((Nbonds[, "Nbondrule"] + Nbonds[, "charge"]) -
                  Nbonds[, "Nbondcount"])
      if (NH < 0)
        NH <- 0
      return(NH)
    }
  }
  if (class(x) == "SDF") {
    bonds <- .bonds(x, type)
    return(bonds)
  }
  if (class(x) == "SDFset") {
    bonds_set <- lapply(seq(along = x), function(y) .bonds(x[[y]],
                                                           type))
    names(bonds_set) <- cid(x)
    if (type == "bonds") {
      return(bonds_set)
    }
    if (type == "charge") {
      return(bonds_set)
    }
    if (type == "addNH") {
      return(unlist(bonds_set))
    }
  }
}
