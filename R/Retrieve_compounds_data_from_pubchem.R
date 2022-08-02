#' @title Retrieve compound data
#' @desription retrieve compound information by pubchem.cid in compound.record.file, this sheet should be the first sheet in xlsx file
#' @return compound.record
#' @Import  webchem
#' @import enviPat
#' @export
#' @examples
Retrieve_compounds_data_from_pubchem <-
  function(compound.record.file) {
    ### import
    {
      #compound.record.file <- "Standard.record.2021.12.18.STD_01.xlsx"
      compound.record.workbook <- openxlsx::loadWorkbook(compound.record.file)
      if (!"Compound_info" %in% names(compound.record.workbook)) {
        openxlsx::renameWorksheet(compound.record.workbook,1,"Compound_info")
      }

      compound.record <- openxlsx::read.xlsx(compound.record.workbook)

    }


    ### check key
    {
      if (is.null(compound.record$pubchem.cid)) {
        message("pubchem.cid not exist")
      }

    }

    ### get compound information from pubchem
    {
      compound.record$pubchem.cid <- gsub(pattern = "[^A-z0-9]",
                                          x = compound.record$pubchem.cid,
                                          replacement = "")
      pubchem.info <-
       pc_prop(
          compound.record$pubchem.cid,
          c(
            "Title",
            "InChi",
            "ExactMass",
            "MolecularFormula",
            "InChIKey",
            "IsotopeAtomCount",
            "Charge"
          )
        )




      compound.record$name <- pubchem.info$Title
      compound.record$inchi <- pubchem.info$InChI
      compound.record$inchikey <- pubchem.info$InChIKey
      compound.record$pubchem.cid <- pubchem.info$CID
      if (is.null(compound.record$formula )) {
        compound.record$formula <- pubchem.info$MolecularFormula
      }
      compound.record$IsotopeAtomCount <-
        pubchem.info$IsotopeAtomCount
      compound.record$Charge <- pubchem.info$Charge
      compound.record$exact.mass <-
        pubchem.info$ExactMass %>% as.numeric()

    }
    ### check
    {
      data("isotopes")
      formula.checked <-
        check_chemform(isotopes = isotopes, chemforms = compound.record$formula)
      is.contain.na <-
        check_ded(formulas = formula.checked$new_formula , deduct = "Na1") %>% as.logical()
      is.contain.cl <-
        check_ded(formulas = formula.checked$new_formula , deduct = "Cl1") %>% as.logical()
      compound.record$formula <- formula.checked$new_formula
      compound.record$is.salt <- (!is.contain.na)|(!is.contain.cl)
      compound.record$salt.to.remove <- NA
      compound.record$is.isotope <- !compound.record$exact.mass-formula.checked$monoisotopic_mass < 1e-6
      compound.record$isotope.to.replace <- NA
      if (sum(compound.record$is.salt)!= 0) {

        paste0(sum(compound.record$is.salt), " formula contain salt, plsease check ")%>%
          crayon::red()%>%
          message()
      }

      if (sum(compound.record$is.formula.not.equal.mass)!= 0) {
        paste0(sum(compound.record$is.formula.not.equal.mass),
               " formula and mass not match, plsease check ")%>%
          crayon::red()%>%
          message()
      }
    }

    ### sort and save xlsx
    {
      col.head <-
        c("name", "formula","is.salt","salt.to.remove",
          "is.formula.not.equal.mass" ,"exact.mass", "pubchem.cid", "inchikey")
      compound.record <- compound.record %>%
        select(col.head, everything())

      openxlsx::writeData(wb = compound.record.workbook,sheet = "Compound_info",x = compound.record)
      openxlsx::saveWorkbook(compound.record.workbook ,
                             file = compound.record.file,
                             overwrite = T)

    }
    return(compound.record)




  }
