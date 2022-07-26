#' @title Retrieve compound data
#' @desription retrieve compound information by pubchem.cid in compound.record.file
#' @return compound.record
#' @export
#' @import webchem
#' @examples
Retrieve_compounds_data_from_pubchem <-
  function(compound.record.file) {
    ### import
    {
      #compound.record.file <- "Standard.record.2021.12.18.STD_01.xlsx"
      compound.record <- readxl::read_excel(compound.record.file)

    }


    ### check key
    {
      if (is.null(compound.record$pubchem.cid)) {
        message("pubchem.cid not exist")
      }
      # if (is.null(compound.record$mzML.positive)) {
      #   message("Select positive mzML file")
      #   compound.record$mzML.positive <-
      #     choose.files(default = getwd(), caption = "Select positive")
      # }
      # if (is.null(compound.record$mzML.negative)) {
      #   message("Select negative mzML file")
      #   compound.record$mzML.negative <-
      #     choose.files(default = getwd(), caption = "Select negative")
      # }
    }

    ### get compound information from pubchem
    {
      pubchem.info <-
        webchem::pc_prop(
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
      compound.record$formula <- pubchem.info$MolecularFormula
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

      compound.record$formula <- formula.checked$new_formula
      compound.record$is.salt <- !is.contain.na


    }

    ### sort
    {
      col.head <-
        c("name", "formula", "exact.mass", "pubchem.cid", "inchikey")
      compound.record <- compound.record %>%
        select(col.head, everything())
      openxlsx::write.xlsx(compound.record , file = compound.record.file)

    }
    return(compound.record)




  }
