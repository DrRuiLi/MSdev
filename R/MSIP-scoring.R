#' @title Calculate MSIP Structural Elucidation Score
#' @description Evaluates how well experimental MS2 spectra can elucidate the structural
#'   and isotopic labeling positions of a metabolite by mapping them to in-silico
#'   theoretical fragments (CFM-ID).
#'
#' @param spectra A data.frame of experimental MS2 peaks with columns \code{mz} and
#'   \code{intensity}. Alternatively, a \code{Spectra} object.
#' @param cfmd A \code{CFM_data} object containing theoretical fragment information,
#'   typically obtained from \code{\link{CFM_predict}} or \code{\link{CFM_annotate_by_predict}}.
#' @param ppm Numeric, m/z tolerance for peak matching in parts per million (default: 10).
#' @param energy Character, which collision energy level to use from CFMD peak_assignment.
#'   One of "energy0", "energy1", "energy2", or "all" (default: "all").
#'
#' @details
#' The function performs three main steps:
#' \enumerate{
#'   \item \strong{Spectra Matching}: Matches experimental m/z values against theoretical
#'     fragment m/z from CFMD using the specified ppm tolerance.
#'   \item \strong{Atomic Position Localization}: Extracts atom indices from matched
#'     fragments and computes the union of all elucidated atomic positions.
#'   \item \strong{Scoring}: Calculates fragment coverage and structural elucidation scores.
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{Score}{The overall Structural Elucidation Score (numeric 0-1), defined as
#'     the ratio of unique atoms covered to total atoms in the molecule.}
#'   \item{Fragment_Coverage}{The fragment matching ratio (numeric 0-1), defined as
#'     the number of matched CFMD fragments divided by total CFMD fragments.}
#'   \item{Matched_Fragments}{A data.frame of the matched subset from \code{cfmd},
#'     including fragment_id, fragment_mz, smiles, and matched atom indices.}
#'   \item{Covered_Atoms}{A unique character vector of the elucidated atomic identifiers
#'     (e.g., "C1", "C2", "O1").}
#'   \item{Total_Atoms}{Total number of atoms in the intact parent molecule.}
#'   \item{Matched_Peaks}{A data.frame of matched experimental peaks with fragment assignments.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example workflow using MSdb
#' # 1. Load compound database and get spectra
#' cpdb <- MSdb::get_CompoundDB()
#' cid_filter <- CompoundDb::CompoundIdFilter("HMDB0000125")
#' sp <- CompoundDb::Spectra(cpdb, filter = cid_filter)
#' cpd <- CompoundDb::compounds(cpdb, filter = cid_filter)
#'
#' # 2. Generate CFM data from SMILES
#' cfmd <- get_CFM_data_from_smiles(cpd$smiles)
#'
#' # 3. Calculate structural elucidation score
#' result <- get_MSIP_structural_elucidation_score(sp, cfmd)
#'
#' # Access results
#' result$Score
#' result$Fragment_Coverage
#' result$Covered_Atoms
#' }
get_MSIP_structural_elucidation_score <- function(spectra,
                                                   cfmd,
                                                   ppm = 10,
                                                   energy = "all") {

  ### Input validation and preprocessing
  {
    if (!is(cfmd, "CFM_data")) {
      stop("'cfmd' must be a CFM_data S4 object.")
    }

    # Convert Spectra object to data.frame if needed
    if (is(spectra, "Spectra")) {
      sp_data <- get_Spectra_data(spectra)
      spectra <- data.frame(
        mz = sp_data$mz,
        intensity = sp_data$intensity
      )
    }

    if (!is.data.frame(spectra)) {
      stop("'spectra' must be a data.frame with columns 'mz' and 'intensity', or a Spectra object.")
    }

    if (!all(c("mz", "intensity") %in% colnames(spectra))) {
      stop("'spectra' must contain columns 'mz' and 'intensity'.")
    }

    if (nrow(spectra) == 0) {
      warning("Empty spectra provided.")
      return(list(
        Score = 0,
        Fragment_Coverage = 0,
        Matched_Fragments = data.frame(),
        Covered_Atoms = character(0),
        Total_Atoms = 0,
        Matched_Peaks = data.frame()
      ))
    }
  }


  ### Step 1: Get CFMD theoretical fragment m/z values
  {
    cfm_peaks <- cfmd@peak_assignment

    # Filter by energy if specified
    if (energy != "all") {
      cfm_peaks <- cfm_peaks %>% dplyr::filter(energy == !!energy)
    }

    # Get unique fragment definitions
    fragment_define <- cfmd@fragment_define

    # Use fragment_define for m/z matching (unique fragments)
    if (nrow(fragment_define) == 0) {
      warning("No fragments found in CFMD fragment_define.")
      return(list(
        Score = 0,
        Fragment_Coverage = 0,
        Matched_Fragments = data.frame(),
        Covered_Atoms = character(0),
        Total_Atoms = 0,
        Matched_Peaks = data.frame()
      ))
    }

    fragment_mz <- fragment_define$fragment_mz
    fragment_ids <- fragment_define$fragment_id
    total_fragments <- length(fragment_ids)
  }


  ### Step 2: Match experimental spectra to theoretical fragments
  {
    exp_mz <- spectra$mz

    # Use match_mz to find closest matches within ppm tolerance
    match_idx <- match_mz(mz1 = exp_mz,
                          mz2 = fragment_mz,
                          mz.ppm = ppm)

    # Create matched peaks data.frame
    matched_peaks <- data.frame(
      exp_mz = exp_mz,
      exp_intensity = spectra$intensity,
      fragment_idx = match_idx,
      fragment_id = fragment_ids[match_idx],
      fragment_mz = fragment_mz[match_idx],
      mz_error_ppm = ifelse(is.na(match_idx), NA,
                            abs(exp_mz - fragment_mz[match_idx]) / exp_mz * 1e6),
      stringsAsFactors = FALSE
    )

    # Filter to only matched peaks
    matched_peaks <- matched_peaks %>% dplyr::filter(!is.na(fragment_idx))

    if (nrow(matched_peaks) == 0) {
      warning("No experimental peaks matched to CFMD fragments within ", ppm, " ppm tolerance.")
      return(list(
        Score = 0,
        Fragment_Coverage = 0,
        Matched_Fragments = data.frame(),
        Covered_Atoms = character(0),
        Total_Atoms = 0,
        Matched_Peaks = data.frame()
      ))
    }

    # Get unique matched fragment IDs
    matched_fragment_ids <- unique(matched_peaks$fragment_id)
    n_matched_fragments <- length(matched_fragment_ids)
  }


  ### Step 3: Extract atom indices from matched fragments
  {
    # Get total atoms from parent molecule (first fragment in fragment_igraph)
    if (length(cfmd@fragment_igraph) == 0) {
      stop("MSIPAtomMap does not contain fragment_igraph. Run MSIPAtomMap_get_igraph() first.")
    }

    parent_igraph <- cfmd@fragment_igraph[[1]]
    all_atoms <- get_sdf_igraph_atom(parent_igraph)
    total_atoms <- length(all_atoms)

    # Extract covered atoms from fragment_atom_map
    covered_atoms_list <- list()

    if (length(cfmd@fragment_atom_map) > 0) {
      # fragment_atom_map is available - use it for precise atom mapping
      fragment_atom_map <- cfmd@fragment_atom_map

      for (frag_id in matched_fragment_ids) {
        if (frag_id %in% names(fragment_atom_map)) {
          atom_matrix <- fragment_atom_map[[frag_id]]
          if (!is.null(atom_matrix) && !identical(atom_matrix, NA)) {
            # Row names represent parent atoms; atoms with any probability > 0 are present
            atoms_present <- rownames(atom_matrix)[rowSums(atom_matrix, na.rm = TRUE) > 0]
            covered_atoms_list[[frag_id]] <- atoms_present
          }
        }
      }
    } else if (length(cfmd@fragment_igraph) > 0) {
      # Fallback: extract atoms directly from fragment igraph objects
      for (frag_id in matched_fragment_ids) {
        if (frag_id %in% names(cfmd@fragment_igraph)) {
          frag_ig <- cfmd@fragment_igraph[[frag_id]]
          frag_atoms <- get_sdf_igraph_atom(frag_ig)
          covered_atoms_list[[frag_id]] <- frag_atoms
        }
      }
    }

    # Compute union of all covered atoms
    covered_atoms <- unique(unlist(covered_atoms_list))
    n_covered_atoms <- length(covered_atoms)
  }


  ### Step 4: Calculate scores
  {
    fragment_coverage <- n_matched_fragments / total_fragments
    structural_elucidation_score <- n_covered_atoms / total_atoms

    # Handle edge cases
    if (total_fragments == 0) fragment_coverage <- 0
    if (total_atoms == 0) structural_elucidation_score <- 0
  }


  ### Step 5: Build matched fragments data.frame
  {
    matched_fragments_df <- fragment_define %>%
      dplyr::filter(fragment_id %in% matched_fragment_ids) %>%
      dplyr::mutate(
        matched_atoms = sapply(fragment_id, function(fid) {
          if (fid %in% names(covered_atoms_list)) {
            paste(covered_atoms_list[[fid]], collapse = ";")
          } else {
            NA_character_
          }
        }),
        n_matched_atoms = sapply(fragment_id, function(fid) {
          if (fid %in% names(covered_atoms_list)) {
            length(covered_atoms_list[[fid]])
          } else {
            0L
          }
        })
      )
  }


  ### Return results
  result <- list(
    Score = structural_elucidation_score,
    Fragment_Coverage = fragment_coverage,
    Matched_Fragments = matched_fragments_df,
    Covered_Atoms = covered_atoms,
    Total_Atoms = total_atoms,
    Matched_Peaks = matched_peaks
  )

  return(result)
}
