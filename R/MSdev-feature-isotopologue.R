#' Natural isotopologue catalog built from \code{MSCC::elem_table}.
#'
#' One row per (isotope, k). \code{prio} orders competing explanations:
#' 1 = 13C, 2 = 15N / 18O / 2H, 3 = 34S / 37Cl / 81Br, 4 = 41K.
#' @noRd
.msdev_isotopologue_catalog <- function(max_atoms = c(S = 3, Cl = 3, Br = 3, K = 1),
                                        polarity = "Positive") {
  spec <- data.frame(
    iso = c("[13]C", "[2]H", "[15]N", "[18]O", "[34]S", "[37]Cl", "[81]Br", "[41]K"),
    parent = c("C", "H", "N", "O", "S", "Cl", "Br", "K"),
    kmax = c(5L, 3L, 3L, 3L, NA, NA, NA, NA),
    hetero = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
    prio = c(1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L),
    stringsAsFactors = FALSE
  )
  het <- spec$hetero
  cap <- as.integer(max_atoms[spec$parent[het]])
  cap[is.na(cap)] <- 0L
  spec$kmax[het] <- cap
  if (!identical(polarity, "Positive")) {
    spec <- spec[spec$parent != "K", , drop = FALSE]
  }
  spec <- spec[spec$kmax > 0L, , drop = FALSE]
  et <- MSCC::elem_table
  spec$mass_diff <- et$Mass_Dif[match(spec$iso, et$element)]
  spec$abundance <- et$abundance[match(spec$iso, et$element)]
  spec$major_mass <- et$mass[match(spec$parent, et$element)]
  spec$ord <- seq_len(nrow(spec))
  spec <- spec[is.finite(spec$mass_diff) & is.finite(spec$abundance) &
                 is.finite(spec$major_mass), , drop = FALSE]
  rows <- lapply(seq_len(nrow(spec)), function(i) {
    k <- seq_len(spec$kmax[[i]])
    data.frame(spec[rep(i, length(k)), , drop = FALSE], k = k,
               mass = spec$mass_diff[[i]] * k, row.names = NULL,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$kmax <- NULL
  out
}

#' Standard ppm of an m/z gap against a theoretical mass difference.
#' @noRd
.msdev_isotopologue_gap_ppm <- function(mz1, mz2, iso_mass) {
  ((mz2 - mz1) - iso_mass) / ((mz1 + mz2) / 2) * 1e6
}

#' Atom count used for the abundance ceiling: CHON scale with m/z,
#' heteroatoms are capped by \code{max_atoms}.
#' @noRd
.msdev_isotopologue_n_atom <- function(mz, parent, major_mass, max_atoms) {
  het <- parent %in% c("S", "Cl", "Br", "K")
  n <- floor(mz / major_mass)
  cap <- as.numeric(max_atoms[parent])
  cap[is.na(cap)] <- 0
  n[het] <- cap[het]
  n
}

#' Largest plausible \code{I_heavy / I_mono} for counts \code{cnt}
#' (named by isotope, e.g. \code{c("[13]C" = 2)}) at \code{mz}.
#'
#' Product over isotopes of \code{choose(n, k) * abundance^k}.
#' @noRd
.msdev_isotopologue_max_rel <- function(mz, cnt, catalog, max_atoms) {
  if (!length(cnt) || !is.finite(mz)) {
    return(NA_real_)
  }
  sp <- catalog[!duplicated(catalog$iso), , drop = FALSE]
  i <- match(names(cnt), sp$iso)
  if (anyNA(i)) {
    return(NA_real_)
  }
  n <- .msdev_isotopologue_n_atom(mz, sp$parent[i], sp$major_mass[i], max_atoms)
  k <- as.numeric(cnt)
  if (any(n < k)) {
    return(0)
  }
  exp(sum(lchoose(n, k) + k * log(sp$abundance[i])))
}

#' Counts as text, e.g. \code{[13]C2C-2[34]S1S-1}, in catalog order.
#' @noRd
.msdev_isotopologue_formula <- function(cnt, catalog) {
  if (!length(cnt)) {
    return(NA_character_)
  }
  sp <- catalog[!duplicated(catalog$iso), , drop = FALSE]
  i <- match(names(cnt), sp$iso)
  o <- order(sp$ord[i])
  k <- as.integer(cnt)[o]
  paste0(names(cnt)[o], k, sp$parent[i][o], "-", k, collapse = "")
}

#' Median \code{into} ratio heavy / light across samples where both are positive.
#' @noRd
.msdev_isotopologue_into_ratio <- function(into, num, den) {
  if (is.null(into) || !length(num)) {
    return(rep(NA_real_, length(num)))
  }
  a <- into[num, , drop = FALSE]
  b <- into[den, , drop = FALSE]
  r <- a / b
  r[!is.finite(r) | r <= 0] <- NA_real_
  out <- apply(r, 1L, stats::median, na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  as.numeric(out)
}

.msdev_isotopologue_into <- function(xcms.xcms) {
  fv <- tryCatch(
    xcms::featureValues(xcms.xcms, value = "into", filled = FALSE),
    error = function(e) NULL
  )
  if (is.null(fv)) {
    fv <- tryCatch(
      xcms::featureValues(xcms.xcms, value = "into", filled = TRUE),
      error = function(e) NULL
    )
  }
  if (is.null(fv)) {
    return(NULL)
  }
  m <- as.matrix(fv)
  storage.mode(m) <- "double"
  m[!is.finite(m) | m <= 0] <- NA_real_
  if (!any(is.finite(m))) {
    return(NULL)
  }
  m
}

#' Candidate edges: similar pairs whose m/z gap fits a catalog row.
#'
#' Rows within \code{3 * ppm} are kept so near misses stay visible;
#' \code{pass} requires mass, \code{min_rel} and intensity.
#' @noRd
.msdev_isotopologue_pairs <- function(sim, fid, mz, into, catalog, sim_threshold,
                                      ppm, max_atoms, intensity_factor_chon,
                                      intensity_factor_hetero, min_rel) {
  empty <- data.frame(
    lo = integer(), hi = integer(), anchor = character(), heavy = character(),
    iso = character(), parent = character(), k = integer(), prio = integer(),
    hetero = logical(), ppm = numeric(), sim = numeric(), ratio = numeric(),
    ceiling = numeric(), factor = numeric(), pass_mass = logical(),
    pass_min_rel = logical(), pass_intensity = logical(), pass = logical(),
    stringsAsFactors = FALSE
  )
  sm <- Matrix::summary(methods::as(sim, "TsparseMatrix"))
  sm <- sm[sm$i < sm$j & sm$x > sim_threshold, , drop = FALSE]
  if (!nrow(sm) || !nrow(catalog)) {
    return(empty)
  }
  i <- match(rownames(sim)[sm$i], fid)
  j <- match(colnames(sim)[sm$j], fid)
  ok <- !is.na(i) & !is.na(j)
  i <- i[ok]
  j <- j[ok]
  x <- sm$x[ok]
  swap <- mz[i] > mz[j]
  lo <- ifelse(swap, j, i)
  hi <- ifelse(swap, i, j)
  gap <- mz[hi] - mz[lo]
  keep <- is.finite(gap) & gap >= 0.5 & gap <= max(catalog$mass) + 0.05
  lo <- lo[keep]
  hi <- hi[keep]
  x <- x[keep]
  if (!length(lo)) {
    return(empty)
  }
  parts <- lapply(seq_len(nrow(catalog)), function(r) {
    p <- .msdev_isotopologue_gap_ppm(mz[lo], mz[hi], catalog$mass[[r]])
    w <- which(abs(p) <= 3 * ppm)
    if (!length(w)) {
      return(NULL)
    }
    data.frame(
      lo = lo[w], hi = hi[w], iso = catalog$iso[[r]],
      parent = catalog$parent[[r]], k = catalog$k[[r]],
      prio = catalog$prio[[r]], hetero = catalog$hetero[[r]],
      ppm = p[w], sim = x[w], stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, parts)
  if (is.null(out) || !nrow(out)) {
    return(empty)
  }
  out$anchor <- fid[out$lo]
  out$heavy <- fid[out$hi]
  out$ceiling <- vapply(seq_len(nrow(out)), function(r) {
    .msdev_isotopologue_max_rel(
      mz[[out$lo[[r]]]], stats::setNames(out$k[[r]], out$iso[[r]]),
      catalog, max_atoms
    )
  }, numeric(1))
  out$ratio <- .msdev_isotopologue_into_ratio(into, out$hi, out$lo)
  out$factor <- ifelse(out$hetero, intensity_factor_hetero, intensity_factor_chon)
  out$pass_mass <- abs(out$ppm) <= ppm
  out$pass_min_rel <- is.finite(out$ceiling) & out$ceiling >= min_rel
  out$pass_intensity <- if (is.null(into)) {
    rep(TRUE, nrow(out))
  } else {
    is.finite(out$ratio) & out$ratio <= out$ceiling * out$factor
  }
  out$pass <- out$pass_mass & out$pass_min_rel & out$pass_intensity
  out <- out[order(out$hi, out$lo, out$prio, abs(out$ppm)), names(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

#' One parent per heavy feature, then walk to the monoisotopic root.
#' @noRd
.msdev_isotopologue_resolve <- function(pairs, sim, fid, mz, into, catalog, ppm,
                                        max_atoms, intensity_factor_chon,
                                        intensity_factor_hetero, min_rel) {
  e <- pairs[pairs$pass, , drop = FALSE]
  empty_links <- data.frame(
    feature_id = character(), mzmed = numeric(), status = character(),
    anchor = character(), parent = character(), formula = character(),
    step_formula = character(), ppm = numeric(), sim = numeric(),
    ratio = numeric(), n_parent = integer(), stringsAsFactors = FALSE
  )
  if (!nrow(e)) {
    return(list(links = empty_links, conflicts = e[0, , drop = FALSE]))
  }
  e <- e[order(e$hi, e$prio, abs(e$ppm), -e$sim), , drop = FALSE]
  best <- e[!duplicated(e$hi), , drop = FALSE]
  linked <- paste(e$lo, e$hi)
  ambiguous <- vapply(best$hi, function(h) {
    pp <- unique(e$lo[e$hi == h])
    if (length(pp) < 2L) {
      return(FALSE)
    }
    cmb <- utils::combn(pp, 2L)
    !all(paste(cmb[1, ], cmb[2, ]) %in% linked | paste(cmb[2, ], cmb[1, ]) %in% linked)
  }, logical(1))
  n_parent <- vapply(best$hi, function(h) length(unique(e$lo[e$hi == h])), integer(1))
  parent_of <- stats::setNames(best$lo, best$hi)
  step_cnt <- lapply(seq_len(nrow(best)), function(r) stats::setNames(best$k[[r]], best$iso[[r]]))
  names(step_cnt) <- as.character(best$hi)
  add_cnt <- function(a, b) {
    for (nm in names(b)) {
      a[nm] <- if (nm %in% names(a)) a[[nm]] + b[[nm]] else b[[nm]]
    }
    a
  }
  sim_of <- function(a, b) {
    v <- as.numeric(sim[fid[a], fid[b]])
    if (is.finite(v) && v > 0) v else NA_real_
  }
  rows <- vector("list", nrow(best))
  for (r in seq_len(nrow(best))) {
    h <- best$hi[[r]]
    cnt <- step_cnt[[as.character(h)]]
    cur <- best$lo[[r]]
    path_sim <- best$sim[[r]]
    while (!is.na(nxt <- parent_of[as.character(cur)])) {
      cnt <- add_cnt(cnt, step_cnt[[as.character(cur)]])
      path_sim <- min(path_sim, e$sim[e$hi == cur & e$lo == nxt][1L])
      cur <- unname(nxt)
    }
    root <- cur
    step_form <- .msdev_isotopologue_formula(step_cnt[[as.character(h)]], catalog)
    status <- "isotopologue"
    anchor <- root
    form <- .msdev_isotopologue_formula(cnt, catalog)
    if (root != best$lo[[r]]) {
      sp <- catalog[match(names(cnt), catalog$iso), , drop = FALSE]
      mass <- sum(as.numeric(cnt) * sp$mass_diff)
      p_root <- .msdev_isotopologue_gap_ppm(mz[[root]], mz[[h]], mass)
      ceil <- .msdev_isotopologue_max_rel(mz[[root]], cnt, catalog, max_atoms)
      fac <- if (any(sp$hetero)) intensity_factor_hetero else intensity_factor_chon
      rat <- .msdev_isotopologue_into_ratio(into, h, root)
      ok <- abs(p_root) <= ppm && is.finite(ceil) && ceil >= min_rel &&
        (is.null(into) || (is.finite(rat) && rat <= ceil * fac))
      if (ok) {
        p_use <- p_root
        s_use <- sim_of(h, root)
        if (is.na(s_use)) s_use <- path_sim
      } else {
        status <- "chain_fallback"
        anchor <- best$lo[[r]]
        form <- step_form
        p_use <- best$ppm[[r]]
        s_use <- best$sim[[r]]
      }
    } else {
      p_use <- best$ppm[[r]]
      s_use <- best$sim[[r]]
    }
    if (ambiguous[[r]]) {
      status <- "ambiguous"
    }
    rows[[r]] <- data.frame(
      feature_id = fid[[h]], mzmed = mz[[h]], status = status,
      anchor = fid[[anchor]], parent = fid[[best$lo[[r]]]],
      formula = form, step_formula = step_form, ppm = p_use, sim = s_use,
      ratio = .msdev_isotopologue_into_ratio(into, h, anchor),
      n_parent = n_parent[[r]], stringsAsFactors = FALSE
    )
  }
  links <- do.call(rbind, rows)
  amb_h <- best$hi[ambiguous]
  conflicts <- e[e$hi %in% amb_h, , drop = FALSE]
  rownames(links) <- NULL
  rownames(conflicts) <- NULL
  list(links = links, conflicts = conflicts)
}


#' @title Natural isotopologue ladders for xcms features
#' @description Finds natural-abundance isotopologues (13C, 2H, 15N, 18O, 34S,
#'   37Cl, 81Br, 41K) among all features of one polarity, without chemical
#'   annotation or feature groups. Only isotopologues (different isotopic
#'   composition, different mass) are detected; isotopomers are not
#'   distinguishable by MS1 mass.
#'
#'   A feature pair is a candidate when its aggregated EIC similarity
#'   (\code{otherData(xcms)$EIC_Similarity_agg}) is above
#'   \code{sim_threshold}. The lighter feature is the parent; the m/z gap must
#'   match \code{k} times an isotope mass difference within \code{ppm}
#'   (standard ppm of the mean m/z). The abundance ceiling is
#'   \code{prod(choose(n, k) * abundance^k)} with \code{n = floor(mz / major mass)}
#'   for C/H/N/O and \code{n = max_atoms} for S/Cl/Br/K. A row is used only if
#'   the ceiling is at least \code{min_rel} (detectable k for that m/z) and the
#'   observed \code{into} ratio is at most ceiling times the factor
#'   (\code{intensity_factor_chon} or \code{intensity_factor_hetero}).
#'   41K is used only in positive mode.
#'
#'   Each heavy feature keeps one parent, ordered by isotope priority
#'   (13C, then 15N/18O/2H, then 34S/37Cl/81Br, then 41K), |ppm| and
#'   similarity. Parents that are not linked to each other make the feature
#'   \code{ambiguous}. Parents are followed to the monoisotopic root; the
#'   cumulative composition is rechecked against the root (mass, \code{min_rel},
#'   intensity), otherwise the direct parent is kept (\code{chain_fallback}).
#'
#'   Constants come from \code{MSCC::elem_table}.
#' @param xcms.xcms xcms object with \code{otherData(xcms)$EIC_Similarity_agg}
#'   (see \code{\link{xcms_feature_EIC_similarity}}).
#' @param polarity \code{"Positive"} or \code{"Negative"}.
#' @param sim_threshold numeric(1). Pairs need similarity strictly above this.
#' @param ppm numeric(1). Tolerance of the m/z gap, ppm of the mean m/z.
#' @param max_atoms Named integer. Maximum S, Cl, Br, K atoms per ion.
#' @param intensity_factor_chon numeric(1). Allowed factor over the ceiling for
#'   C/H/N/O-only compositions.
#' @param intensity_factor_hetero numeric(1). Allowed factor over the ceiling
#'   when S/Cl/Br/K is involved.
#' @param min_rel numeric(1). Minimum theoretical ceiling for a k to be tried.
#' @return Updated \code{xcms.xcms}. \code{featureDefinitions} gains
#'   \code{isotopologue_status} (\code{"monoisotopic"}, \code{"isotopologue"},
#'   \code{"ambiguous"}, \code{"chain_fallback"} or \code{NA}),
#'   \code{isotopologue_anchor} (own id for monoisotopic roots),
#'   \code{isotopologue_formula}, \code{isotopologue_ppm},
#'   \code{isotopologue_sim}. \code{otherData(xcms)$Isotopologue_Global} holds
#'   \code{params}, \code{summary}, \code{pairs}, \code{links},
#'   \code{conflicts}.
#' @export
xcms_feature_natural_isotopologues <- function(xcms.xcms,
                                               polarity = c("Positive", "Negative"),
                                               sim_threshold = 0.9,
                                               ppm = 5,
                                               max_atoms = c(S = 3, Cl = 3, Br = 3, K = 1),
                                               intensity_factor_chon = 6,
                                               intensity_factor_hetero = 1.5,
                                               min_rel = 0.001) {
  polarity <- match.arg(polarity)
  od <- MsExperiment::otherData(xcms.xcms)
  sim <- od$EIC_Similarity_agg
  if (is.null(sim)) {
    stop("otherData(xcms)$EIC_Similarity_agg is missing; ",
         "run MSdev_feature_EIC_similarity() first")
  }
  fd <- as.data.frame(xcms::featureDefinitions(xcms.xcms))
  fid <- if ("feature_id" %in% names(fd)) as.character(fd$feature_id) else rownames(fd)
  mz <- as.numeric(fd$mzmed)
  into <- .msdev_isotopologue_into(xcms.xcms)
  if (is.null(into)) {
    message("Isotopologues ", polarity, ": no feature into; intensity check skipped")
  }
  catalog <- .msdev_isotopologue_catalog(max_atoms, polarity)
  args <- list(
    fid = fid, mz = mz, into = into, catalog = catalog, ppm = ppm,
    max_atoms = max_atoms, intensity_factor_chon = intensity_factor_chon,
    intensity_factor_hetero = intensity_factor_hetero, min_rel = min_rel
  )
  pairs <- do.call(.msdev_isotopologue_pairs,
                   c(list(sim = sim, sim_threshold = sim_threshold), args))
  res <- do.call(.msdev_isotopologue_resolve, c(list(pairs = pairs, sim = sim), args))
  links <- res$links

  n <- nrow(fd)
  status <- rep(NA_character_, n)
  anchor <- rep(NA_character_, n)
  form <- rep(NA_character_, n)
  ppm_v <- rep(NA_real_, n)
  sim_v <- rep(NA_real_, n)
  roots <- setdiff(unique(links$anchor), links$feature_id)
  ri <- match(roots, fid)
  status[ri] <- "monoisotopic"
  anchor[ri] <- roots
  hi <- match(links$feature_id, fid)
  status[hi] <- links$status
  anchor[hi] <- links$anchor
  form[hi] <- links$formula
  ppm_v[hi] <- links$ppm
  sim_v[hi] <- links$sim
  fd$isotopologue_status <- status
  fd$isotopologue_anchor <- anchor
  fd$isotopologue_formula <- form
  fd$isotopologue_ppm <- ppm_v
  fd$isotopologue_sim <- sim_v
  xcms.xcms <- .xcms_featureDefinitions_replace(xcms.xcms, fd)

  near <- pairs$pass_mass
  summary <- data.frame(
    polarity = polarity,
    n_features = n,
    n_mass_pairs = sum(near),
    n_reject_min_rel = sum(near & !pairs$pass_min_rel),
    n_reject_intensity = sum(near & pairs$pass_min_rel & !pairs$pass_intensity),
    n_isotopologue = nrow(links),
    n_monoisotopic = length(roots),
    n_ambiguous = sum(links$status == "ambiguous"),
    n_chain_fallback = sum(links$status == "chain_fallback"),
    stringsAsFactors = FALSE
  )
  od <- MsExperiment::otherData(xcms.xcms)
  od$Isotopologue_Global <- list(
    params = list(
      polarity = polarity, sim_threshold = sim_threshold, ppm = ppm,
      ppm_definition = "standard ppm of mean m/z",
      max_atoms = max_atoms, intensity_factor_chon = intensity_factor_chon,
      intensity_factor_hetero = intensity_factor_hetero, min_rel = min_rel,
      catalog = catalog
    ),
    summary = summary,
    pairs = pairs,
    links = links,
    conflicts = res$conflicts
  )
  MsExperiment::otherData(xcms.xcms) <- od
  message_with_time(
    "Isotopologues ", polarity, ": ", summary$n_isotopologue, " isotopologues on ",
    summary$n_monoisotopic, " monoisotopic roots; ", summary$n_ambiguous,
    " ambiguous, ", summary$n_chain_fallback, " chain fallback; rejected ",
    summary$n_reject_min_rel, " by min_rel, ", summary$n_reject_intensity,
    " by intensity"
  )
  xcms.xcms
}


#' @title Natural isotopologue ladders per polarity
#' @description Per-polarity wrapper around
#'   \code{\link{xcms_feature_natural_isotopologues}}. Run after
#'   \code{\link{MSdev_feature_EIC_similarity}}.
#' @param object MSdev object.
#' @inheritParams xcms_feature_natural_isotopologues
#' @return MSdev object; each polarity's MS1 xcms object carries the
#'   \code{isotopologue_*} columns and \code{otherData(xcms)$Isotopologue_Global}.
#' @export
MSdev_feature_natural_isotopologues <- function(object,
                                                sim_threshold = 0.9,
                                                ppm = 5,
                                                max_atoms = c(S = 3, Cl = 3, Br = 3, K = 1),
                                                intensity_factor_chon = 6,
                                                intensity_factor_hetero = 1.5,
                                                min_rel = 0.001) {
  stopifnot(inherits(object, "MSdev"))
  for (pol in c("Negative", "Positive")) {
    ms1_key <- paste0(pol, "MS1")
    xcms.xcms <- object@xcmsData[[ms1_key]]
    if (is.null(xcms.xcms) || identical(xcms.xcms, NA)) {
      next
    }
    object@xcmsData[[ms1_key]] <- xcms_feature_natural_isotopologues(
      xcms.xcms,
      polarity = pol,
      sim_threshold = sim_threshold,
      ppm = ppm,
      max_atoms = max_atoms,
      intensity_factor_chon = intensity_factor_chon,
      intensity_factor_hetero = intensity_factor_hetero,
      min_rel = min_rel
    )
  }
  object
}
