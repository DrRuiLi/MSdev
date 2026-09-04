#' Feature-group EIC similarity, grouping, and reporting
#'
#' Helpers and exported APIs for pairwise EIC similarity, feature grouping
#' by complete-linkage similarity, and PDF reporting of mirror EIC plots.
#'
#' @name MSdev-feature-group-EIC
NULL

#' Resolve selected_sample to chromatogram column indices
#' @noRd
.resolve_selected_sample <- function(selected_sample, chrom_names, pdata_names = NULL) {
  n <- length(chrom_names)
  if (is.null(selected_sample)) {
    return(seq_len(n))
  }
  if (is.numeric(selected_sample)) {
    idx <- as.integer(selected_sample)
    if (any(is.na(idx)) || any(idx < 1L) || any(idx > n)) {
      stop("'selected_sample' indices must be between 1 and ", n)
    }
    return(idx)
  }
  if (is.character(selected_sample)) {
    idx <- match(selected_sample, chrom_names)
    if (anyNA(idx) && !is.null(pdata_names)) {
      miss <- which(is.na(idx))
      idx2 <- match(selected_sample[miss], pdata_names)
      if (anyNA(idx2)) {
        stop(
          "'selected_sample' not found in chromatogram colnames or sample.name: ",
          paste(selected_sample[miss][is.na(idx2)], collapse = ", ")
        )
      }
      idx[miss] <- idx2
    } else if (anyNA(idx)) {
      stop(
        "'selected_sample' not found in chromatogram colnames: ",
        paste(selected_sample[is.na(idx)], collapse = ", ")
      )
    }
    return(idx)
  }
  stop("'selected_sample' must be NULL, integer index/indices, or sample name(s)")
}


#' Feature index pairs with |rt_i - rt_j| < rt_tol (ordered sliding window)
#' @noRd
.rt_neighbor_pairs <- function(rt, rt_tol) {
  n <- length(rt)
  if (n < 2L) {
    return(matrix(integer(), ncol = 2L))
  }
  ord <- order(rt)
  rt_ord <- rt[ord]
  out_i <- integer(0)
  out_j <- integer(0)
  j <- 2L
  for (i in seq_len(n - 1L)) {
    if (j <= i) {
      j <- i + 1L
    }
    while (j <= n && (rt_ord[j] - rt_ord[i]) < rt_tol) {
      j <- j + 1L
    }
    end <- j - 1L
    if (end >= i + 1L) {
      js <- seq.int(i + 1L, end)
      out_i <- c(out_i, rep.int(ord[i], length(js)))
      out_j <- c(out_j, ord[js])
    }
  }
  cbind(out_i, out_j)
}


#' Split 1:n into n_chunks nearly equal integer index vectors
#' @noRd
.eic_sim_split_job_indices <- function(n, n_chunks) {
  if (n <= 0L) {
    return(list())
  }
  n_chunks <- max(1L, min(as.integer(n_chunks), n))
  base <- n %/% n_chunks
  rem <- n %% n_chunks
  sizes <- rep.int(base, n_chunks)
  if (rem > 0L) {
    sizes[seq_len(rem)] <- sizes[seq_len(rem)] + 1L
  }
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  lapply(seq_len(n_chunks), function(i) seq.int(starts[i], ends[i]))
}


#' Build bplapply payloads: sample-major jobs; only chrom columns used by each chunk
#' @noRd
.eic_sim_make_chunk_payloads <- function(chroms, sample_idx, pairs, n_chunks) {
  ns <- length(sample_idx)
  np <- nrow(pairs)
  n_jobs <- ns * np
  job_groups <- .eic_sim_split_job_indices(n_jobs, n_chunks)
  chrom_by_sample <- lapply(seq_along(sample_idx), function(si) {
    chroms[, sample_idx[si], drop = FALSE]
  })
  lapply(job_groups, function(T) {
    si <- as.integer(((T - 1L) %/% np) + 1L)
    pk <- as.integer(((T - 1L) %% np) + 1L)
    u <- sort(unique(si))
    chroms_sub <- chrom_by_sample[u]
    names(chroms_sub) <- as.character(u)
    list(si = si, pk = pk, chroms = chroms_sub)
  })
}


#' Score one chunk of pooled (sample, pair) jobs; chroms list keyed by si
#' @noRd
.eic_sim_chunk <- function(payload,
                           pairs,
                           ALIGNFUN = MSnbase::alignRt,
                           ALIGNFUNARGS = list(tolerance = 0, method = "closest"),
                           FUN = stats::cor,
                           FUNARGS = list(use = "pairwise.complete.obs")) {
  si <- payload$si
  pk <- payload$pk
  chroms_map <- payload$chroms
  n <- length(si)
  if (n == 0L) {
    return(data.frame(
      si = integer(0), i = integer(0), j = integer(0), score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  i_vec <- pairs[pk, 1L]
  j_vec <- pairs[pk, 2L]
  ## Chromatograms are already zero-filled over their windows on the shared
  ## sample RT grid (see XChromatograms_subset_feature), so RTs align across
  ## features and compareChromatograms can be applied directly.
  scores <- rep(NA_real_, n)
  for (k in seq_len(n)) {
    chroms_col <- chroms_map[[as.character(si[k])]]
    scores[k] <- tryCatch(
      MSnbase::compareChromatograms(
        chroms_col[i_vec[k], 1L],
        chroms_col[j_vec[k], 1L],
        ALIGNFUN = ALIGNFUN,
        ALIGNFUNARGS = ALIGNFUNARGS,
        FUN = FUN,
        FUNARGS = FUNARGS
      ),
      error = function(e) NA_real_
    )
  }
  ok <- is.finite(scores)
  data.frame(
    si = si[ok], i = i_vec[ok], j = j_vec[ok], score = scores[ok],
    stringsAsFactors = FALSE
  )
}


#' Build diagonal-only or scored sparse similarity matrix for one sample
#' @noRd
.eic_sim_sparse_from_edges <- function(i, j, score, nft, fids) {
  if (nft < 1L) {
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(0L, 0L), dimnames = list(fids, fids)
    ))
  }
  if (length(score) == 0L) {
    return(Matrix::sparseMatrix(
      i = seq_len(nft), j = seq_len(nft), x = rep(1, nft),
      dims = c(nft, nft), dimnames = list(fids, fids)
    ))
  }
  Matrix::sparseMatrix(
    i = c(i, j, seq_len(nft)),
    j = c(j, i, seq_len(nft)),
    x = c(score, score, rep(1, nft)),
    dims = c(nft, nft),
    dimnames = list(fids, fids)
  )
}


#' Densify sparse similarity; fill absent (non-neighbor) entries
#'
#' @param sp sparse similarity matrix.
#' @param absent numeric(1). Value for pairs with no stored edge (outside
#'   \code{rt_tol}). Use \code{0} (default) to treat non-overlap as dissimilar,
#'   or \code{NA_real_} to leave them unknown for
#'   \code{MsFeatures::groupSimilarityMatrix}.
#' @noRd
.sparse_sim_to_dense <- function(sp, absent = 0) {
  if (length(absent) != 1L || !(is.na(absent) || is.numeric(absent))) {
    stop("'absent' must be a single numeric or NA")
  }
  n <- nrow(sp)
  mat <- matrix(as.numeric(absent), nrow = n, ncol = n, dimnames = dimnames(sp))
  s <- Matrix::summary(sp)
  if (nrow(s) > 0L) {
    mat[cbind(s$i, s$j)] <- s$x
  }
  if (n > 0L) {
    diag(mat) <- 1
  }
  mat
}

#' Backward-compatible NA densify (calls \code{.sparse_sim_to_dense})
#' @noRd
.sparse_sim_to_dense_na <- function(sp) {
  .sparse_sim_to_dense(sp, absent = NA_real_)
}


#' Aggregate a list of sparse similarity matrices (75%-quantile by default)
#' @noRd
.aggregate_sparse_sim_list <- function(sim_list, prob = 0.75) {
  sp0 <- sim_list[[1L]]
  nft <- nrow(sp0)
  fids <- rownames(sp0)
  parts <- lapply(sim_list, function(sp) {
    s <- Matrix::summary(sp)
    s[s$i != s$j, , drop = FALSE]
  })
  all <- do.call(rbind, parts)
  if (is.null(all) || nrow(all) == 0L) {
    return(Matrix::sparseMatrix(
      i = seq_len(nft), j = seq_len(nft), x = rep(1, nft),
      dims = c(nft, nft), dimnames = list(fids, fids)
    ))
  }
  agg <- stats::aggregate(
    x ~ i + j,
    data = all,
    FUN = function(z) {
      as.numeric(stats::quantile(z, probs = prob, na.rm = TRUE))
    }
  )
  finite <- is.finite(agg$x)
  agg <- agg[finite, , drop = FALSE]
  Matrix::sparseMatrix(
    i = c(agg$i, seq_len(nft)),
    j = c(agg$j, seq_len(nft)),
    x = c(agg$x, rep(1, nft)),
    dims = c(nft, nft),
    dimnames = list(fids, fids)
  )
}


#' @title Pairwise EIC similarity for xcms features
#' @description Compute per-sample feature-by-feature EIC similarity matrices
#'   for pairs with \code{|rtmed_i - rtmed_j| < rt_tol}. Used by
#'   \code{\link{MSdev_group_feature_EIC}}. Results are stored as symmetric
#'   \code{Matrix} sparse matrices; absent entries are non-neighbor pairs
#'   (filled when densified for grouping; see \code{absent_sim} in
#'   \code{\link{xcms_group_feature_EIC}}). Pairwise scores are pooled across
#'   selected samples and split into \code{n_chunks} BiocParallel jobs (each
#'   chunk ships only the chromatogram columns it needs). EICs are expected to
#'   be zero-filled over their windows on the shared sample RT grid (via
#'   \code{\link{XChromatograms_subset_feature}}), so correlation reflects the
#'   full peak window rather than only the non-\code{NA} overlap.
#' @param xcms.xcms XCMSnExp / XcmsExperiment with feature definitions.
#' @param chroms Chromatograms (e.g. from \code{\link{get_xcms_feature_chromatogram}})
#'   with feature rownames matching \code{featureDefinitions}.
#' @param rt_tol numeric(1). Maximum absolute RT difference (seconds) for which
#'   EIC similarity is computed. Default 5.
#' @param expandRt numeric(1). Seconds added on each side of each feature's
#'   \code{peakRtMin}/\code{peakRtMax} window when cropping EICs via
#'   \code{\link{XChromatograms_subset_feature}} before correlation (requires
#'   \code{peakRtMin}/\code{peakRtMax} in \code{chroms@featureDefinitions}).
#'   Default \code{2}.
#' @param min_width numeric(1). Minimum RT window width (seconds) after
#'   \code{expandRt}; shorter windows are padded equally on both sides.
#'   Default \code{20}.
#' @param selected_sample NULL, integer index/indices, or sample name(s)
#'   (\code{sample.name} / chromatogram colnames). NULL uses all samples.
#' @param n_chunks NULL or positive integer. Number of BiocParallel chunks into
#'   which \code{n_samples * n_rt_pairs} similarity jobs are split. When
#'   \code{NULL} (default), uses \code{min(bpnworkers(BPPARAM), n_jobs)} for
#'   backward compatibility. Set explicitly (e.g. \code{200L}) on large batch
#'   projects to keep Snow payloads small; \code{BPPARAM} workers then control
#'   parallelism only.
#' @param BPPARAM BiocParallel backend. Default \code{SerialParam()}.
#' @return Named list of feature-by-feature \code{dgCMatrix} similarity matrices
#'   (one per selected sample).
#' @importFrom Matrix sparseMatrix summary
#' @export
get_xcms_feature_EIC_similarity <- function(xcms.xcms,
                                            chroms,
                                            rt_tol = 5,
                                            expandRt = 2,
                                            min_width = 20,
                                            selected_sample = NULL,
                                            n_chunks = NULL,
                                            BPPARAM = BiocParallel::SerialParam()) {
  if (!is.numeric(rt_tol) || length(rt_tol) != 1L || rt_tol <= 0) {
    stop("'rt_tol' must be a positive numeric(1)")
  }

  chroms <- XChromatograms_subset_feature(chroms,
                                          expandRt = expandRt,
                                          min_width = min_width)

  fdf <- as.data.frame(xcms::featureDefinitions(xcms.xcms))
  fids <- rownames(chroms)
  if (is.null(fids) || anyNA(fids)) {
    stop("Chromatograms lack feature rownames")
  }
  if ("feature_id" %in% colnames(fdf)) {
    rt <- fdf$rtmed[match(fids, fdf$feature_id)]
  } else {
    rt <- fdf$rtmed[match(fids, rownames(fdf))]
  }
  if (anyNA(rt)) {
    stop("Could not map chromatogram features to rtmed")
  }

  chrom_names <- colnames(chroms)
  if (is.null(chrom_names)) {
    chrom_names <- paste0("sample_", seq_len(ncol(chroms)))
    colnames(chroms) <- chrom_names
  }
  pdata_names <- tryCatch(
    as.character(Biobase::pData(xcms.xcms)$sample.name),
    error = function(e) NULL
  )
  sample_idx <- .resolve_selected_sample(selected_sample, chrom_names, pdata_names)
  sample_names <- chrom_names[sample_idx]

  pairs <- .rt_neighbor_pairs(rt, rt_tol)
  ns <- length(sample_idx)
  np <- nrow(pairs)
  n_jobs <- ns * np
  n_workers <- max(1L, as.integer(BiocParallel::bpnworkers(BPPARAM)))
  n_chunks <- .resolve_eic_n_chunks(n_jobs, n_workers, n_chunks)
  message(
    "  features=", length(fids),
    " rt_neighbor_pairs=", np,
    " samples=", ns,
    " jobs=", n_jobs,
    " chunks=", if (n_jobs == 0L) 0L else n_chunks,
    " workers=", n_workers,
    " (rt_tol=", rt_tol, ")"
  )

  nft <- length(fids)
  ALIGNFUNARGS <- list(tolerance = 0, method = "closest")
  FUNARGS <- list(use = "pairwise.complete.obs")

  if (n_jobs == 0L) {
    sim_list <- vector("list", ns)
    names(sim_list) <- sample_names
    for (si in seq_len(ns)) {
      sim_list[[si]] <- .eic_sim_sparse_from_edges(
        i = integer(0), j = integer(0), score = numeric(0),
        nft = nft, fids = fids
      )
    }
    return(sim_list)
  }

  payloads <- .eic_sim_make_chunk_payloads(
    chroms = chroms,
    sample_idx = sample_idx,
    pairs = pairs,
    n_chunks = n_chunks
  )
  res_chunks <- BiocParallel::bplapply(
    payloads,
    FUN = .eic_sim_chunk,
    pairs = pairs,
    ALIGNFUNARGS = ALIGNFUNARGS,
    FUNARGS = FUNARGS,
    BPPARAM = BPPARAM
  )
  edges <- do.call(rbind, res_chunks)

  sim_list <- vector("list", ns)
  names(sim_list) <- sample_names
  for (si in seq_len(ns)) {
    if (is.null(edges) || nrow(edges) == 0L) {
      sub_i <- integer(0)
      sub_j <- integer(0)
      sub_x <- numeric(0)
    } else {
      keep <- edges$si == si
      sub_i <- edges$i[keep]
      sub_j <- edges$j[keep]
      sub_x <- edges$score[keep]
    }
    sim_list[[si]] <- .eic_sim_sparse_from_edges(
      i = sub_i, j = sub_j, score = sub_x,
      nft = nft, fids = fids
    )
  }
  sim_list
}


#' Resolve BiocParallel chunk count for EIC similarity jobs
#' @noRd
.resolve_eic_n_chunks <- function(n_jobs, n_workers, n_chunks) {
  if (is.null(n_chunks)) {
    return(max(1L, min(n_workers, max(n_jobs, 1L))))
  }
  n_chunks <- as.integer(n_chunks)[[1L]]
  if (!is.finite(n_chunks) || n_chunks < 1L) {
    stop("'n_chunks' must be NULL or a positive integer")
  }
  max(1L, min(n_chunks, max(n_jobs, 1L)))
}


#' Group features by complete-linkage similarity threshold
#'
#' Reimplementation of \code{MsFeatures::groupSimilarityMatrix} with a fixed
#' group-ID lookup when joining an existing group (named keys instead of
#' integer positional indexing). Join rule is complete-linkage: a feature may
#' join a group only if similarity to \emph{all} current members is finite and
#' \code{>= threshold}.
#'
#' Available arguments: \code{x}, \code{threshold}, \code{full}, \code{...}.
#'
#' @param x numeric square similarity matrix (symmetric when \code{full = TRUE}).
#' @param threshold numeric(1). Minimum similarity to join / stay linked.
#'   Default \code{0.9}.
#' @param full logical(1). If \code{FALSE}, only the upper triangle is used.
#' @param ... Ignored; kept for call compatibility with MsFeatures.
#' @return Integer vector of group IDs, length \code{nrow(x)}.
#' @seealso \code{\link{groupSimilarityMatrix_hclustAverage}},
#'   \code{\link{xcms_group_feature_EIC}}
#' @export
groupSimilarityMatrix_completeLinkage <- function(x,
                                                  threshold = 0.9,
                                                  full = TRUE,
                                                  ...) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' should be a matrix")
  }
  x <- as.matrix(x)
  nr <- nrow(x)
  if (nr != ncol(x)) {
    stop("'x' should be a symmetric matrix")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L) {
    stop("'threshold' must be numeric(1)")
  }
  if (!is.logical(full) || length(full) != 1L) {
    stop("'full' must be logical(1)")
  }

  if (!full) {
    x[lower.tri(x)] <- NA
  }
  res <- rep(NA_integer_, nr)
  sl <- seq_len(nr)
  x[cbind(sl, sl)] <- NA

  idx_pairs <- which(x >= threshold, arr.ind = TRUE)
  if (nrow(idx_pairs) == 0L) {
    return(seq_len(nr))
  }
  idx_pairs <- idx_pairs[order(x[idx_pairs], decreasing = TRUE), , drop = FALSE]
  grp_id <- 1L

  .passes_complete <- function(cors) {
    !(any(is.na(cors)) || any(cors < threshold))
  }

  for (i in seq_len(nrow(idx_pairs))) {
    got_grp <- res[idx_pairs[i, ]]
    nas <- is.na(got_grp)
    if (any(nas)) {
      if (sum(nas) == 2L) {
        grps <- unique(res[!is.na(res)])
        mean_cor_to_grp <- setNames(
          rep(NA_real_, length(grps)),
          as.character(grps)
        )
        idx <- idx_pairs[i, ]
        for (grp in grps) {
          idx_grp <- which(res == grp)
          cor_to_grp <- x[idx, idx_grp]
          if (full) {
            cor_to_grp <- c(cor_to_grp, x[idx_grp, idx])
          }
          if (.passes_complete(cor_to_grp)) {
            mean_cor_to_grp[[as.character(grp)]] <- mean(cor_to_grp)
          }
        }
        mean_cor_to_grp <- mean_cor_to_grp[is.finite(mean_cor_to_grp)]
        if (length(mean_cor_to_grp)) {
          best <- names(sort(mean_cor_to_grp, decreasing = TRUE))[[1L]]
          res[idx] <- as.integer(best)
        } else {
          res[idx] <- grp_id
          grp_id <- grp_id + 1L
        }
      } else {
        idx <- idx_pairs[i, nas]
        idx_grp <- which(res == got_grp[!nas])
        cor_to_grp <- x[idx, idx_grp]
        if (full) {
          cor_to_grp <- c(cor_to_grp, x[idx_grp, idx])
        }
        if (.passes_complete(cor_to_grp)) {
          res[idx] <- got_grp[!nas]
        }
      }
    } else if (length(unique(got_grp)) > 1L) {
      grp_1 <- which(res == got_grp[1])
      grp_2 <- which(res == got_grp[2])
      cor_1_1 <- x[idx_pairs[i, 1], grp_1]
      cor_1_2 <- x[idx_pairs[i, 1], grp_2]
      cor_2_1 <- x[idx_pairs[i, 2], grp_1]
      cor_2_2 <- x[idx_pairs[i, 2], grp_2]
      if (full) {
        cor_1_1 <- c(cor_1_1, x[grp_1, idx_pairs[i, 1]])
        cor_1_2 <- c(cor_1_2, x[grp_2, idx_pairs[i, 1]])
        cor_2_1 <- c(cor_2_1, x[grp_1, idx_pairs[i, 2]])
        cor_2_2 <- c(cor_2_2, x[grp_2, idx_pairs[i, 2]])
      }
      mcor_1_1 <- mean(cor_1_1, na.rm = TRUE)
      mcor_1_2 <- mean(cor_1_2)
      mcor_2_1 <- mean(cor_2_1)
      mcor_2_2 <- mean(cor_2_2, na.rm = TRUE)
      if (is.finite(mcor_1_2) && is.finite(mcor_1_1) &&
          !any(cor_1_2 < threshold) && mcor_1_2 >= mcor_1_1) {
        res[idx_pairs[i, 1]] <- got_grp[2]
      } else if (is.finite(mcor_2_1) && is.finite(mcor_2_2) &&
                 !any(cor_2_1 < threshold) && mcor_2_1 >= mcor_2_2) {
        res[idx_pairs[i, 2]] <- got_grp[1]
      }
    }
  }

  nas <- is.na(res)
  if (any(nas)) {
    res[nas] <- seq.int(grp_id, length.out = sum(nas))
  }
  res
}


#' Group features by average-linkage hierarchical clustering on similarity
#'
#' Converts a similarity matrix to distance \code{1 - similarity}, runs
#' \code{stats::hclust} with \code{method = "average"} (UPGMA), and cuts the
#' tree at height \code{1 - threshold} so clusters keep average similarity
#' at least \code{threshold}. Non-finite similarities are set to 0 before
#' clustering.
#'
#' Available arguments: \code{x}, \code{threshold}.
#'
#' @param x numeric square similarity matrix.
#' @param threshold numeric(1). Minimum average similarity within a cluster
#'   (cut height \code{1 - threshold}). Default \code{0.9}.
#' @return Integer vector of group IDs, length \code{nrow(x)}.
#' @seealso \code{\link{groupSimilarityMatrix_completeLinkage}},
#'   \code{\link{xcms_group_feature_EIC}}
#' @export
groupSimilarityMatrix_hclustAverage <- function(x, threshold = 0.9) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' should be a matrix")
  }
  x <- as.matrix(x)
  nr <- nrow(x)
  if (nr != ncol(x)) {
    stop("'x' should be a square matrix")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L) {
    stop("'threshold' must be numeric(1)")
  }
  if (nr == 0L) {
    return(integer(0L))
  }
  if (nr == 1L) {
    return(1L)
  }

  sl <- seq_len(nr)
  x[cbind(sl, sl)] <- 1
  x[!is.finite(x)] <- 0
  # Keep matrix dims: pmax(scalar, matrix) drops dim and breaks as.dist.
  dx <- 1 - x
  dx[dx < 0] <- 0
  d <- stats::as.dist(dx)
  hc <- stats::hclust(d, method = "average")
  as.integer(stats::cutree(hc, h = 1 - threshold))
}


#' @title Group xcms features by EIC similarity within RT tolerance
#' @description Compare extracted-ion chromatogram shapes for feature pairs with
#'   \code{|rtmed_i - rtmed_j| < rt_tol} via \code{\link{get_xcms_feature_EIC_similarity}},
#'   aggregate across selected samples (75\% quantile), assign \code{feature_group}
#'   labels with \code{method}, and optionally store per-sample sparse similarity
#'   matrices in \code{otherData(xcms)$EIC_Similarity}.
#' @param xcms.xcms XcmsExperiment / MsExperiment with feature definitions and
#'   an \code{otherData} slot. (all methods)
#' @param chroms Chromatograms with feature rownames matching
#'   \code{featureDefinitions}. (all methods)
#' @param rt_tol numeric(1). Maximum absolute RT difference (seconds) for which
#'   EIC similarity is computed. Default 5. (all methods)
#' @param threshold numeric(1). Similarity cut-off for the chosen \code{method}
#'   (see Details). Default 0.5. (all methods)
#' @param expandRt numeric(1). Seconds added on each side of each feature's
#'   \code{peakRtMin}/\code{peakRtMax} window when cropping EICs via
#'   \code{\link{XChromatograms_subset_feature}} before correlation.
#'   Default \code{2}. (all methods)
#' @param min_width numeric(1). Minimum RT window width (seconds) after
#'   \code{expandRt}; shorter windows are padded equally on both sides.
#'   Default \code{20}. (all methods)
#' @param selected_sample NULL, integer index/indices, or sample name(s).
#'   NULL uses all samples. (all methods)
#' @param keep_Similarity_Matrix logical(1). If TRUE (default), store the named
#'   list of per-sample \code{dgCMatrix} similarity matrices in
#'   \code{otherData(xcms)$EIC_Similarity}. If FALSE, matrices are discarded after
#'   grouping. (all methods)
#' @param absent_sim numeric(1). Fill value for pairs outside \code{rt_tol}
#'   when densifying the sparse similarity matrix before grouping. Default
#'   \code{0} (non-overlap treated as dissimilar). Interpretation depends on
#'   \code{method} (see Details). (all methods)
#' @param method character(1). Grouping algorithm on the dense similarity
#'   matrix. One of \code{"complete_linkage"}, \code{"hclust_average"}.
#'   Default \code{"complete_linkage"}. See Details for arguments used by each
#'   method.
#' @param n_chunks NULL or positive integer. Passed to
#'   \code{\link{get_xcms_feature_EIC_similarity}}. (all methods)
#' @param BPPARAM BiocParallel backend passed to
#'   \code{\link{get_xcms_feature_EIC_similarity}}. Default \code{SerialParam()}.
#'   (all methods)
#' @details
#' \strong{Shared arguments} (used before grouping for every \code{method}):
#' \code{xcms.xcms}, \code{chroms}, \code{rt_tol}, \code{expandRt},
#' \code{min_width}, \code{selected_sample}, \code{keep_Similarity_Matrix},
#' \code{absent_sim}, \code{n_chunks}, \code{BPPARAM}.
#'
#' \strong{Arguments by \code{method}} (how grouping uses densified similarity):
#' \describe{
#'   \item{\code{complete_linkage}}{
#'     Calls \code{\link{groupSimilarityMatrix_completeLinkage}(x, threshold)}.
#'     Available args: \code{threshold} (minimum finite similarity to every
#'     current group member to join); \code{absent_sim} (\code{NA} leaves
#'     pairs unknown and blocks joins that require those pairs;
#'     numeric fill is used as the similarity value).
#'     Helper-only args not exposed here: \code{full}, \code{...}.
#'   }
#'   \item{\code{hclust_average}}{
#'     Calls \code{\link{groupSimilarityMatrix_hclustAverage}(x, threshold)}.
#'     Available args: \code{threshold} (tree cut at height
#'     \code{1 - threshold}, i.e. average similarity \eqn{\ge} threshold);
#'     \code{absent_sim} (non-finite entries are set to similarity 0 before
#'     converting to distance \code{1 - sim}).
#'   }
#' }
#' @return Updated \code{xcms.xcms} with \code{featureGroups} set and a
#'   \code{feature_group_rt} column added to \code{featureDefinitions} (the
#'   median \code{rtmed} of each group's member features); when
#'   \code{keep_Similarity_Matrix} is TRUE, also
#'   \code{otherData(xcms)$EIC_Similarity}.
#' @export
xcms_group_feature_EIC <- function(xcms.xcms,
                                   chroms,
                                   rt_tol = 5,
                                   threshold = 0.5,
                                   expandRt = 2,
                                   min_width = 20,
                                   selected_sample = NULL,
                                   keep_Similarity_Matrix = TRUE,
                                   absent_sim = 0,
                                   method = c("complete_linkage", "hclust_average"),
                                   n_chunks = NULL,
                                   BPPARAM = BiocParallel::SerialParam()) {
  method <- match.arg(method)
  if (!is.numeric(rt_tol) || length(rt_tol) != 1L || rt_tol <= 0) {
    stop("'rt_tol' must be a positive numeric(1)")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L) {
    stop("'threshold' must be numeric(1)")
  }
  if (!is.logical(keep_Similarity_Matrix) || length(keep_Similarity_Matrix) != 1L) {
    stop("'keep_Similarity_Matrix' must be logical(1)")
  }
  if (length(absent_sim) != 1L || !(is.na(absent_sim) || is.numeric(absent_sim))) {
    stop("'absent_sim' must be a single numeric or NA")
  }

  sim_list <- get_xcms_feature_EIC_similarity(
    xcms.xcms = xcms.xcms,
    chroms = chroms,
    rt_tol = rt_tol,
    expandRt = expandRt,
    min_width = min_width,
    selected_sample = selected_sample,
    n_chunks = n_chunks,
    BPPARAM = BPPARAM
  )

  fids <- rownames(sim_list[[1L]])
  fdf <- as.data.frame(xcms::featureDefinitions(xcms.xcms))

  ns <- length(sim_list)
  sim_sp <- if (ns == 1L) {
    sim_list[[1L]]
  } else {
    .aggregate_sparse_sim_list(sim_list, prob = 0.75)
  }
  sim2d <- .sparse_sim_to_dense(sim_sp, absent = absent_sim)

  grp <- switch(
    method,
    complete_linkage = groupSimilarityMatrix_completeLinkage(
      sim2d, threshold = threshold
    ),
    hclust_average = groupSimilarityMatrix_hclustAverage(
      sim2d, threshold = threshold
    )
  )
  f_new <- paste0("FG.", MsFeatures:::.format_id(grp))
  names(f_new) <- fids

  if ("feature_id" %in% colnames(fdf)) {
    lab <- f_new[as.character(fdf$feature_id)]
  } else {
    lab <- f_new[rownames(fdf)]
  }
  if (anyNA(lab)) {
    lab[is.na(lab)] <- paste0("FG.", MsFeatures:::.format_id(seq_len(sum(is.na(lab)))))
  }
  xcms::featureGroups(xcms.xcms) <- as.character(lab)

  fd <- as.data.frame(xcms::featureDefinitions(xcms.xcms))
  fd$feature_group_rt <- stats::ave(
    as.numeric(fd$rtmed),
    as.character(lab),
    FUN = function(x) stats::median(x, na.rm = TRUE)
  )
  xcms.xcms <- .xcms_featureDefinitions_replace(xcms.xcms, fd)

  if (isTRUE(keep_Similarity_Matrix)) {
    od <- MsExperiment::otherData(xcms.xcms)
    od$EIC_Similarity <- sim_list
    MsExperiment::otherData(xcms.xcms) <- od
  }

  message("  feature groups: ", length(unique(lab)))
  xcms.xcms
}


#' @title Group features by EIC similarity within RT tolerance
#' @description Per-polarity wrapper around \code{\link{xcms_group_feature_EIC}}.
#'   Reuses chromatograms from \code{\link{MSdev_get_feature_chrom}}, updates MS1
#'   \code{feature_group} labels, and optionally stores per-sample sparse
#'   similarity matrices in \code{otherData(xcms)$EIC_Similarity}.
#' @param object MSdev object
#' @param rt_tol numeric(1). Maximum absolute RT difference (seconds) for which
#'   EIC similarity is computed. Default 5.
#' @param threshold numeric(1). Similarity cut-off for the chosen
#'   \code{method} (see \code{\link{xcms_group_feature_EIC}}). Default 0.5.
#' @param expandRt numeric(1). Seconds added on each side of each feature's
#'   \code{peakRtMin}/\code{peakRtMax} window when cropping EICs via
#'   \code{\link{XChromatograms_subset_feature}} before correlation.
#'   Default \code{2}.
#' @param min_width numeric(1). Minimum RT window width (seconds) after
#'   \code{expandRt}; shorter windows are padded equally on both sides.
#'   Default \code{20}.
#' @param selected_sample NULL, integer index/indices, or sample name(s)
#'   (\code{sample.name} / chromatogram colnames). NULL uses all samples.
#' @param forceExtractChrom logical(1). If TRUE, (re)extract chromatograms via
#'   \code{MSdev_get_feature_chrom} even if already stored.
#' @param keep_Similarity_Matrix logical(1). Passed to
#'   \code{xcms_group_feature_EIC}. Default TRUE.
#' @param absent_sim numeric(1). Passed to \code{xcms_group_feature_EIC}.
#'   Default \code{0}; set \code{NA} for unknown absents.
#' @param method character(1). Grouping method passed to
#'   \code{\link{xcms_group_feature_EIC}}: \code{"complete_linkage"} or
#'   \code{"hclust_average"}. Default \code{"complete_linkage"}.
#' @param n_chunks NULL or positive integer. Passed to
#'   \code{\link{xcms_group_feature_EIC}} / \code{\link{get_xcms_feature_EIC_similarity}}.
#'   Use on large batch projects to avoid oversized Snow payloads (e.g. \code{200L}).
#' @param BPPARAM BiocParallel backend for chromatogram extraction and EIC
#'   similarity scoring.
#' @return MSdev object with updated MS1 \code{feature_group} labels; when
#'   \code{keep_Similarity_Matrix} is TRUE, also
#'   \code{otherData(xcms)$EIC_Similarity} on each polarity's MS1 object.
#' @export
MSdev_group_feature_EIC <- function(object,
                                    rt_tol = 5,
                                    threshold = 0.5,
                                    expandRt = 2,
                                    min_width = 20,
                                    selected_sample = NULL,
                                    forceExtractChrom = FALSE,
                                    keep_Similarity_Matrix = TRUE,
                                    absent_sim = 0,
                                    method = c("complete_linkage", "hclust_average"),
                                    n_chunks = NULL,
                                    BPPARAM = SnowParam(
                                      workers = max(1L, floor(parallel::detectCores() / 2)),
                                      progressbar = TRUE)) {
  stopifnot(inherits(object, "MSdev"))
  method <- match.arg(method)
  if (!is.numeric(rt_tol) || length(rt_tol) != 1L || rt_tol <= 0) {
    stop("'rt_tol' must be a positive numeric(1)")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L) {
    stop("'threshold' must be numeric(1)")
  }
  if (length(absent_sim) != 1L || !(is.na(absent_sim) || is.numeric(absent_sim))) {
    stop("'absent_sim' must be a single numeric or NA")
  }

  polarities <- c("Negative", "Positive")
  need_chrom <- forceExtractChrom
  if (!need_chrom) {
    for (pol in polarities) {
      ms1_key <- paste0(pol, "MS1")
      chrom_key <- paste0(pol, "_Chromatograms")
      xcms_pol <- object@xcmsData[[ms1_key]]
      if (!is.null(xcms_pol) && !identical(xcms_pol, NA) &&
          is.null(object@xcmsData[[chrom_key]])) {
        need_chrom <- TRUE
        break
      }
    }
  }
  if (need_chrom) {
    message_with_time("Extracting feature chromatograms via MSdev_get_feature_chrom")
    object <- MSdev_get_feature_chrom(object, BPPARAM = BPPARAM)
  }

  for (pol in polarities) {
    ms1_key <- paste0(pol, "MS1")
    chrom_key <- paste0(pol, "_Chromatograms")
    xcms.xcms <- object@xcmsData[[ms1_key]]
    if (is.null(xcms.xcms) || identical(xcms.xcms, NA)) {
      next
    }
    if (is.null(object@xcmsData[[chrom_key]])) {
      warning("No chromatograms for ", pol, "; skip.")
      next
    }

    message_with_time("EIC similarity grouping: ", pol)
    xcms.xcms <- xcms_group_feature_EIC(
      xcms.xcms = xcms.xcms,
      chroms = onDiskData_retrieve(object@xcmsData[[chrom_key]]),
      rt_tol = rt_tol,
      threshold = threshold,
      expandRt = expandRt,
      min_width = min_width,
      selected_sample = selected_sample,
      keep_Similarity_Matrix = keep_Similarity_Matrix,
      absent_sim = absent_sim,
      method = method,
      n_chunks = n_chunks,
      BPPARAM = BPPARAM
    )
    object@xcmsData[[ms1_key]] <- xcms.xcms
  }

  object
}


#' @title Report feature-group EIC mirror plots as PDFs
#' @description For each polarity with MS1 feature groups, generate
#'   \code{\link{plot_xcms_feature_group_EIC_comparasion}} for every group with
#'   at least two features and append pages into
#'   \code{Feature_group_EIC_Positive.pdf} /
#'   \code{Feature_group_EIC_Negative.pdf} under \code{projectDir}.
#'   Reuses stored chromatograms when available
#'   (\code{Positive_Chromatograms} / \code{Negative_Chromatograms}).
#' @param object An \code{MSdev} object with \code{featureGroups} on MS1 xcms
#'   objects (e.g. after \code{\link{MSdev_group_feature_EIC}}).
#' @param expandRt numeric(1). Seconds added on each side of each feature's
#'   peak RT window, passed to
#'   \code{plot_xcms_feature_group_EIC_comparasion} (default \code{2}).
#' @param min_width numeric(1). Minimum RT window width (seconds), passed to
#'   \code{plot_xcms_feature_group_EIC_comparasion} (default \code{20}).
#' @param max_features Maximum features per group plot; \code{NULL} keeps all
#'   (default).
#' @param sample_index Sample index into stored \code{EIC_Similarity}
#'   (default \code{1L}).
#' @param min_features Minimum group size to plot (default \code{2L}).
#' @param width,height PDF page size in inches. If \code{NULL}, scaled from
#'   the number of features in each group.
#' @return Invisible character vector of written PDF path(s).
#' @export
Report_MSdev_feature_group_EIC <- function(object,
                                           expandRt = 2,
                                           min_width = 20,
                                           max_features = NULL,
                                           sample_index = 1L,
                                           min_features = 2L,
                                           width = NULL,
                                           height = NULL) {
  stopifnot(inherits(object, "MSdev"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for Report_MSdev_feature_group_EIC()")
  }
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("Package 'xcms' is required for Report_MSdev_feature_group_EIC()")
  }

  project_dir <- object@projectInfo$projectDir
  if (is.null(project_dir) || !nzchar(project_dir)) {
    stop("`object@projectInfo$projectDir` is missing")
  }
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)

  polarities <- c("Positive", "Negative")
  out_files <- character(0)

  for (pol in polarities) {
    ms1_key <- paste0(pol, "MS1")
    chrom_key <- paste0(pol, "_Chromatograms")
    xcms.xcms <- object@xcmsData[[ms1_key]]
    if (is.null(xcms.xcms) || identical(xcms.xcms, NA)) {
      next
    }

    fg_all <- as.character(xcms::featureGroups(xcms.xcms))
    fg_tab <- table(fg_all[!is.na(fg_all) & nzchar(fg_all)])
    fg_keep <- names(fg_tab)[as.integer(fg_tab) >= as.integer(min_features)]
    if (!length(fg_keep)) {
      message("No feature groups with >= ", min_features,
              " features for ", pol, "; skip.")
      next
    }
    ## Prefer stable order by group size (largest first), then name
    fg_keep <- fg_keep[order(-as.integer(fg_tab[fg_keep]), fg_keep)]

    chroms <- NULL
    chrom_store <- object@xcmsData[[chrom_key]]
    if (!is.null(chrom_store) && !identical(chrom_store, NA)) {
      chroms <- tryCatch(
        onDiskData_retrieve(chrom_store),
        error = function(e) {
          warning("Could not retrieve ", chrom_key, ": ", conditionMessage(e))
          NULL
        }
      )
    }

    out_file <- file.path(project_dir, paste0("Feature_group_EIC_", pol, ".pdf"))
    if (file.exists(out_file)) {
      unlink(out_file)
    }

    message_with_time(
      "Reporting ", length(fg_keep), " feature-group EIC plots for ", pol,
      " -> ", basename(out_file)
    )

    n_ok <- 0L
    for (i in seq_along(fg_keep)) {
      g <- fg_keep[[i]]
      n_feat <- as.integer(fg_tab[[g]])
      page_w <- if (is.null(width)) {
        max(6, min(0.55 * n_feat + 2, 36))
      } else {
        width
      }
      page_h <- if (is.null(height)) {
        max(6, min(0.55 * n_feat + 2, 36))
      } else {
        height
      }

      p <- tryCatch(
        plot_xcms_feature_group_EIC_comparasion(
          xcms = xcms.xcms,
          feature_group = g,
          chroms = chroms,
          expandRt = expandRt,
          min_width = min_width,
          max_features = max_features,
          sample_index = sample_index
        ),
        error = function(e) {
          warning("Skip ", pol, " / ", g, ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(p)) next

      export_graph2pdf(
        p,
        file_path = out_file,
        append = n_ok > 0L,
        width = page_w,
        height = page_h
      )
      n_ok <- n_ok + 1L
    }

    if (n_ok > 0L) {
      out_files <- c(out_files, out_file)
      message_with_time("Wrote ", n_ok, " pages: ", out_file)
    } else {
      message("No pages written for ", pol)
    }
  }

  invisible(out_files)
}
