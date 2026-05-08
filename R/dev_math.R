#' @title Wrapper for t.test
#' @description A convenience wrapper that returns the p-value from t.test, or 1 on error.
#' @param ... arguments passed to `stats::t.test`
#'
#' @return numeric p-value (or 1 if the test fails)
#' @export
t.test_dev <- function(...){

  try.catch <-try(p.value <- t.test(...)$p.value,silent = T)
  if (grepl("Error " , try.catch)) {
    return(1)
  }
  return(p.value)
}





groupMz <- function(x,ppm = 10,return.type = c("vector","data.frame")){

 # x.table <- data.frame( mz = x)%>%
 #   dplyr::mutate(raw.order = 1:length(x))%>%
 #   dplyr::arrange(mz)%>%
 #   dplyr::mutate(mz.diff = c(diff(mz),0),
 #                 mz.ppm = mz.diff/mz*1e6,
 #                 mz.group = "")
 # i <- 1
 # i.group <- 1
 # this.group.idx <- c()
 # while(i <= nrow(x.table)){
 #   x.table$mz.ppm[i]
 #   if (x.table$mz.ppm[i] <ppm.thresh) {
 #     this.group.idx <- c(this.group.idx,i)
 #     x.table$mz.group[this.group.idx] <- paste0("ion_group",sprintf("%06d",i.group))
 #   }else{
 #     this.group.idx <-  c(this.group.idx,i)
 #     x.table$mz.group[this.group.idx] <- paste0("ion_group",sprintf("%06d",i.group))
 #     i.group <- i.group+1
 #     i <- i+1
 #     this.group.idx <- c()
 #     next
 #   }
#
 #   i <- i+1
 #   next
#
 # }
#
  x.na <- is.na(x)
  x.raw <- x
  x <- x[!x.na]

  return.type <- match.arg(return.type)
  x.group <- MsCoreUtils::group(x,ppm = ppm)

  x.group.na <- seq_along(x.raw)
  x.group.na[x.na] <- NA
  x.group.na[!x.na] <- x.group
  if (return.type == "vector")
    return(x.group.na)
  x.table <-data.frame(mz = x.raw,
                       mz.group = x.group.na) %>%
    dplyr::group_by(mz.group)%>%
    dplyr::mutate(mz.center = median(mz),
                  mz.diff = abs(mz.center - mz),
                  mz.ppm = mz.diff/mz*1e6,
                  mz.width = max(mz)-min(mz),
                  mz.width.ppm = mz.width/mz*1e6)%>%
    dplyr::ungroup()

  return(x.table)


}


cluster_rt <- function(rt,rt.tol){
  if (length(rt)==1) {
    return(1)
  }
  dist(rt)%>%
    hclust()%>%
    cutree(h = rt.tol)

}

cluster_ion <- function(mz,rt , ppm.thresh =10, rt.tol = 15){



  ion.df <- data.frame(mz,rt)%>%
    dplyr::mutate(mz.group = groupMz(mz,ppm = ppm.thresh))%>%
    dplyr::group_by(mz.group)%>%
    dplyr::mutate(rt.cluster = cluster_rt(rt,rt.tol))%>%
    dplyr::group_by(mz.group,rt.cluster)%>%
    dplyr::mutate(ion.cluster = paste0("Ion_cluster",sprintf("%08d",cur_group_id())))%>%
    dplyr::ungroup()%>%
    dplyr::pull(ion.cluster)
  return(ion.df)




}



#' @title Matrix Subtraction of Two Vectors
#' @description Expands two vectors into matrices and computes element-wise subtraction,
#'   creating a full difference matrix where each element `(i,j)` equals `v1[i] - v2[j]`.
#' @param v1 Numeric vector. Optionally can have names which will be used as row names.
#' @param v2 Numeric vector. Optionally can have names which will be used as column names.
#' @return A matrix of dimensions length(v1) x length(v2) containing all pairwise differences.
#'   Row and column names are preserved if present in the input vectors.
#' @export
#'
#' @examples
#' a <- 3:8
#' b <- 1:2
#' matrixSub(a,b)
matrixSub <- function(v1,v2){

  m1 <- matrix(rep(v1,length(v2)),nrow = length(v1),byrow = F)
  m1
  m2 <- matrix(rep(v2,length(v1)),ncol = length(v2),byrow = T)
  m2
  x <- m1-m2
  rownames(x) <- names(v1)
  colnames(x) <- names(v2)
  return(x)
}


cor.mtest <- function(x){

  index.map <- expand.grid(1:ncol(x),1:ncol(x))

  cor.test.x <- function(y,x){

    p <- cor.test(x[,index.map[y,1]],x[,index.map[y,2]])$p.value
    r <-cor(x[,index.map[y,1]],x[,index.map[y,2]])
    return(list(p=p,r=r))

  }

  cor.df<- BiocParallel::bplapply(1:nrow(index.map),cor.test.x,x)%>%
    data.table::rbindlist()

  cor.r.m <- cor.p.m <-
    matrix(ncol = ncol(x),nrow = ncol(x))%>%
    `rownames<-`(colnames(x))%>%
    `colnames<-`(colnames(x))

  cor.r.m[1:25] <- cor.df$r
  cor.p.m[1:25] <- cor.df$p
  return(list(p = cor.p.m,
              r = cor.r.m))



}


#' @title Extract Middle Portion of Vector
#' @description Similar to head() and tail(), but returns elements from the middle portion of a vector
#'   centered around the median position. Useful for examining central tendency of ordered data.
#' @param x A numeric or character vector from which to extract elements.
#' @param n Integer specifying the number of elements to return. Default is 10.
#' @return A vector containing up to n elements from the middle portion of the input vector.
#' @export
#'
median_part <- function(x, n = 10){

  x.indice <- seq_along(x)
  x.median <-x[floor((median(x.indice)-n/2) : (median(x.indice)+n/2))]
  if (length(x.median)>n) {
    x.median <- x.median[1:n]
  }
  x.median
}


mean_f <- function(x,f,...){
  sapply(split(x,f),mean,...)
}

median_f <- function(x,f,...){
  sapply(split(x,f),median,...)
}

sum_matrix <- function(...){

  mat.list <- list(...)
  if(length(mat.list)<2){
    return(mat.list[[1]])
  }

  mat <- mat.list[[1]]
  for (i in 2:length(mat.list)) {
    mat <- mat + mat.list[[i]]
  }

  return(mat)



}

add_matrix <- function(mat1,mat2){

  rn <- unique(c(rownames(mat1),rownames(mat2)))%>%
    groupStringFactor()%>%levels()
  cn <- unique(c(colnames(mat1),colnames(mat2)))%>%
    groupStringFactor()%>%levels()
  m1 <- get_matrix_value_fill_with_NA(mat1,rn,cn)
  m2 <- get_matrix_value_fill_with_NA(mat2,rn,cn)
  combined_array <- array(c(m1, m2), dim = c(nrow(m1), ncol(m2), 2))
  sum_matrix <- apply(combined_array, c(1, 2),
                       function(x) sum(x, na.rm = TRUE))
  dimnames(sum_matrix) <- dimnames(m1)
  return(sum_matrix)
}

#' @title Normalize Values Using Min-Max Scaling
#' @description Applies min-max normalization to scale values to the range `[0, 1]`.
#'   Works on vectors, matrices, and data frames (normalizes each row for 2D structures).
#' @param x A numeric vector, matrix, or data frame to be normalized.
#' @return The normalized data with the same structure as the input.
#'   For vectors, returns a numeric vector. For matrices/data frames, returns a matrix with values
#'   normalized row-wise.
#' @export
#'

normalize_max_min<- function(x){

  f <- function(z){
    (z-min(z,na.rm = T))/(max(z,na.rm = T)-min(z,na.rm = T))
  }
  if (is.vector(x)) {
    return(f(x))
  }
  if (is.matrix(x)|is.data.frame(x)) {
    apply(x, 1, f)%>%t%>%return()
  }

}



mz.range.ppm <- function(mz = 200 ,ppm = 5){

  cbind(mz-mz*ppm/1e6,  mz+mz*ppm/1e6 )


}


expand_range <- function(x= c(5,10),add = 0,multi = 0){

  x <- c(x[1]-diff(x)*multi-add,
         x[2]+diff(x)*multi+add)
  return(x)
}

sum_dup <- function(x){

  sum(duplicated(x))

}




gaussian_functioin <- function(x , a =1,b = 0,c = 0.5){

  a * exp(-(x-b)^2/2/c^2)
}




#' @title Match Ions by M/Z and Retention Time
#' @description Matches two lists of ions based on mass-to-charge ratio (m/z) and retention time (RT).
#'   Returns all matched pairs including multiple matches, with error values for each match.
#' @param mz1 Numeric vector of m/z values for the first ion list.
#' @param rt1 Numeric vector of retention times for the first ion list. Default is NA values of same length as mz1.
#' @param mz2 Numeric vector of m/z values for the second ion list.
#' @param rt2 Numeric vector of retention times for the second ion list. Default is NA values of same length as mz2.
#' @param mz.ppm Numeric tolerance for m/z matching in parts per million (ppm). Default is 10.
#' @param rt.tol Numeric tolerance for retention time matching. Default is Inf (no RT filtering).
#' @return A data frame with columns: ion1 (index in first list), ion2 (index in second list),
#'   mz.error (relative m/z error), and rt.error (absolute RT difference).
#' @export
#'
match_mz_rt <- function(mz1,rt1 =rep( NA,length(mz1)),
                        mz2,rt2 =rep( NA,length(mz2)),
                        mz.ppm = 10,
                        rt.tol = Inf){



  .f <- function(mz,rt,mz2,rt2,mz.ppm){

    mz.error <- abs(mz - mz2)/mz
    mz.matched <- which(mz.error < mz.ppm*1e-6)
    mz.error <- mz.error[mz.matched]
    rt.error <- abs(rt2[mz.matched]-rt)
    return(data.frame(ion2 =mz.matched,
                      mz.error = mz.error,
                      rt.error = rt.error))

  }

  match.list<- bpmapply(mz1,rt1,
                        FUN = .f,MoreArgs = list(mz2,rt2,mz.ppm),
                        BPPARAM = SerialParam(#workers = 31,
                                            progressbar = F),
                        SIMPLIFY=F)
  match.df <- data.table::rbindlist(match.list,idcol = "ion1")

  match.df <- match.df%>%
    dplyr::filter(rt.error < rt.tol |is.na(rt.error))
  return(match.df)

}


#' @title Match Ions by M/Z Only
#' @description Matches two lists of ions based solely on mass-to-charge ratio (m/z).
#'   Returns the closest match for each ion in the first list, handling multiple potential matches.
#' @param mz1 Numeric vector of m/z values for the first ion list to match.
#' @param mz2 Numeric vector of m/z values for the second ion list to match against.
#' @param mz.ppm Numeric tolerance for m/z matching in parts per million (ppm). Default is 10.
#' @return A numeric vector of the same length as mz1, containing indices of matched ions in mz2.
#'   NA values indicate no match was found within the specified tolerance.
#' @export
#'

match_mz <- function(mz1,mz2,mz.ppm = 10){

  if (length(mz1)==0) {
    return(NULL)
  }
  match.df <- match_mz_rt(mz1= mz1,
                          mz2= mz2,
                          mz.ppm = mz.ppm)
  match.df <- match.df%>%
    dplyr::group_by(ion1)%>%
    dplyr::slice_min(mz.error)

  mz.match <- rep(NA,length(mz1))
  mz.match[match.df$ion1] <- match.df$ion2

  return(mz.match)

}

match_mz_grid <- function(mz1, mz2 , ppm = 10 ){

  match.df <-  expand.grid(
    ion1 = seq_along(mz1),
    ion2 = seq_along(mz2)
  )
  match.df <- match.df[match.df$ion1!=match.df$ion2,]
  match.df$mz.diff <- mz2[match.df$ion2]-mz1[match.df$ion1]
  match.df$mz.mean <- (mz1[match.df$ion1]+mz2[match.df$ion2])/2
  match.df$mz.ppm <- abs( match.df$mz.diff /match.df$mz.mean)*1e6
  match.df <- match.df[match.df$mz.ppm < ppm,]

  return(match.df)

  ion1 <- rep(seq_along(mz1),each = length(mz2))
  ion2 <- rep(seq_along(mz2),times = length(mz1))
  mz.diff <- mz2[ion2]-mz1[ion1]
  mz.mean <- (mz1[ion1]+mz2[ion2])/2
  mz.ppm <- abs( mz.diff / mz.mean)*1e6
  idx <- mz.ppm < ppm
  data.frame(
    ion1 = ion1[idx],
    ion2 = ion2[idx]
  )


}


groupHclust <- function (x, maxDiff = 5){

  if (length(x)<2) {
    return(1)
  }

  d <- dist(x)
  hc <- hclust(d)
  fg <- cutree(hc,h = maxDiff)
  return(fg)

}





plot_density <- function(x){
  plot(density(x,na.rm =T))
}

which.mid <- function(x){

  x.mid <- median(x)
  which.min(abs((x-x.mid)))
}


get_formula_from_lm <- function(lm_fit) {
  # Extract coefficients
  coefficients <- coef(lm_fit)

  # Construct the formula string
  formula_string <- paste0(
    "y = ",
    round(coefficients[1], 3),  # Intercept
    paste(
      sapply(2:length(coefficients), function(i) {
        coef_value <- round(coefficients[i], 3)
        sign <- ifelse(coef_value < 0, " - ", " + ")  # Handle signs
        paste0(sign, abs(coef_value), " * ", names(coefficients)[i])
      }),
      collapse = ""  # Collapse terms into one string
    )
  )

  return(formula_string)
}


rdot_product <- function(x,y,...){

  ids <- intersect(which(!is.na(x[,1])),
                   which(!is.na(y[,1])))
  if(length(ids)){
    a <- x[ids,2]
    b <- y[ids,2]
    return(sum(a*b,na.rm = T)/
             sqrt(sum(a^2,na.rm = T))/sqrt(sum(b^2,na.rm = T)))
  }else{
    return(0)
  }
}


weighted_icc_lmer <- function(score_matrix, weights) {


  # Input checks
  if (!is.matrix(score_matrix)) stop("Input 'score_matrix' must be a matrix.")
  if (ncol(score_matrix) != length(weights)) stop("Length of 'weights' must match number of columns (raters).")


  # if all same
  if (all(apply(score_matrix, 1, function(row) length(unique(row))) == 1)) {
    if (length(unique(rowSums(score_matrix))) == 1) {
      #message("All raters agree and all subjects have the same score. ICC is undefined (0/0), returning NA.")
      return(NA)
    } else {
      #message("All raters agree. Perfect agreement. ICC = 1")
      return(1)
    }
  }

  # Add subject ID
  colnames(score_matrix) <- paste0("V",1:ncol(score_matrix))
  df <- as.data.frame(score_matrix)
  df$Subject <- factor(seq_len(nrow(df)))

  # Long format
  long_df <- tidyr::pivot_longer(df,
                          cols = -Subject,
                          names_to = "Rater",
                          values_to = "Score")

  # Map weights to each rater
  weights <- weights / sum(weights)
  weight_df <- data.frame(Rater = colnames(score_matrix), Weight = weights)
  long_df <- left_join(long_df, weight_df, by = "Rater")

  # Fit weighted linear mixed model
  model <- lme4::lmer(Score ~ 1 + (1 | Subject), data = long_df, weights = Weight)

  # Extract variance components and compute ICC
  vc <- as.data.frame(lme4::VarCorr(model))
  icc <- vc$vcov[vc$grp == "Subject"] / sum(vc$vcov)

  return(icc)
}



weighted_icc <- function(score_mat, weights) {

  if (ncol(score_mat) != length(weights)) {
    stop("Length of weights must match number of raters (columns in score_mat).")
  }


  if (all(apply(score_mat, 1, function(row) length(unique(row))) == 1)) {

    if (length(unique(rowSums(score_mat))) == 1) {
     # message("All raters agree and all subjects have the same score. ICC is undefined (0/0), returning NA.")
      return(NA)
    } else {
      #message("All raters agree. Perfect agreement. ICC = 1")
      return(1)
    }
  }


  weights <- weights / sum(weights)


  subj_mean <- as.numeric(score_mat %*% weights)
  grand_mean <- mean(subj_mean)


  MS_between <- sum((subj_mean - grand_mean)^2) / (nrow(score_mat) - 1)


  MS_within <- sum(apply(score_mat, 1, function(x) sum(weights * (x - sum(weights * x))^2))) / nrow(score_mat)


  icc <- MS_between / (MS_between + MS_within)

  return(icc)
}

rsd <- function(x, na.rm = TRUE) {
   sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
}

match_mz_foverlaps <- function(mz1, mz2, ppm.base = mz1, ppm = 10 ) {

  if (length(mz1 )>length(mz2)) {

    mz.max.ppm <- max(ppm.base)*ppm/1e6

    dt_x <- data.table(x = mz2,
                       ion2 = seq_along(mz2))
    dt_x[, `:=`(
      xmin = x - mz.max.ppm,
      xmax = x + mz.max.ppm
    )]

    data.table::setkey(dt_x, xmin, xmax)

    dt_y <- data.table(y = mz1,
                       ion1 = seq_along(mz1),
                       y_start = mz1, y_end = mz1)
    data.table::setkey(dt_y, y_start, y_end)

    res <- data.table::foverlaps(dt_y, dt_x, by.x = c("y_start", "y_end"),
                                 type = "within", nomatch = 0L)


    res[, ppmb := ppm.base[ion1] ]
    res[, mz.ppm := (y - x) / ppmb * 1e6]
    res <- res[abs(mz.ppm) < ppm]


    res <- res[,c("ion1","ion2","mz.ppm")]



    return(res[])
  }

  dt_x <- data.table(x = mz1,
                     ion1 = seq_along(mz1))
  dt_x[, `:=`(
    xmin = x - ppm.base * ppm / 1e6,
    xmax = x + ppm.base * ppm / 1e6
  )]
  data.table::setkey(dt_x, xmin, xmax)

  dt_y <- data.table(y = mz2,
                     ion2 = seq_along(mz2),
                     y_start = mz2, y_end = mz2)
  data.table::setkey(dt_y, y_start, y_end)

  res <- data.table::foverlaps(dt_y, dt_x, by.x = c("y_start", "y_end"),
                               type = "within", nomatch = 0L)


  res[, mz.ppm := abs((y - x) / ppm.base[ion1] * 1e6)]



  res <- res[,c("ion1","ion2","mz.ppm")]
  res[]
}


ttwdfs <- function(a,p){
  if (p > 1) {
    p <- p/100
  }
  a/(1-p)
}

irange.intersect <- function(IR1,IR2){

  hits <- IRanges::findOverlaps(IR1,IR2)
  ov <- IRanges::pintersect(
    IR1[queryHits(hits)],
    IR2[subjectHits(hits)]
  )
  IRanges::reduce(ov)
}
