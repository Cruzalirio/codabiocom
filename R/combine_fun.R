
#' Approximate variance calculation for the separability index \(S_k\)
#'
#' Computes an approximate variance of the separability index \(S_k\) defined over
#' subsets of OTUs, accounting for a constant correlation \code{rho} between
#' pairwise AUC estimates.
#'
#' This function is intended for internal use only.
#'
#' @param y Value fo function
#' @param min_k Integer. Maximum number of OTUs \eqn{k} for which to calculate the variance.
#'              Defaults to \code{ncol(V_matrix)}.
#' @param enforce_shape description
#'
#' @return A numeric vector of length \code{max_k}, where the entry at position \eqn{k}
#'         contains the approximate variance of the separability index \(S_k\).
#'         The value for \eqn{k = 1} is \code{NA} as the index is undefined there.
#'
#' @keywords internal
#'
#' @importFrom stats lm resid coef
detect_knee_piecewise <- function(y, min_k = 5, enforce_shape = TRUE){

  x <- seq_along(y)
  n <- length(y)

  best_k <- NA
  best_rss <- Inf

  for(k in min_k:(n - min_k)){

    fit1 <- lm(y[1:k] ~ x[1:k])
    fit2 <- lm(y[(k+1):n] ~ x[(k+1):n])

    if(enforce_shape){
      slope1 <- coef(fit1)[2]
      slope2 <- coef(fit2)[2]

      # Require decreasing second segment
      if(slope2 >= 0) next
    }

    rss_total <- sum(resid(fit1)^2) + sum(resid(fit2)^2)

    if(rss_total < best_rss){
      best_rss <- rss_total
      best_k <- k
    }
  }

  if(is.na(best_k))
    best_k <- which.max(y)

  return(best_k)
}


