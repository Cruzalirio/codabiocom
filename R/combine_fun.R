utils::globalVariables(c("i"))
#' Combine lists of AUC and VAR matrices by row-binding
#'
#' This function combines multiple lists, each containing matrices named \code{AUC} and \code{VAR},
#' by row-binding the corresponding matrices together.
#'
#' @param ... Two or more list objects to combine. Each list must contain elements named \code{AUC} and \code{VAR},
#'   which are matrices of the same number of columns.
#'
#' @return \code{AUC} A matrix formed by row-binding the \code{AUC} matrices from the input lists.
#' @return \code{VAR} A matrix formed by row-binding the \code{VAR} matrices from the input lists.
#' @keywords internal

combine_fun <- function(...) {
  args <- list(...)
  n <- length(args)

  if (n == 0) stop("No arguments provided")
  if (n == 1) return(args[[1]])

  for (i in seq_len(n)) {
    if (!is.list(args[[i]]) ||
        !all(c("AUC", "VAR") %in% names(args[[i]]))) {
      stop(sprintf("Argument %d must be a list with elements 'AUC' and 'VAR'", i))
    }
  }

  res <- list(
    AUC = args[[1]]$AUC,
    VAR = args[[1]]$VAR
  )

  for (i in 2:n) {
    res$AUC <- rbind(res$AUC, args[[i]]$AUC)
    res$VAR <- rbind(res$VAR, args[[i]]$VAR)
  }

  return(res)
}



#' Approximate variance calculation for the separability index \(S_k\)
#'
#' Computes an approximate variance of the separability index \(S_k\) defined over
#' subsets of OTUs, accounting for a constant correlation \code{rho} between
#' pairwise AUC estimates.
#'
#' This function is intended for internal use only.
#'
#' @param V_matrix A symmetric matrix (m x m) of variances for pairwise OTU AUCs.
#'                 Only the lower triangle (pairs \(j < j'\)) is used.
#' @param max_k Integer. Maximum number of OTUs \eqn{k} for which to calculate the variance.
#'              Defaults to \code{ncol(V_matrix)}.
#' @param rho Numeric. Assumed constant correlation between pairs of AUCs. Default is 0.1.
#'
#' @return A numeric vector of length \code{max_k}, where the entry at position \eqn{k}
#'         contains the approximate variance of the separability index \(S_k\).
#'         The value for \eqn{k = 1} is \code{NA} as the index is undefined there.
#'
#' @keywords internal
var_separability_index <- function(V_matrix, max_k = ncol(V_matrix), rho = 0.1) {
  # Inicializar vector de varianzas
  var_vec <- numeric(max_k)
  var_vec[1] <- NA

  # Precalcular la matriz de desviaciones
  sqrt_V <- sqrt(V_matrix)
  V_upper <- V_matrix
  V_upper[lower.tri(V_upper, diag = TRUE)] <- 0
  sqrt_V_upper <- sqrt_V
  sqrt_V_upper[lower.tri(sqrt_V_upper, diag = TRUE)] <- 0

  # Iterar de manera vectorizada
  for (k in 2:max_k) {
    tri_idx <- upper.tri(V_matrix[1:k, 1:k])
    var_pairs <- V_matrix[1:k, 1:k][tri_idx]
    S1 <- sum(var_pairs)
    S <- sum(sqrt_V[1:k, 1:k][tri_idx])
    coef <- (2 / (k * (k - 1)))^2
    var_vec[k] <- coef * (S1 + rho * (S^2 - S1))
  }

  return(var_vec)
}

