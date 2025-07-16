#' @title calcAUClr
#' @description
#' Compute pairwise AUCs for all OTUs in compositional microbiome data.
#'
#' @param data Abundance matrix or data frame (rows are samples, columns are taxa/variables).
#' @param group A vector indicating the group of each sample.
#' @param cores Number of cores for parallelization. If \code{NULL}, uses \code{parallel::detectCores()-1}.
#' @param X An optional \eqn{n \times p} matrix of covariates for each sample.
#' @param conf.level The confidence level for intervals, not in percent (e.g., 0.95 for 95\% CI).
#' @param method The method to use: \code{"hanley"}, \code{"delong"}, or \code{"bootstrap"}.
#' @param rho Mean correlation assumed between pairs of AUCs.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{AUC}}{The upper triangular matrix of pairwise AUCs between OTUs.}
#'   \item{\code{VAR}}{The upper triangular matrix of the variances of each AUC.}
#' }
#'
#' @examples
#' data(HIV)
#' x_HIVImp <- zCompositions::cmultRepl(x_HIV, method = "GBM",
#'   output = "p-counts", suppress.print = TRUE, z.warning = 0.99)
#' Xnp <- model.matrix(y_HIV ~ MSM_HIV)
#' AUC <- calcAUClr(data = x_HIVImp, group = y_HIV, cores = 2, X = Xnp,
#' method = "hanley")
#' AUC$AUC[1:10, 1:10]
#'
#' @export
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ foreach %dopar%


calcAUClr <- function(data, group, cores = NULL, X = NULL,  conf.level = 0.95,
                      method = c("hanley", "delong", "bootstrap"),
                      rho = 0.1) {
  if (is.null(cores)) {
    cores <- parallel::detectCores() - 1
  }
  cl <- parallel::makePSOCKcluster(cores)
  doParallel::registerDoParallel(cl)

  res <- foreach::foreach(
    i = seq_len(ncol(data)),
    .combine = combine_fun,
    .init = list(
      AUC = matrix(nrow = 0, ncol = ncol(data)),
      VAR = matrix(nrow = 0, ncol = ncol(data))
    ),
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = c('doParallel', 'pROC', 'HandTill2001', 'nnet'),
    .export = c("rowlogratios")
  ) %dopar% {
    codabiocom::rowlogratios(data = data, group= group, col = i,
                             X = X, conf.level = conf.level,
                             method = method ,
                             rho = rho)
  }

  parallel::stopCluster(cl)

  # Asegura que nombres de filas y columnas coincidan
  colnames(res$AUC) <- colnames(data)
  rownames(res$AUC) <- colnames(data)
  colnames(res$VAR) <- colnames(data)
  rownames(res$VAR) <- colnames(data)

  # Completar triangulo inferior (matrices simétricas)
  res$AUC[lower.tri(res$AUC)] <- t(res$AUC)[lower.tri(res$AUC)]
  res$VAR[lower.tri(res$VAR)] <- t(res$VAR)[lower.tri(res$VAR)]

  diag(res$AUC) <- 0
  diag(res$VAR) <- 0

  return(list=c(AUC=res$AUC, VAR=res$VAR))
}

