#' @title calcAUClr
#' @description
#' parallel AUC for data
#'
#' @param data abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param group a vector with the sample groups
#' @param cores a number of cores fo paralelization, if \code{cores=NULL}, \code{parallel::detectCores()-1} will be used
#' @param X a \eqn{n\times p} matrix of p covariates observed in each sample
#' @param conf.level the width of the confidence interval as [0,1], never in percent. Default: 0.95, resulting in a 95% CI
#' @param method 	the method to use: 'hanley', 'delong' or 'bootstrap'
#' @param rho Mean of correlation between pairs of AUC
#' @return \code{AUC} the upper triangular of AUC between OTUS
#' @return \code{VAR} the upper triangular of the variance of each AUC between OTUS
#'
#' @examples
#' data(HIV)
#' x_HIVImp = zCompositions::cmultRepl(x_HIV, method="GBM",
#' output="p-counts",suppress.print=TRUE,z.warning=0.99)
#' Xnp <- model.matrix(y_HIV~MSM_HIV)
#' AUC <- calcAUClr(data = x_HIVImp, group = y_HIV, cores=2, X=Xnp)
#' AUC[1:10,1:10]
#' @export
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ foreach
#' @importFrom foreach %dopar%



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

  return(res)
}

