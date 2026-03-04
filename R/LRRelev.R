#' @title LRRelev
#' @description
#' Compute parallel log-ratio separability indices for microbiome data.
#'
#' @param data A data frame or matrix with counts: OTUs in columns, samples in rows.
#' @param group A vector indicating the group of each sample.
#' @param taxa A vector with the taxonomic classification for each OTU.
#' @param otus A vector with the name of each OTU.
#' @param sample A vector with the name or ID of each sample.
#' @param threshold Minimum total count for an OTU to be included. OTUs below this are excluded.
#' @param X An optional \eqn{n \times p} matrix of covariates observed in each sample.
#' @param conf.level The confidence level for intervals, not in percent (e.g., 0.95 for 95\% CI).
#' @param method The method to use: \code{"hanley"}, \code{"delong"}, or \code{"bootstrap"}.
#' @param rho Mean correlation assumed between pairs of AUCs.
#' @param cores Number of cores for parallelization. If \code{NULL}, uses \code{parallel::detectCores()-1}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{data1Imp}}{Data with zeros imputed using GBM replacement.}
#'   \item{\code{OTUS}}{A data frame with OTUs and their association index.}
#'   \item{\code{Misery}}{A vector of OTUs with total counts less than the threshold.}
#'   \item{\code{uniqueOTUS}}{A vector of OTUs present in only one sample.}
#'   \item{\code{AUCs}}{The matrix of pairwise AUCs for OTUs.}
#'   \item{\code{OTUSRelev}}{A vector with the most relevant OTUs based on the index.}
#' }
#'
#' @examples
#' # Example with HIV dataset:
#' data(HIV)
#' Xnp <- model.matrix(y_HIV ~ MSM_HIV)
#' output1 <- LRRelev(data = x_HIV, sample = rownames(x_HIV), group = y_HIV,
#'   taxa = colnames(x_HIV), otus = colnames(x_HIV), cores = 2,
#'    X = Xnp, method = "hanley")
#'
#' @export
#'
#' @importFrom zCompositions cmultRepl
#' @importFrom stats qnorm quantile var
#' @importFrom utils combn



LRRelev <- function (data, sample, group, taxa, otus, threshold=2,
                     cores=NULL, X=NULL, conf.level = 0.95,
                     method = c("hanley", "delong"),
                     rho = 0.1){
  Misery <- as.vector(which(colSums(data) <= threshold))
  if (length(Misery) > 0) {
    data1 <- data[, -Misery]
    taxa1 <- taxa[-Misery]
    otus1 <- otus[-Misery]
  } else {
    data1 <- data
    taxa1 <- taxa
    otus1 <- otus
  }

  # Delete OTUs that appear in only one sample
  hit <- function(x) { min(c(x, 1)) }
  data0 <- as.matrix(data1)
  data0[] <- vapply(data0, hit, numeric(1))
  uniqueOTUs <- which(colSums(data0) == 1)
  if (length(uniqueOTUs) > 0) {
    data2 <- data1[, -uniqueOTUs]
    taxa2 <- taxa1[-uniqueOTUs]
    otus2 <- otus1[-uniqueOTUs]
  } else {
    data2 <- data1
    taxa2 <- taxa1
    otus2 <- otus1
  }

  # Zero-imputation
  if(sum(data2==0)>0){
  data1ZI <- zCompositions::cmultRepl(
    data2, method = "GBM", output = "p-counts",
    suppress.print = TRUE, z.warning = 0.99
  )}else{
    data1ZI <- data2
  }
    # Calculate AUC and VAR matrices
    res <- codabiocom::calcAUClr(data= data1ZI, group = group, cores = cores,
                                 X=X,conf.level = conf.level,
                                 method = method)

    # Order OTUs by AUC importance
    o <- order(colSums(abs(res$AUC)), decreasing = TRUE)
    M <- res$AUC[o, o]
    V <- res$VAR[o, o]

    MTemp <- M
    VTemp <- V

    maxrow <- ncol(M)
    colnames(M) <- o
    rownames(M) <- o
    colnames(V) <- o
    rownames(V) <- o

    LRS <- list(max_log_ratio = colnames(M)[which(M == max(abs(M)),
                                                  arr.ind = TRUE)[(2:1)]],
                names_max_log_ratio = colnames(data1ZI)[as.numeric(colnames(M)[which(M == max(abs(M)),
                                                                                  arr.ind = TRUE)[(2:1)]])],
                order_importance = o,
                name_most_import_variables = colnames(data1ZI)[o[1:maxrow]],
                association_logratio_y = M)

    assoc <- rep(0, ncol(data1ZI))
    var_assoc <- rep(0, ncol(data1ZI))
    for (m in 1:ncol(data1ZI)) {
      assoc[m] <- sum(M[1:m, 1:m]) / (m^2 - m)
    }
    assoc[1] <- NA
    var_assoc <- var_separability_index(V, max_k = maxrow, rho = rho)

    CISup <- assoc + qnorm((1+conf.level)/2) * sqrt(var_assoc)
    CIInf <- assoc - qnorm((1+conf.level)/2) * sqrt(var_assoc)
    maxim <- which.max(assoc)
  return(list(
    dataImp = data1ZI,
    OTUS = data.frame(
      otus = otus2[LRS$order_importance],
      assoc = assoc,
      assoc_var = var_assoc,
      ciLower = CIInf,
      ciUpper = CISup),
    Misery = otus[Misery],
    uniqueOTUS = otus1[uniqueOTUs],
    AUCs = MTemp,
    VARs = VTemp,
    OTUSRelev = otus2[LRS$order_importance[1:maxim]]
  ))
}



