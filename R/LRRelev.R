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



LRRelev <- function(data, sample, group, taxa, otus,
                    threshold = 2,
                    cores = NULL, X = NULL,
                    conf.level = 0.95,
                    method = c("hanley", "delong"),
                    rho = 0.1) {

  method <- match.arg(method)
  group <- as.factor(group)

  # --------------------------------------------------
  # 1️⃣ Filtrar OTUs con bajo conteo
  # --------------------------------------------------
  lowOTUs <- which(colSums(data) <= threshold)
  if(length(lowOTUs) > 0){
    data1 <- data[, -lowOTUs, drop = FALSE]
    taxa1 <- taxa[-lowOTUs]
    otus1 <- otus[-lowOTUs]
  } else {
    data1 <- data
    taxa1 <- taxa
    otus1 <- otus
  }

  # --------------------------------------------------
  # 2️⃣ Eliminar OTUs presentes solo en una muestra
  # --------------------------------------------------
  data0 <- pmax(data1, 1)
  uniqueOTUs <- which(colSums(data0 == 1) == nrow(data0))
  if(length(uniqueOTUs) > 0){
    data2 <- data1[, -uniqueOTUs, drop = FALSE]
    taxa2 <- taxa1[-uniqueOTUs]
    otus2 <- otus1[-uniqueOTUs]
  } else {
    data2 <- data1
    taxa2 <- taxa1
    otus2 <- otus1
  }

  # --------------------------------------------------
  # 3️⃣ Imputar ceros
  # --------------------------------------------------
  if(sum(data2 == 0) > 0){
    data_imp <- zCompositions::cmultRepl(
      data2, method = "GBM",
      output = "p-counts",
      suppress.print = TRUE, z.warning = 0.99
    )
  } else {
    data_imp <- data2
  }

  # --------------------------------------------------
  # 4️⃣ Calcular matrices AUC y VAR usando calcAUClr
  # --------------------------------------------------
  res <- calcAUClr(
    data = data_imp,
    group = group,
    cores = cores,
    X = X,
    conf.level = conf.level,
    method = method,
    rho = rho
  )

  AUCmat <- res$AUC
  VARmat <- res$VAR

  # --------------------------------------------------
  return(list(
    dataImp = data_imp,
    Misery = otus[lowOTUs],
    uniqueOTUS = otus1[uniqueOTUs],
    AUCs = AUCmat,
    VARs = VARmat,
    OTUS2 = otus2))
}


