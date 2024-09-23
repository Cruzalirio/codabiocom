#' @title calcAUClr
#' @description
#' parallel AUC for data
#'
#' @param data abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param group a vector with the sample groups
#' @return \code{max log-ratio}
#' @return \code{names max log-ratio}
#' @return \code{order of importance}
#' @return \code{name of most important variables}
#' @return \code{"association log-ratio with y}
#'
#' @examples
#' data(HIV)
#' AUC <- calcAUClr(x_HIV, y_HIV)
#' AUC$`association log-ratio with y`
#' @export
#' @importFrom coda4microbiome impute_zeros
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ foreach
#' @importFrom foreach %dopar%



calcAUClr <- function(data, group){
  data <- coda4microbiome::impute_zeros(data)
  cores <- parallel::makeCluster(parallel::detectCores()-1, type='PSOCK') # grabs max available
  cl <- parallel::makeCluster(getOption('cl.cores', cores))
  doParallel::registerDoParallel(cl)
  foreach::registerDoSEQ()
  res <- foreach::foreach(i = seq_len(ncol(data)),
                 .combine = rbind,
                 .multicombine = TRUE,
                 .inorder = FALSE,
                 .packages = c('doParallel', "pROC"),
                  .export = c("rowlogratios")) %dopar% {
                   codabiocom::rowlogratios(data, i, group)
                 }

  parallel::stopCluster(cl)
  res[lower.tri(res) ] <- t(res)[lower.tri(res) ]
  o <- order(colSums(abs(res)), decreasing = TRUE)
  M <- res[o, o]
  maxrow <- ncol(M)
  colnames(M) <- o
  rownames(M) <- colnames(M)
  results <- list(`max log-ratio` = colnames(M)[which(M == max(abs(M)),
                                                      arr.ind = TRUE)[(2:1)]],
                  `names max log-ratio` = colnames(data)[as.numeric(colnames(M)[which(M == max(abs(M)),
                                                                                      arr.ind = TRUE)[(2:1)]])],
                  `order of importance` = o,
                  `name of most important variables` = colnames(data)[o[1:maxrow]],
                  `association log-ratio with y` = M)
  return(results)
}
