#' @title calcAUClr
#' @description
#' parallel AUC for data
#'
#' @param data abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param group a vector with the sample groups
#' @param cores a number of cores fo paralelization, if \code{cores=NULL}, \code{parallel::detectCores()-1} will be used
#' @return \code{res} the upper triangular of AUC between OTUS
#' @param binary \code{logical} TRUE if group are binary, FALSE if group are multinomial
#'
#' @examples
#' data(HIV)
#' AUC <- calcAUClr(x_HIV, y_HIV, cores=2)
#' AUC[1:10,1:10]
#' @export
#' @importFrom coda4microbiome impute_zeros
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ foreach
#' @importFrom foreach %dopar%



calcAUClr <- function(data, group, cores=NULL, binary=TRUE){
  # cores <- parallel::makeCluster(parallel::detectCores()-1, type='PSOCK') # grabs max available
  if(is.null(cores)){
    cores <- parallel::detectCores()-1
  }
  #cl <- parallel::makeCluster(getOption('cl.cores', cores))
  cl <- parallel::makePSOCKcluster(cores)
  doParallel::registerDoParallel(cl)

  #foreach::registerDoSEQ()
  res <- foreach::foreach(i = seq_len(ncol(data)),
                 .combine = rbind,
                 .multicombine = TRUE,
                 .inorder = FALSE,
                 .packages = c('doParallel', "pROC"),
                  .export = c("rowlogratios", "LRRelev")) %dopar% {
                   codabiocom::rowlogratios(data, i, group)
                 }

  parallel::stopCluster(cl)
  i <- NULL
  rm(i)
  return(res)
}
