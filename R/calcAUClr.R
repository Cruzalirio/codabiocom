#' @title calcAUClr
#' @description
#' parallel AUC for data
#'
#' @param data abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param group a vector with the sample groups
#' @param cores a number of cores fo paralelization, if \code{cores=NULL}, \code{parallel::detectCores()-1} will be used
#' @return \code{res} the upper triangular of AUC between OTUS
#' @param X a \eqn{n\times p} matrix of p covariates observed in each sample
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



calcAUClr <- function(data, group, cores=NULL, X =NULL){
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
                   codabiocom::rowlogratios(data, i, group, X)
                 }

  parallel::stopCluster(cl)
  i <- NULL
  rm(i)
  return(res)
}
