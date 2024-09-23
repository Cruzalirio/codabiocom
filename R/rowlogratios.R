#' @title rowlogratios
#' @description
#' parallel logarithmic ratios
#'
#' @param data data with counts with OTUs in the columns and samples in rows
#' @param group a vector with the sample groups
#' @param col the number of column for calculate the logratios
#'
#'@return \code{matlrAUC} A vector of size \eqn{1\times n} with the OTU's attainments in the given column
#' @examples
#' data(HIV)
#' AUC <- rowlogratios(x_HIV, 2, y_HIV)
#' AUC[1:10]
#'
#' @export
#' @importFrom pROC auc
#' @importFrom pROC roc
#' @importFrom coda4microbiome impute_zeros



rowlogratios <- function(data, col, group){
  #data = data.frame(data)
  data = coda4microbiome::impute_zeros(data)
  matlrAUC = matrix(0, nrow=1, ncol=ncol(data))
  for (j in 1:ncol(data)) {
    if(j>col){
      matlrcoljTemp <- log(data[, col]) - log(data[, j])
      matlrAUC[1, j] <- pROC::auc(pROC::roc(group,matlrcoljTemp, quiet = TRUE))[[1]]
    }
  }
  return(matlrAUC)
}

