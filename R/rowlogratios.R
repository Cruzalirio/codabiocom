#' @title rowlogratios
#' @description
#' parallel logarithmic ratios
#'
#' @param data data with counts with OTUs in the columns and samples in rows
#' @param group a vector with the sample groups
#' @param col the number of column for calculate the logratios
#' @param binary \code{logical} TRUE if group are binary, FALSE if group are multinomial
#'
#'@return \code{matlrAUC} A vector of size \eqn{1\times n} with the OTU's attainments in the given column
#' @examples
#' data(HIV)
#' AUC <- rowlogratios(x_HIV, 2, y_HIV)
#' AUC[1:10]
#'
#' @export
#' @importFrom pROC auc roc
#' @importFrom coda4microbiome impute_zeros
#' @importFrom HandTill2001 multcap
#' @importFrom nnet multinom
#' @importFrom stats predict


rowlogratios <- function(data, col, group, binary=TRUE){
  #data = data.frame(data)
  data = coda4microbiome::impute_zeros(data)
  matlrAUC = matrix(0, nrow=1, ncol=ncol(data))
  if(binary){
    for (j in 1:ncol(data)) {
      if(j>col){
        matlrcoljTemp <- log(data[, col]) - log(data[, j])
        matlrAUC[1, j] <- pROC::auc(pROC::roc(group,matlrcoljTemp, quiet = TRUE))[[1]]
      }
    }
  }else{
    for (j in 1:ncol(data)) {
      if(j>col){
        matlrcoljTemp <- log(data[, col]) - log(data[, j])
        model1=nnet::multinom(group ~ matlrcoljTemp, trace = F)
        result1= stats::predict(model1, matlrcoljTemp, type='probs')
        matlrAUC[1, j] =
          HandTill2001::auc(HandTill2001::multcap(response = group,
                                                  predicted = as.matrix(result1)))
      }
    }
  }

  return(matlrAUC)
}

