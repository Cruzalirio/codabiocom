#' @title rowlogratios
#' @description
#' parallel logarithmic ratios
#'
#' @param data data with counts with OTUs in the columns and samples in rows
#' @param group a vector with the sample groups
#' @param col the number of column for calculate the logratios
#' @param X a \eqn{n\times p} matrix of p covariates observed in each sample
#'
#'@return \code{matlrAUC} A vector of size \eqn{1\times n} with the OTU's attainments in the given column
#' @examples
#' data(HIV)
#' x_HIVImp <- zCompositions::cmultRepl(x_HIV, method="GBM",
#' output="p-counts",suppress.print=TRUE,z.warning=0.99)
#' Xnp <- model.matrix(y_HIV~MSM_HIV)
#' AUC <- rowlogratios(data= x_HIVImp, col= 2, group=y_HIV, X =Xnp)
#' AUC[1:10]
#'
#' @export
#' @importFrom pROC auc roc
#' @importFrom HandTill2001 multcap
#' @importFrom nnet multinom
#' @importFrom stats predict
#' @importFrom zCompositions  cmultRepl


rowlogratios <- function(data, col, group, X=NULL){
  #data = data.frame(data)
  matlrAUC = matrix(0, nrow=1, ncol=ncol(data))
    for (j in 1:ncol(data)) {
      if(j>col){
        matlrcoljTemp <- log(data[, col]) - log(data[, j])
        matlrcoljTemp <- cbind(X, matlrcoljTemp)
        model1 <- nnet::multinom(group ~ matlrcoljTemp, trace = F)
        result1 <- stats::predict(model1, type='probs')
        if(ncol(as.matrix(result1))==1){
          matlrAUC[1, j] <-
            pROC::auc(pROC::roc(group,result1, quiet = TRUE))[[1]]
        }else{
        matlrAUC[1, j] <-
          HandTill2001::auc(HandTill2001::multcap(response = group,
                                                  predicted = as.matrix(result1)))
        }
      }
    }
  return(matlrAUC)
}

