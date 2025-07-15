#' @title rowlogratios
#' @description
#' parallel logarithmic ratios
#'
#' @param data data with counts with OTUs in the columns and samples in rows
#' @param group a vector with the sample groups
#' @param col the number of column for calculate the logratios
#' @param X a \eqn{n\times p} matrix of p covariates observed in each sample
#' @param conf.level the width of the confidence interval as [0,1], never in percent. Default: 0.95, resulting in a 95% CI.
#' @param method 	the method to use: \code{hanley} , \code{delong} or \code{bootstrap}.
#' @param rho Mean of correlation between pairs of AUC.
#'
#'@return \code{AUC} A vector of size \eqn{1\times m} with the AUC between the OTU's attainments in the given column
#'@return \code{VAR} A vector of size \eqn{1\times m} with the variances of AUC between the OTU's attainments in the given column
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


rowlogratios <- function(data, group, col, X = NULL, conf.level = 0.95,
                         method = c("hanley", "delong", "bootstrap"),
                         rho = 0.1) {
  # Inicializa matrices
  n_OTUs <- ncol(data)
  matlrAUC <- matrix(NA, nrow = 1, ncol = n_OTUs)
  matlrVAR <- matrix(NA, nrow = 1, ncol = n_OTUs)

  for (j in 1:n_OTUs) {
    if (j > col) {
      # Construye log-ratio
      matlrcoljTemp <- log(data[, col]) - log(data[, j])
      matlrcoljTemp <- cbind(X, matlrcoljTemp)

      # Ajusta modelo multinomial
      model1 <- nnet::multinom(group ~ matlrcoljTemp, trace = FALSE)
      result1 <- stats::predict(model1, type = "probs")

      n_classes <- length(unique(group))

      # AUC binario
      if (n_classes == 2) {
        roc_obj <- pROC::roc(group, result1, quiet = TRUE)
        auc_val <- pROC::auc(roc_obj)[[1]]
        if(method=="hanley"){
          n0 <- sum(group == unique(group)[1])
          n1 <- sum(group == unique(group)[2])

          Q1 <- auc_val / (2 - auc_val)
          Q2 <- 2 * auc_val^2 / (1 + auc_val)
          var_val <- (auc_val * (1 - auc_val) + (n1 - 1) * (Q1 - auc_val^2) +
                        (n0 - 1) * (Q2 - auc_val^2)) / (n1 * n0)

        }else if(method=="delong"){
          auc_temp <- pROC::ci.auc(roc_obj, method="delong", conf.level=conf.level)
          var_val <- (auc_temp[3]-auc_temp[1])/(2*qnorm((1+conf.level)/2))
        }else{
          # auc_temp <- pROC::ci.auc(roc_obj, method="bootstrap", conf.level=conf.level)
          # var_val <- (auc_temp[3]-auc_temp[1])/(2*qnorm((1+conf.level)/2))
          var_val <- NA
        }

      } else {
        # AUC multiclass Hand & Till
        auc_obj <- HandTill2001::multcap(response = group, predicted = as.matrix(result1))
        auc_val <- HandTill2001::auc(auc_obj)

        # Approx var: get binary pairs
        class_pairs <- combn(levels(as.factor(group)), 2)
        auc_pairs <- numeric(ncol(class_pairs))
        var_pairs <- numeric(ncol(class_pairs))

        for (p in 1:ncol(class_pairs)) {
          c1 <- class_pairs[1, p]
          c2 <- class_pairs[2, p]
          binary <- group[group %in% c(c1,c2)]
          binary <- factor(binary, levels = c(c1, c2))
          pred_binary <- result1[group %in% c(c1,c2), which(levels(as.factor(group))==c1)]
          roc_bin <- pROC::roc(binary, pred_binary, quiet = TRUE)
          auc_bin <- pROC::auc(roc_bin)[[1]]
          auc_pairs[p] <- auc_bin
          if(method=="hanley"){
            n0 <- sum(binary == c1)
            n1 <- sum(binary == c2)
            Q1 <- auc_bin / (2 - auc_bin)
            Q2 <- 2 * auc_bin^2 / (1 + auc_bin)
            var_bin <- (auc_bin * (1 - auc_bin) + (n1 - 1) * (Q1 - auc_bin^2) + (n0 - 1) * (Q2 - auc_bin^2)) / (n1 * n0)
            var_pairs[p] <- var_bin
          }else if(method=="delong"){
            auc_temp <- pROC::ci.auc(roc_bin, method="delong", conf.level=conf.level)
            var_pairs[p] <- (auc_temp[3]-auc_temp[1])/(2*qnorm((1+conf.level)/2))
          }else{
            #auc_temp <- pROC::ci.auc(roc_bin, method="bootstrap", conf.level=conf.level)
            #var_pairs[p] <- (auc_temp[3]-auc_temp[1])/(2*qnorm((1+conf.level)/2))
            var_pairs[p] <- NA
          }
        }

        # Combine: mean AUC and variance with cov approx
        C <- n_classes
        auc_val <- 2/(C*(C-1)) * sum(auc_pairs)
        cov_sum <- 0
        for (p in 1:(length(var_pairs)-1)) {
          for (q in (p+1):length(var_pairs)) {
            cov_sum <- cov_sum + rho * sqrt(var_pairs[p] * var_pairs[q])
          }
        }
        var_val <- (4/(C*(C-1))^2) * (sum(var_pairs) + 2*cov_sum)
      }

      matlrAUC[1, j] <- auc_val
      matlrVAR[1, j] <- var_val
    }
  }
  return(list(AUC = matlrAUC, VAR = matlrVAR))
}


