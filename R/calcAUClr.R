#' @title calcAUClr
#' @description
#' Compute pairwise AUCs for all OTUs in compositional microbiome data.
#'
#' @param data Abundance matrix or data frame (rows are samples, columns are taxa/variables).
#' @param group A vector indicating the group of each sample.
#' @param cores Number of cores for parallelization. If \code{NULL}, uses \code{parallel::detectCores()-1}.
#' @param X An optional \eqn{n \times p} matrix of covariates for each sample.
#' @param conf.level The confidence level for intervals, not in percent (e.g., 0.95 for 95\% CI).
#' @param method The method to use: \code{"hanley"}, \code{"delong"}, or \code{"bootstrap"}.
#' @param rho Mean correlation assumed between pairs of AUCs.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{AUC}}{The upper triangular matrix of pairwise AUCs between OTUs.}
#'   \item{\code{VAR}}{The upper triangular matrix of the variances of each AUC.}
#' }
#'
#' @examples
#' data(HIV)
#' x_HIVImp <- zCompositions::cmultRepl(x_HIV, method = "GBM",
#'   output = "p-counts", suppress.print = TRUE, z.warning = 0.99)
#' Xnp <- model.matrix(y_HIV ~ MSM_HIV)
#' AUC <- calcAUClr(data = x_HIVImp, group = y_HIV, cores = 2, X = Xnp,
#' method = "hanley")
#' AUC$AUC[1:10, 1:10]
#'
#' @export
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ foreach %dopar%
#' @importFrom glmnet cv.glmnet
#' @importFrom stats predict



calcAUClr <- function(data, group, cores = NULL, X = NULL,
                      conf.level = 0.95,
                      method = c("hanley", "delong"),
                      rho = 0.1) {

  method <- match.arg(method)

  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }

  group <- as.factor(group)
  C <- length(levels(group))
  p <- ncol(data)

  logdata <- log(data)

  # --------------------------------------------------
  # 1️⃣ Ajuste penalizado UNA sola vez
  # --------------------------------------------------
  if (!is.null(X)) {

    Xmat <- as.matrix(X)

    if(C>2){
      cvfit <- glmnet::cv.glmnet(
        x = Xmat,
        y = group,
        family = "multinomial"
      )

      prob_X <- stats::predict(
        cvfit,
        newx = Xmat,
        s = "lambda.1se",
        type = "response"
      )[,,1]

    }else{
      cvfit <- glmnet::cv.glmnet(
        x = Xmat,
        y = group,
        family = "binomial"
      )

      prob_X <- predict(
        cvfit,
        newx = Xmat,
        s = "lambda.1se",
        type = "response"
      )[,1]

    }

  } else {
    prob_X <- NULL
  }

  # --------------------------------------------------
  # 2️⃣ Cluster multiplataforma
  # --------------------------------------------------
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # --------------------------------------------------
  # 3️⃣ Paralelización correcta
  # --------------------------------------------------
  results <- foreach::foreach(
    i = 1:(p-1),
    .packages = c("pROC"),
    .combine = 'c'
  ) %dopar% {

    AUC_row <- numeric(p)
    VAR_row <- numeric(p)

    for (j in (i+1):p) {

      score <- logdata[, i] - logdata[, j]

      if (C == 2) {

        roc_obj <- pROC::roc(group, score, quiet = TRUE)
        auc_val <- as.numeric(pROC::auc(roc_obj))

        if (method == "hanley") {

          n0 <- sum(group == levels(group)[1])
          n1 <- sum(group == levels(group)[2])

          Q1 <- auc_val / (2 - auc_val)
          Q2 <- 2 * auc_val^2 / (1 + auc_val)

          var_val <- (auc_val * (1 - auc_val) +
                        (n1 - 1) * (Q1 - auc_val^2) +
                        (n0 - 1) * (Q2 - auc_val^2)) /
            (n1 * n0)

        } else {

          ci_auc <- pROC::ci.auc(
            roc_obj,
            method = "delong",
            conf.level = conf.level
          )

          se <- (ci_auc[3] - ci_auc[1]) /
            (2 * qnorm((1 + conf.level)/2))

          var_val <- se^2
        }

      } else {

        class_pairs <- combn(levels(group), 2)
        auc_pairs <- numeric(ncol(class_pairs))
        var_pairs <- numeric(ncol(class_pairs))

        for (k in 1:ncol(class_pairs)) {

          idx <- group %in% class_pairs[,k]
          g_bin <- droplevels(group[idx])

          roc_bin <- pROC::roc(g_bin,
                               score[idx],
                               quiet = TRUE)

          auc_bin <- as.numeric(pROC::auc(roc_bin))
          auc_pairs[k] <- auc_bin

          if (method == "hanley") {

            n0 <- sum(g_bin == levels(g_bin)[1])
            n1 <- sum(g_bin == levels(g_bin)[2])

            Q1 <- auc_bin / (2 - auc_bin)
            Q2 <- 2 * auc_bin^2 / (1 + auc_bin)

            var_pairs[k] <- (auc_bin * (1 - auc_bin) +
                               (n1 - 1) * (Q1 - auc_bin^2) +
                               (n0 - 1) * (Q2 - auc_bin^2)) /
              (n1 * n0)

          } else {

            ci_auc <- pROC::ci.auc(
              roc_bin,
              method = "delong",
              conf.level = conf.level
            )

            se <- (ci_auc[3] - ci_auc[1]) /
              (2 * qnorm((1 + conf.level)/2))

            var_pairs[k] <- se^2
          }
        }

        auc_val <- 2 / (C * (C - 1)) * sum(auc_pairs)

        cov_sum <- 0
        for (a in 1:(length(var_pairs)-1)) {
          for (b in (a+1):length(var_pairs)) {
            cov_sum <- cov_sum +
              rho * sqrt(var_pairs[a] * var_pairs[b])
          }
        }

        var_val <- (4 / (C * (C - 1))^2) *
          (sum(var_pairs) + 2 * cov_sum)
      }

      AUC_row[j] <- auc_val
      VAR_row[j] <- var_val
    }

    list(list(i = i,
              AUC_row = AUC_row,
              VAR_row = VAR_row))
  }

  # --------------------------------------------------
  # 4️⃣ Reconstrucción segura
  # --------------------------------------------------
  AUCmat <- matrix(0, p, p)
  VARmat <- matrix(0, p, p)

  for (res in results) {
    i <- res$i
    AUCmat[i, ] <- res$AUC_row
    VARmat[i, ] <- res$VAR_row
  }

  AUCmat[lower.tri(AUCmat)] <- t(AUCmat)[lower.tri(AUCmat)]
  VARmat[lower.tri(VARmat)] <- t(VARmat)[lower.tri(VARmat)]

  colnames(AUCmat) <- colnames(data)
  rownames(AUCmat) <- colnames(data)
  colnames(VARmat) <- colnames(data)
  rownames(VARmat) <- colnames(data)

  return(list(AUC = AUCmat, VAR = VARmat))
}
