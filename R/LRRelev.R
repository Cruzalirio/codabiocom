#' @title LRRelev
#' @description
#' parallel logarithmic ratios
#'
#' @param data data with counts with OTUs in the columns and samples in rows
#' @param group a vector with the sample groups
#' @param taxa a vector with the taxon classification of each OTU
#' @param otus a vector with the name of each OTUS
#' @param sample name of each sample
#' @param threshold minimum count for the OTU entered for comparison
#' @param X a \eqn{n\times p} matrix of p covariates observed in each sample
#' @param cores a number of cores fo paralelization, if \code{cores=NULL}, \code{parallel::detectCores()-1} will be used
#' @return \code{data1Imp} data with imputed zeros
#' @return \code{OTUS} a dataframe with otus and their association index Cesc calculate
#' @return \code{Misery} a list of OTUS with counts less than the threshold
#' @return \code{uniqueOTUS} a list of OTUS with counts only in one sample
#' @return \code{AUCs} The matrix of AUCs by OTUS
#' @return \code{OTUSRelev} The list of relevant otus
#' @examples
#' data(HIV)
#' Xnp <- model.matrix(y_HIV~MSM_HIV)
#' output1 <- LRRelev(data=x_HIV, sample = rownames(x_HIV), group = y_HIV,
#'  taxa = colnames(x_HIV),otus = colnames(x_HIV), cores=2, X=Xnp)
#'
#'
#' @export
#'
#' @importFrom zCompositions  cmultRepl



LRRelev <- function (data, sample, group, taxa, otus, threshold=2, cores=NULL, X=NULL){
  #data = data.frame(data)
  # Delete OTUS with total count less of umbral
  Misery <- as.vector(which(colSums(data)<=threshold))
  if(length(Misery)>0){
    data1 <- data[,-Misery]
    taxa1 <- taxa[-Misery]
    otus1 <- otus[-Misery]
  }else{
    data1 <- data
    taxa1 <- taxa
    otus1 <- otus
  }
  # Delete OTUS whic appears only in an Sample
  hit <- function(x){min(c(x,1))}
  data0 <- as.matrix(data1)
  data0[] <- vapply(data0, hit, numeric(1))
  uniqueOTUS <- which(colSums(data0)==1)
  if(length(uniqueOTUS)>0){
    data2 <- data1[,-uniqueOTUS]
    taxa2 <- taxa1[-uniqueOTUS]
    otus2 <- otus1[-uniqueOTUS]
  }else{
    data2 <- data1
    taxa2 <- taxa1
    otus2 <- otus1
  }
  data1ZI = zCompositions::cmultRepl(data2, method="GBM",output="p-counts",
                                     suppress.print=TRUE,z.warning=0.99)
  res <- codabiocom::calcAUClr(data1ZI,group, cores, X)
  res[lower.tri(res) ] <- t(res)[lower.tri(res) ]
  colnames(res) <- names(data2)
  rownames(res) <- names(data2)
  diag(res) <- 0
  o <- order(colSums(abs(res)), decreasing = TRUE)
  M <- res[o, o]
  MTemp <- M
  maxrow <- ncol(M)
  colnames(M) <- o
  rownames(M) <- colnames(M)
  LRS <- list(`max log-ratio` = colnames(M)[which(M == max(abs(M)),
                                                      arr.ind = TRUE)[(2:1)]],
                  `names max log-ratio` = colnames(data)[as.numeric(colnames(M)[which(M == max(abs(M)),
                                                                                      arr.ind = TRUE)[(2:1)]])],
                  `order of importance` = o,
                  `name of most important variables` = colnames(data)[o[1:maxrow]],
                  `association log-ratio with y` = M)
  assoc <- rep(0,ncol(data1ZI))
  for (m in 1:ncol(data1ZI)){
    assoc[m] <- sum(LRS$`association log-ratio with y`[1:m,1:m])/(m^2-m)
  }
  maxim <- which.max(assoc)
  #LRImp <- sort(as.numeric(LRS$`order of importance`[1:maxim]))
  #data1Imp <- data1ZI[,LRImp]
  return(list(dataImp = data1ZI,
              OTUS = data.frame(otus=otus[LRS$`order of importance`],
                                assoc = assoc),
              Misery =otus[Misery], uniqueOTUS = otus1[uniqueOTUS],
              AUCs = MTemp,
              OTUSRelev = otus[LRS$`order of importance`[1:maxim]]))
}


