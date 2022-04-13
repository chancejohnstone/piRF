## ---------------------------
##
## Script name:
##
## Purpose of script: code to implement RF methods from Roy, Larocque 2019
##
## Author: Chancellor Johnstone
##
## Date Created: 2019-09-03
##
## Copyright (c) Chancellor Johnstone, 2019
## Email: cjohnsto@iastate.edu
##
## ---------------------------
##
## Notes:
##
## using original random forest implementation through ranger to implement methods proposed in paper.
## alternative splitting methods and prediction interval methods outlined in paper, not currently implemented
##
## --------------------------

#' implements RF prediction interval method in Roy, Larocque 2019. Helper function.
#'
#' Currently implemented is the quantile method with BOP intervals. Used inside rfint().
#' @param formula Object of class formula or character describing the model to fit. Interaction terms supported only for numerical variables.
#' @param train_data Training data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Matches ranger() requirements.
#' @param pred_data Test data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Utilizes ranger::predict() to get prediction intervals for test data.
#' @param num_trees Number of trees.
#' @param min_node_size Minimum number of observations before split at a node.
#' @param m_try Number of variables to randomly select from at each split.
#' @param keep_inbag Saves matrix of observations and which tree(s) they occur in. Required to be true to generate variance estimates for Ghosal, Hooker 2018 method. *Should not be an option...
#' @param intervals Generate prediction intervals or not.
#' @param interval_method which prediction interval type to generate. Several outlined in paper; currently only one method implemented.
#' @param calibrate calibrate prediction intervals based on out-of-bag performance. Adjusts alpha to get nominal coverage.
#' @param alpha Significance level for prediction intervals.
#' @param num_threads The number of threads to use in parallel. Default is the current number of cores.
#' @param interval_type Type of prediction interval to generate.
#' Options are \code{method = c("two-sided", "lower", "upper")}. Default is  \code{method = "two-sided"}.
#' @keywords internal
RoyRF <- function(formula = NULL, train_data = NULL, pred_data = NULL, num_trees = NULL,
                  min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                  intervals = TRUE, interval_method = "quantile", calibrate = FALSE, alpha = NULL, num_threads = NULL,
                  tolerance = NULL, step_percent = NULL, under = NULL, method = NULL,
                  max_iter = NULL, interval_type = NULL){

  #parse formula
  if (!is.null(formula)) {
    train_data <- parse.formula(formula, data = train_data, env = parent.frame())
  } else {
    stop("Error: Please give formula!")
  }

  #get dependent variable
  dep <- names(train_data)[1]

  rf <- ranger::ranger(formula, data = train_data, num.trees = num_trees,
                   min.node.size = min_node_size, mtry = m_try,
                   keep.inbag = keep_inbag, num.threads = num_threads)

  #predictions for test data
  rf_preds <- predict(rf, pred_data)$predictions

  if (intervals) {

    #get BOP for each new prediction
    all_BOP <- genBOP(rf, train_data = train_data, pred_data = pred_data,
                      alpha = alpha, num_threads = num_threads, calibrate = calibrate)

    BOP <- all_BOP$BOP
    oobBOP <- all_BOP$oobBOP
    dep <- all_BOP$dep

    if (calibrate) {
      #adjust alpha based on calibration; use calibrate() function
      alpha <- calibrate(oobBOP, alpha = alpha, response_data = train_data[,dep], tolerance = tolerance,
                         step_percent = step_percent, under = under, method = method,
                         max_iter = max_iter)
    }

    if(interval_method == "quantile") {
      int <- genqInt(BOP, alpha, interval_type = interval_type)
    } else if (interval_method == "HDI") {
      int <- genHDInt(BOP, alpha)
    } else if (interval_method == "CHDI") {
      int <- genCHDInt(BOP, alpha)
    }

    rf$int <- int
    rf$BOP <- BOP
    rf$oobBOP <- oobBOP
  }

  return(list(preds = rf_preds, pred_intervals = rf$int, alpha = alpha))
}

#' Generate BOP sets from Roy, Larocque 2019
#'
#' This function is primarily meant to be used within the RoyRF() function. All parameters are same as in RoyRF().
#' @keywords internal
genBOP <- function(rf, inbag = rf$inbag.counts, alpha = alpha,
                   pred_data, train_data, num_threads = num_threads,
                   calibrate = calibrate){

  #define %dopar% locally
  `%dopar%` <- foreach::`%dopar%`

  #declare index variables for foreach loop
  i <- j <- NULL

  #pred_data <- train_data
  npred <- dim(pred_data)[1]
  ntotal <- npred + dim(train_data)[1]
  B <- rf$num.trees

  #get dependent variable name
  dep <- names(train_data)[1]

  #terminal nodes for all new preds and training data
  term_nodes <- predict(rf, rbind(pred_data, train_data), type = "terminalNodes", num.threads = num_threads)$predictions

  #getting oob index for each training data point
  oob <- unlist(rf$inbag.counts)
  dim(oob) <- c(B, dim(train_data)[1])
  oob <- t(oob)
  oob <- oob == 0
  oob[oob == FALSE] <- NA

  #getting BOP list for each x value
  BOP <- list()
  oobBOP <- list()

  #done in parallel
  BOP <- foreach::foreach(i = 1:npred) %dopar% {
    #terminal node within each tree; compare to trees of training data only
    pred_nodes <- term_nodes[i,]

    #chance to matrix form
    pred_nodes <- t(matrix(pred_nodes, ncol = ntotal - npred, nrow = B))

    pred_same <- term_nodes[(npred+1):ntotal,] == pred_nodes
    pred_same[pred_same == FALSE] <- NA

    #get matching dependent variable values; doing oob as well
    hold <- as.vector(train_data[dep][,1])
    #dep_rep <- matrix(hold, ncol=B, nrow=length(hold), byrow=FALSE)
    dep_rep <- t(matrix(hold, nrow=B, ncol=length(hold), byrow = TRUE))
    dep_vals <- dep_rep * pred_same

    #get unique node values from each bag;
    #combine into one collection of BOP values for each new data point
    unq_dep_vals <- unlist(apply(dep_vals, FUN = unique, MARGIN = 2))
    BOP <- unq_dep_vals[is.na(unq_dep_vals) == FALSE]
    BOP
  }

  ntrain <- dim(train_data)[1]

  #computationally expensive with large train dataset
  #allow for a subsample to be random selected instead of oobBOP being generated for each train observation
  if(calibrate){
    oobBOP <- foreach::foreach(j = 1:ntrain) %dopar% {
      #terminal node within each tree; compare to trees of training data only
      pred_nodes <- term_nodes[(npred+j),]

      #chance to matrix form
      pred_nodes <- t(matrix(pred_nodes, ncol = ntotal - npred, nrow = B))

      pred_same <- term_nodes[(npred+1):ntotal,] == pred_nodes
      pred_same[pred_same == FALSE] <- NA

      #get matching dependent variable values; doing oob as well
      hold <- as.vector(train_data[dep][,1])
      dep_rep <- t(matrix(hold, nrow=B, ncol=length(hold), byrow = TRUE))
      oob_rep <- matrix(oob[j,], ncol=B, nrow=length(hold), byrow=TRUE)
      oob_dep_vals <- dep_rep * pred_same * oob_rep

      #get unique node values from each bag;
      #combine into one collection of BOP values for each new data point
      oob_unq_dep_vals <- unlist(apply(oob_dep_vals, FUN = unique, MARGIN = 2))
      oobBOP <- oob_unq_dep_vals[is.na(oob_unq_dep_vals) == FALSE]
      oobBOP
    }
  } else {oobBOP <- NULL}

  return(list(BOP = BOP, oobBOP = oobBOP, dep = dep))
}

#' Generates BOP quantile prediction intervals from Roy, Larocque 2019.
#'
#' This function is primarily meant to be used within the RoyRF() function.
#' @keywords internal
genqInt <- function(BOP, alpha = alpha, interval_type = interval_type){

  #one sided intervals
  #if(interval_type == "two-sided"){
  #  alpha <- alpha
  #} else {
  #  alpha <- alpha*2
  #}

  #one sided intervals
  if(interval_type == "two-sided"){
    alpha1 <- alpha/2
    alpha2 <- 1-alpha/2
  } else if(interval_type == "upper"){
    alpha1 <- 0
    alpha2 <- 1-alpha
  } else {
    alpha1 <- alpha
    alpha2 <- 1
  }

  q <- unlist(lapply(BOP, FUN = quantile, probs = c(alpha1, alpha2)))
  dim(q) <- c(2, length(q)/2)
  q <- t(q)

  return(q)

}

#' generates BOP HDI prediction intervals from Roy, Larocque 2019
#'
#' This function is primarily meant to be used within the RoyRF() function. Could potentially result in non-contiguous intervals.
#' @keywords internal
genHDInt <- function(BOP, alpha = alpha){

  #getting hdr function from hdrcde package
  hdr <- hdrcde::hdr

  hdi <- list()
  hdi_list <- lapply(BOP, FUN = hdr, prob = 1*100-alpha*100)

  #extracting hdi from list
  for(i in 1:length(BOP)) {
    hdi[[i]] <- hdi_list[[i]]$hdr
  }

  return(hdi)
}

#' generates BOP contiguous HDI prediction intervals from Roy, Larocque 2019
#'
#' This function is primarily meant to be used within the RoyRF() function.
#' @keywords internal
genCHDInt <- function(BOP, alpha = alpha){
  #prediction intervals based on BOP HDI; connects HDI; uses density estimation...
  hdi <- genHDInt(BOP, alpha = alpha)
  right <- lapply(hdi, FUN = max)
  left <- lapply(hdi, FUN = min)

  int <- cbind(unlist(left), unlist(right))

  return(int)
}

#' @keywords internal
fit_bop <- function(formula = NULL,
                    train_data = NULL,
                    num_trees = NULL,
                    min_node_size = NULL,
                    m_try = NULL,
                    keep_inbag = TRUE,
                    intervals = TRUE,
                    num_threads = NULL){

  #parse formula
  if (!is.null(formula)) {
    train_data <- parse.formula(formula, data = train_data, env = parent.frame())
  } else {
    stop("Error: Please give formula!")
  }

  rf <- ranger::ranger(formula, data = train_data, num.trees = num_trees,
                       min.node.size = min_node_size, mtry = m_try,
                       keep.inbag = keep_inbag, num.threads = num_threads)

  rf$train_data <- train_data
  rf$num_threads <- num_threads

  return(rf)

}

#' @keywords internal
predict_bop <- function(model,
                        pred_data,
                        interval_method = "quantile",
                        calibrate = FALSE,
                        alpha = .1,
                        tolerance = NULL,
                        step_percent = NULL,
                        under = NULL,
                        max_iter = NULL,
                        interval_type = NULL,
                        num_threads = NULL){

  #predictions for test data
  rf_preds <- predict(model, pred_data)$predictions

  #get BOP for each new prediction
  all_BOP <- genBOP(model, train_data = model$train_data, pred_data = pred_data,
                    alpha = alpha, num_threads = model$num_threads, calibrate = calibrate)

  BOP <- all_BOP$BOP
  oobBOP <- all_BOP$oobBOP
  dep <- all_BOP$dep

  if (calibrate) {
    #adjust alpha based on calibration; use calibrate() function
    alpha <- calibrate(oobBOP, alpha = alpha, response_data = train_data[,dep], tolerance = tolerance,
                       step_percent = step_percent, under = under, method = method,
                       max_iter = max_iter)
  }

  if(interval_method == "quantile") {
    int <- genqInt(BOP, alpha, interval_type = interval_type)
  } else if (interval_method == "HDI") {
    int <- genHDInt(BOP, alpha)
  } else if (interval_method == "CHDI") {
    int <- genCHDInt(BOP, alpha)
  }

  model$int <- int
  model$BOP <- BOP
  model$oobBOP <- oobBOP

  return(list(preds = rf_preds, pred_intervals = model$int, alpha = alpha))

}


