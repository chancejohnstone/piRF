## ---------------------------
##
## Script name: Romano, Patterson, Candes 2018
##
## Purpose of script: implement quantile conformal regression described in RPC 2018
##
## Author: Chancellor Johnstone
##
## Date Created: 2019-09-27
##
## Copyright (c) Chancellor Johnstone, 2019
## Email: cjohnsto@iastate.edu
##
## ---------------------------
##
## Notes:
##
##
## --------------------------

#' implements RF prediction interval using split conformal prediction as outlined in Romano, Patterson, Candes 2018. Helper function.
#'
#' This function implements split conformal prediction intervals for RFs. Currently used in rfint().
#' @param formula Object of class formula or character describing the model to fit. Interaction terms supported only for numerical variables.
#' @param train_data Training data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Matches ranger() requirements.
#' @param pred_data Test data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Utilizes ranger::predict() to get prediction intervals for test data.
#' @param num_trees Number of trees.
#' @param min_node_size Minimum number of observations before split at a node.
#' @param m_try Number of variables to randomly select from at each split.
#' @param keep_inbag Saves matrix of observations and which tree(s) they occur in. Required to be true to generate variance estimates for Ghosal, Hooker 2018 method. *Should not be an option...
#' @param intervals Generate prediction intervals or not.
#' @param alpha Significance level for prediction intervals.
#' @param forest_type Determines what type of forest: regression forest vs. quantile regression forest. *Should not be an option...
#' @param num_threads The number of threads to use in parallel. Default is the current number of cores.
#' @param interval_type Type of prediction interval to generate.
#' Options are \code{method = c("two-sided", "lower", "upper")}. Default is  \code{method = "two-sided"}.
#' @keywords internal
#' @import stats
CQRF <- function(formula = NULL, train_data = NULL, pred_data = NULL, num_trees = NULL,
                 min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                 intervals = TRUE, alpha = NULL, forest_type = "RF", num_threads = NULL,
                 interval_type = NULL){

  #one sided intervals
  #if(interval_type == "two-sided"){
  #  alpha <- alpha
  #} else {
  #  alpha <- alpha*2
  #}

  #parse formula
  if (!is.null(formula)) {
    train_data <- parse.formula(formula, data = train_data, env = parent.frame())
  } else {
    stop("Error: Please give formula!")
  }

  #get dependent variable
  dep <- names(train_data)[1]

  if (forest_type == "QRF"){
    quantreg <- TRUE
  } else {quantreg <- FALSE}

  #partition train data into I1 and I2
  train <- sample(1:nrow(train_data), floor(nrow(train_data)/2))
  not_test <- 1:nrow(train_data) %in% train

  I1 <- train_data[train,]
  I2 <- train_data[!not_test,]

  #train on I1
  rf <- ranger::ranger(formula, data = I1, num.trees = num_trees,
               min.node.size = min_node_size, mtry = m_try,
               keep.inbag = keep_inbag, quantreg = TRUE,
               num.threads = num_threads)

  #get conformity scores for everything in I2
  quant_preds <- predict(rf, I2, type = 'quantiles', quantiles = c(alpha1, alpha2), num.threads = num_threads)

  E <- cbind(quant_preds$predictions[,1] - I2[,dep], I2[,dep] - quant_preds$predictions[,2])

  #collection of conformity scores
  maxE <- apply(E, FUN = max, MARGIN = 1)

  #get median as well...
  preds <- predict(rf, pred_data, type = "quantiles", quantiles = c(alpha1, 0.5, alpha2), num.threads = num_threads)
  Q <- quantile(maxE, probs = 1 - (alpha/2 + alpha/2))
  intervals <- cbind(preds$predictions[,1] - Q, preds$predictions[,3] + Q)

  return(list(preds = preds$predictions[,2], pred_intervals = intervals))
}

#' @keywords internal
fit_cqrf <- function(formula = NULL,
                     train_data = NULL,
                     alpha = .1,
                     interval_type = "two-sided",
                     split = .5,
                     num_trees = NULL,
                     min_node_size = NULL,
                     m_try = NULL,
                     keep_inbag = TRUE,
                     forest_type = "RF",
                     num_threads = NULL){

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

  #parse formula
  if (!is.null(formula)) {
    train_data <- parse.formula(formula, data = train_data, env = parent.frame())
  } else {
    stop("Error: Please give formula!")
  }

  #get dependent variable
  dep <- names(train_data)[1]

  if (forest_type == "QRF"){
    quantreg <- TRUE
  } else {quantreg <- FALSE}

  #partition train data into I1 and I2
  train <- sample(1:nrow(train_data), floor(nrow(train_data)*split))
  not_test <- 1:nrow(train_data) %in% train

  I1 <- train_data[train,]
  I2 <- train_data[!not_test,]

  #train on I1
  rf <- ranger::ranger(formula, data = I1, num.trees = num_trees,
                       min.node.size = min_node_size, mtry = m_try,
                       keep.inbag = keep_inbag, quantreg = TRUE,
                       num.threads = num_threads)

  #get conformity scores for everything in I2
  quant_preds <- predict(rf, I2, type = 'quantiles', quantiles = c(alpha1, alpha2), num.threads = num_threads)

  E <- cbind(quant_preds$predictions[,1] - I2[,dep], I2[,dep] - quant_preds$predictions[,2])

  #collection of conformity scores
  maxE <- apply(E, FUN = max, MARGIN = 1)

  return(list(rf = rf, conf = maxE))
}

#' @keywords internal
predict_cqrf <- function(model,
                         pred_data = NULL,
                         alpha = NULL,
                         num_threads = NULL,
                         intervals = TRUE,
                         interval_type = NULL){

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

  #browser()
  #get median as well...
  preds <- predict(model$rf, pred_data, type = "quantiles", quantiles = c(alpha1, 0.5, alpha2), num.threads = num_threads)
  Q <- quantile(model$conf, probs = 1 - (alpha/2 + alpha/2))
  intervals <- cbind(preds$predictions[,1] - Q, preds$predictions[,3] + Q)

  return(list(preds = preds$predictions[,2], pred_intervals = intervals))

}
