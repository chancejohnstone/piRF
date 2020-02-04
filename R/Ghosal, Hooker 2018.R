## ---------------------------
##
## Script name: Boosting Random Forests to Reduce Bia
##
## Purpose of script: Implement boosted random forest introduced in Ghosal, Hooker 2018
##
## Author: Chancellor Johnstone
##
## Date Created: 2019-09-18
##
## Copyright (c) Chancellor Johnstone, 2019
## Email: cjohnsto@iastate.edu
##
## ---------------------------
##
## Notes:
## implementing boosted random forest, not specifically the alternative boosted forest described in Ghosal, Hooker 2018
## working to implementing alternative forest method
## multiple boosts functional
## need to add variant 2 and fix variance issues
## changed genRF function genCombRF
## adding in roxygen2 documentation
## --------------------------

require(ranger)
require(dplyr)
require(reshape2)

#' implements RF prediction interval method in Ghosal, Hooker 2018
#'
#' This function implements variant one of the prediction interval methods in Ghosal, Hooker 2018.
#' @param formula Object of class formula or character describing the model to fit. Interaction terms supported only for numerical variables.
#' @param train_data Training data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Matches ranger() requirements.
#' @param pred_data Test data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Utilizes ranger::predict() to get prediction intervals for test data.
#' @param num_trees Number of trees.
#' @param min_node_size Minimum number of observations before split at a node.
#' @param m_try Number of variables to randomly select from at each split.
#' @param keep_inbag Saves matrix of observations and which tree(s) they occur in. Required to be true to generate variance estimates for Ghosal, Hooker 2018 method.
#' @param intervals Generate prediction intervals or not. Defaults to FALSE.
#' @param alpha Significance level for prediction intervals.
#' @param prop Proportion of training data to sample for each tree. Currently variant 2 not implemented.
#' @param variant Choose which variant to use. Currently variant 2 not implemented.
#' @param num-stages Number of boosting stages. Functional for >= 2; variance estimates need adjustment for variant 2.
#' @param num_threads The number of threads to use in parallel. Default is the current number of cores.
#' @keywords prediction interval, random forest, boosting
#' @export
#' @examples
#' GhosalBoostRF <- function(formula = NULL, train_data = NULL, pred_data = NULL, num_trees = 500,
#' min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
#' intervals = FALSE, alpha = NULL, forest_type = "RF",
#' replace = TRUE, prop = 1, variant = 1,
#' num_threads = num_threads)
#' @noRd
GhosalBoostRF <- function(formula = NULL, train_data = NULL, pred_data = NULL, num_trees = NULL,
                          min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                          intervals = FALSE, alpha = NULL, prop = NULL, variant = 1,
                          num_threads = NULL, num_stages = NULL){

  #parse formula
  if (!is.null(formula)) {
    train_data <- ranger:::parse.formula(formula, data = train_data, env = parent.frame())
    #orders train_data by response
    #train_data <- train_data[order(train_data[,1]), ]
  } else {
    stop("Error: Please give formula!")
  }

  #get dependent variable
  dep <- names(train_data)[1]

  #stage 1 RF; no different than any other RF
  rf <- genCombRF(formula = formula, train_data = train_data,
              pred_data = pred_data, intervals = FALSE,
              num_trees = num_trees, min_node_size = min_node_size,
              m_try = m_try, weights = NULL, importance = "none",
              prop = prop, num_threads = num_threads)

  #call boosting function; adding in variants?
  boostRF <- boostStage(rf, formula = formula, train_data = train_data, pred_data = pred_data, num_trees = num_trees,
                        min_node_size = min_node_size, m_try = m_try, keep_inbag = TRUE,
                        intervals = TRUE,
                        alpha = alpha, weights = NULL, num_stages = 2,
                        prop = prop, num_threads = num_threads)

  #get intervals; do we need to include other stuff? not right now maybe...
  intervals <- GHVar(boostRF, train_data, pred_data, variant, dep, alpha, num_threads = num_threads)
}

#' generate stage 1 RF for Ghosal, Hooker RF implementation
#'
#' This function is primarily meant to be used within the GhosalBoostRF function. All parameters are same as in GhosalBoostRF().
#' @keywords random forest
#' @export
#' @examples
#' genCombRF <- function(formula = NULL, train_data = NULL, pred_data = NULL, num_trees = num_trees,
#' min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
#' intervals = TRUE,
#' alpha = alpha, forest_type = "RF", importance = "none" , weights = NULL,
#' replace = replace, prop = prop, inbag = NULL, num_threads = num_threads)
#' @noRd
genCombRF <- function(formula = NULL, train_data = NULL, pred_data = NULL, num_trees = num_trees,
                  min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                  intervals = TRUE,
                  alpha = NULL, importance = "none" , weights = NULL, prop = NULL, inbag = NULL, num_threads = NULL){

  #replace always set to FALSE for Ghosal, Hooker method...
  replace <- FALSE

  #generate feature weights first
  rf <- ranger(formula, data = train_data, num.trees = num_trees,
               min.node.size = min_node_size, mtry = m_try,
               keep.inbag = keep_inbag, importance = importance, split.select.weights = weights,
               replace = replace, sample.fraction = ifelse(replace, 1, prop), inbag = inbag, num.threads = num_threads)

}

#' generate stage 2 (and more) RF for Ghosal, Hooker RF implementation
#'
#' This function is primarily meant to be used within the GhosalBoostRF() function. All parameters are same as in GhosalBoostRF().
#' @keywords cats
#' @export
#' @examples
#' boostStage <- function(rf, formula = NULL, train_data = NULL, pred_data = NULL, num_trees = num_trees,
#' min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
#' intervals = TRUE,
#' alpha = alpha, forest_type = forest_type, weights = NULL, num_stages = 2,
#' replace = replace, prop = prop, quantreg = TRUE, num_threads = num_threads)
#' @noRd
#boosting function
boostStage <- function(rf, formula = NULL, train_data = NULL, pred_data = NULL, num_trees = num_trees,
                       min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                       intervals = TRUE,
                       alpha = alpha, weights = NULL, num_stages = 2,
                       prop = prop, num_threads = num_threads){

  #get dependent variable
  dep <- names(train_data)[1]

  #get resids from stage 1
  stage1_preds <- predict(rf, data = train_data, num.threads = num_threads)
  resids <- train_data[dep] - stage1_preds$predictions
  colnames(resids) <- "bias"

  #training data with residuals instead of true dependent variable
  drop <- names(train_data) %in% dep

  #next stage; forest with residuals
  boost_list <- list()
  num_boost <- num_stages - 1
  rf_list <- list()
  rf_sum <- 0
  tree_rf_sum <- matrix(0, nrow = nrow(pred_data), ncol = num_trees)
  train_rf_sum <- 0
  rf_sum_intervals <- 0

  #get inbag list from initial rf
  inbag_list <- rf$inbag.counts
  for(i in 1:num_boost){
    aug_train_data <- cbind(resids, train_data[,!drop])

    rf2 <- genCombRF(formula = bias ~. , train_data = aug_train_data,
                 pred_data = pred_data, intervals = TRUE,
                 num_trees = num_trees, min_node_size = min_node_size,
                 m_try = m_try, weights = weights, importance = "none", alpha = alpha,
                 prop = prop, inbag = inbag_list, num_threads = num_threads)

    #add to boost rf list
    boost_list[[i]] <- rf2

    #needs to be adjusted to allow for variant 2; keeps all stage estimates and inbag...

    #train_data predictions for MSE estimate
    train_rf_sum <- train_rf_sum + predict(boost_list[[i]], train_data, num.threads = num_threads)$predictions
    #combine original RF with boosted
    train_preds <- predict(rf, train_data, num.threads = num_threads)$predictions + train_rf_sum

    #need to sum all of the predictions from every boosted forest...
    rf_sum <- rf_sum + predict(boost_list[[i]], pred_data, num.threads = num_threads)$predictions
    preds <- predict(rf, pred_data, num.threads = num_threads)$predictions + rf_sum

    #need to get every tree prediction
    tree_rf_sum <- tree_rf_sum + predict(boost_list[[i]], pred_data, predict.all = TRUE, num.threads = num_threads)$predictions
    tree_preds <- predict(rf, pred_data, predict.all = TRUE, num.threads = num_threads)$predictions + tree_rf_sum

    #get bias of new predictions
    resids <- preds - pred_data[dep]
    colnames(resids) <- "bias"

  }

  #return(list(stage1rf = rf, boostrf = boost_list, preds = preds, pred_intervals = pred_intervals))
  return(list(stage1rf = rf, boostrf = boost_list,
              preds = preds, tree_preds = tree_preds, train_preds = train_preds,
              inbag = rf$inbag.counts))
}

#' generate prediction intervals for Ghosal, Hooker 2018 implementation.
#'
#' This function is primarily meant to be used within the GhosalBoostRF() function. All parameters are same as in GhosalBoostRF().
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords U statistics, random forest, prediction intervals
#' @export
#' @examples
#' GHVar <- function(boostRF, train_data, pred_data, variant, dep, alpha, num_threads = num_threads)
#' @noRd
#get variance estimate
GHVar <- function(boostRF, train_data, pred_data, variant, dep, alpha, num_threads = num_threads){
  #add variance estimate procedure for variant 2; requires estimates, and inbag for each stage...

  n <- nrow(train_data)
  pred_n <- nrow(pred_data)

  sum_cov <- 0
  var_est <- rep(0, times = pred_n)
  tree_var_est <- rep(0, times = pred_n)
  cov_est <- rep(0, times = pred_n)

  #needs to get predictions from boostRF
  tree_preds <- boostRF$tree_preds
  num_trees <- ncol(tree_preds)

  #getting in bag binary for each test point

  #test this; dont know which one is correct...
  in_bag <- unlist(boostRF$inbag)
  #dim(in_bag) <- c(num_trees, dim(train_data)[1])
  #in_bag <- t(in_bag)
  dim(in_bag) <- c(dim(train_data)[1], num_trees)
  in_bag <- in_bag >= 1

  for(i in 1:pred_n){
    sum_cov <- 0
    for(j in 1:n){
      sum_cov <- sum_cov + cov(in_bag[j,], tree_preds[i,])^2
    }

    #get a variance estimate for each prediction point
    cov_est[i] <- sum_cov

  }

  tree_var_est <- apply(tree_preds, FUN = var, MARGIN = 1)
  var_est <- cov_est + tree_var_est/num_trees

  pred_intervals <- NULL

  #need estimate of MSE; need (oob?) predictions of train_data
  mse_est <- sum((boostRF$train_preds - train_data[,dep])^2)/n

  pred_intervals <- cbind(boostRF$preds + qnorm(alpha/2)*sqrt(var_est + mse_est),
                          boostRF$preds - qnorm(alpha/2)*sqrt(var_est + mse_est))

  #return(list(pred_intervals = pred_intervals, tree_var_est = tree_var_est, var_est = var_est, cov_est = cov_est,
  #            mse = mse_est, tree = tree_preds, inbag = in_bag, all_cov = NULL, trees = num_trees))

  return(list(preds = boostRF$preds, pred_intervals = pred_intervals, var_est = var_est,
              mse = mse_est))
}

#variance much lower in trees when compared to boostedForest; check this out...
