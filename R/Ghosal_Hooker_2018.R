## ---------------------------
##
## Script name: Boosting Random Forests to Reduce Bias
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
## multiple boosts functional; testing needs to occur; seems to be an issue with variance estimates
## variant 2 implemented
## --------------------------

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
#' @param num_stages Number of boosting stages. Functional for >= 2; variance estimates need adjustment for variant 2.
#' @param num_threads The number of threads to use in parallel. Default is the current number of cores.
#' @keywords internal
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
                          num_threads = NULL, num_stages = NULL,
                          interval_type = NULL){

  #parse formula
  if (!is.null(formula)) {
    train_data <- parse.formula(formula, data = train_data, env = parent.frame())
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

  #call boosting function
  boostRF <- boostStage(rf, formula = formula, train_data = train_data, pred_data = pred_data, num_trees = num_trees,
                        min_node_size = min_node_size, m_try = m_try, keep_inbag = TRUE,
                        intervals = TRUE, alpha = alpha, weights = NULL, num_stages = num_stages,
                        prop = prop, num_threads = num_threads, variant = variant)

  #get intervals
  intervals <- GHVar(boostRF, train_data, pred_data, variant, dep, alpha, num_threads = num_threads,
                     interval_type = interval_type)
}

#' generate stage 1 RF for Ghosal, Hooker RF implementation
#'
#' This function is primarily meant to be used within the GhosalBoostRF function. All parameters are same as in GhosalBoostRF().
#' @keywords internal
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
  rf <- ranger::ranger(formula, data = train_data, num.trees = num_trees,
               min.node.size = min_node_size, mtry = m_try,
               keep.inbag = keep_inbag, importance = importance, split.select.weights = weights,
               replace = replace, sample.fraction = ifelse(replace, 1, prop), inbag = inbag, num.threads = num_threads)

}

#' generate stage 2 (and more) RF for Ghosal, Hooker RF implementation
#'
#' This function is primarily meant to be used within the GhosalBoostRF() function. All parameters are same as in GhosalBoostRF().
#' @keywords internal
#' @noRd
#boosting function
boostStage <- function(rf, formula = NULL, train_data = NULL, pred_data = NULL, num_trees = num_trees,
                       min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                       intervals = TRUE, alpha = alpha, weights = NULL, num_stages = 2,
                       prop = prop, num_threads = num_threads, variant = NULL){

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
  #change inbag base don variant
  inbag_list <- rf$inbag.counts
  if(variant == 1){
    inbag_list <- rf$inbag.counts
  } else {
    inbag_list <- NULL
  }

  #change so that the original rf is i == 1...
  boost_list[[1]] <- rf
  for(i in 2:num_stages){
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

    #combine original RF predictions with boosted predictions
    train_preds <- predict(rf, train_data, num.threads = num_threads)$predictions + train_rf_sum

    #sum all of the predictions from every boosted forest...
    rf_sum <- rf_sum + predict(boost_list[[i]], pred_data, num.threads = num_threads)$predictions
    preds <- predict(rf, pred_data, num.threads = num_threads)$predictions + rf_sum
    #preds <- rf_sum


    #get every tree prediction
    tree_rf_sum <- tree_rf_sum + predict(boost_list[[i]], pred_data, predict.all = TRUE, num.threads = num_threads)$predictions
    tree_preds <- predict(rf, pred_data, predict.all = TRUE, num.threads = num_threads)$predictions + tree_rf_sum
    #tree_preds <- tree_rf_sum

    #get bias of new predictions
    resids <- preds - pred_data[dep]
    colnames(resids) <- "bias"

  }

  return(list(stage1rf = rf, boostrf = boost_list,
              preds = preds, tree_preds = tree_preds, train_preds = train_preds,
              inbag = rf$inbag.counts))
}

#' generate prediction intervals for Ghosal, Hooker 2018 implementation.
#'
#' This function is primarily meant to be used within the GhosalBoostRF() function. All parameters are same as in GhosalBoostRF().
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords internal
#' @examples
#' GHVar <- function(boostRF, train_data, pred_data, variant, dep, alpha, num_threads = num_threads)
#' @noRd
#get variance estimate
GHVar <- function(boostRF, train_data, pred_data, variant, dep, alpha, num_threads = num_threads, interval_type = interval_type){
  #add variance estimate procedure for variant 2; requires estimates, and inbag for each stage...

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



  #includes original rf
  num_stages <- length(boostRF$boostrf)

  n <- nrow(train_data)
  pred_n <- nrow(pred_data)

  sum_cov <- 0
  var_est <- rep(0, times = pred_n)
  tree_var_est <- rep(0, times = pred_n)
  cov_est <- rep(0, times = pred_n)

  #needs to get predictions from boostRF
  #maybe call this something different
  tree_preds <- boostRF$tree_preds
  num_trees <- ncol(tree_preds)

  #test this; dont know which one is correct...
  in_bag <- unlist(boostRF$inbag)
  dim(in_bag) <- c(dim(train_data)[1], num_trees)
  in_bag <- in_bag >= 1

  if(variant == 1){
    for(i in 1:pred_n){
      sum_cov <- 0
      for(j in 1:n){
        sum_cov <- sum_cov + cov(in_bag[j,], tree_preds[i,])^2
      }

      #get a variance estimate for each prediction point
      cov_est[i] <- sum_cov
    }
  } else {
    tree_preds <- array(0, dim = c(pred_n, num_trees, num_stages))
    for(k in 1:num_stages){
      tree_preds[,,k] <- predict(boostRF$boostrf[[k]], pred_data, predict.all = TRUE, num.threads = num_threads)$predictions
    }

    #pretty slow...
    for(i in 1:pred_n){
      keep_cov <- 0
      for(j in 1:n){
        sum_cov <- 0
        for(k in 1:num_stages){
          #get inbag for each stage; maybe dont have to do this every time?
          in_bag <- unlist(boostRF$boostrf[[k]]$inbag)
          dim(in_bag) <- c(dim(train_data)[1], num_trees)
          in_bag <- in_bag >= 1

          #covariance for each stage
          sum_cov <- sum_cov + cov(in_bag[j,], tree_preds[i,,k])
        }
      keep_cov <- keep_cov + sum_cov^2
      }
      cov_est[i] <- keep_cov
    }
  }

  #tree preds different size matrix depending on variant
  if(variant == 1){
    #variance of tree predictions
    tree_var_est <- apply(tree_preds, FUN = var, MARGIN = 1)
  } else {
    #summing variances of each stages tree predictions...
    tree_var_est <- apply(apply(tree_preds, FUN = var, MARGIN = c(1,3)), FUN = sum, MARGIN = 1)
  }

  #overall variance estimate
  var_est <- cov_est + tree_var_est/num_trees

  pred_intervals <- NULL

  #need estimate of MSE; need (oob?) predictions of train_data
  mse_est <- sum((boostRF$train_preds - train_data[,dep])^2)/n

  pred_intervals <- cbind(boostRF$preds + qnorm(alpha1)*sqrt(var_est + mse_est),
                          boostRF$preds + qnorm(alpha2)*sqrt(var_est + mse_est))

  return(list(preds = boostRF$preds, pred_intervals = pred_intervals, var_est = var_est,
              mse = mse_est))
}

