## ---------------------------
##
## Script name: Bias-Corrected QRF for features and prediction
##
## Purpose of script: Code to implement Tung, Huang, Nyugen, Khan 2014;
##                    bias-correction based on feature selection and prediction; two separate functions
##
## Author: Chancellor Johnstone
##
## Date Created: 2019-09-09
##
## Copyright (c) Chancellor Johnstone, 2019
## Email: cjohnsto@iastate.edu
##
## ---------------------------
##
## Notes:
##
## include options for different number of trees in the featureBias vs predBias step?
## make parallel; optimize when parallel is used vs not, based on R value...
## --------------------------

#' implements RF prediction interval method in Tung, Huang, Nyugen, Khan 2014.
#'
#' This function implements the feature bias and prediction bias methods outlined in Tung 2014.
#' @param formula Object of class formula or character describing the model to fit. Interaction terms supported only for numerical variables.
#' @param train_data Training data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Matches ranger() requirements.
#' @param pred_data Test data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Utilizes ranger::predict() to get prediction intervals for test data.
#' @param num_trees Number of trees.
#' @param min_node_size Minimum number of observations before split at a node.
#' @param m_try Number of variables to randomly select from at each split.
#' @param keep_inbag Saves matrix of observations and which tree(s) they occur in. Required to be true to generate variance estimates for Ghosal, Hooker 2018 method. *Should not be an option...
#' @param intervals Generate prediction intervals or not.
#' @param featureBias perform feature bias step.
#' @param predictionBias perform prediction bias.
#' @param feature_num_tree number of trees to be used in ech random forest generated for feature bias step.
#' @param R number of RFs generated in feature bias stage of Tung 2014 prediction interval. Defualt is 10.
#' @param alpha Significance level for prediction intervals.
#' @param num_threads The number of threads to use in parallel. Default is the current number of cores.
#' @keywords prediction interval, random forest, debiasing, internal
#' @examples
#' TungUbRF <- function(formula = NULL, train_data = NULL, pred_data = NULL, num_trees = NULL,
#' min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
#' intervals = TRUE, feature_num_trees = NULL,
#' alpha = NULL, forest_type = "QRF", featureBias = TRUE, predictionBias = TRUE, R = NULL,
#' num_threads = NULL)
#' @noRd
#bias reduction; choice for feature and/or prediction bias correction
TungUbRF <- function(formula = NULL, train_data = NULL, pred_data = NULL, num_trees = NULL,
                    min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                    intervals = TRUE, feature_num_trees = NULL,
                    alpha = NULL, forest_type = "QRF", featureBias = TRUE, predictionBias = TRUE, R = NULL,
                    num_threads = NULL, interval_type = NULL){

  #garbage collection
  gc()

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

  #get feature weights; calls genRF function
  featureWeights <- NULL
  if(featureBias == TRUE) {
    featureWeights <- genWeights(formula = formula, train_data = train_data, pred_data = pred_data,
                                 feature_num_trees = feature_num_trees, min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                                 forest_type = "QRF", importance = "permutation", R = R, num_threads = num_threads)
  }

  #stage 1 rf generation
  rf <- genRF(formula = formula, train_data = train_data,
              pred_data = pred_data, intervals = FALSE,
              num_trees = num_trees, min_node_size = min_node_size,
              m_try = m_try, forest_type = "QRF", weights = featureWeights, importance = "none", num_threads = num_threads)

  if(predictionBias == TRUE) {
    rf <- predictionUbRF(rf, formula = formula, train_data = train_data, pred_data = pred_data, num_trees = num_trees,
                         min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                         intervals = TRUE, alpha = alpha, forest_type = "QRF", weights = featureWeights,
                         num_threads = num_threads, interval_type = interval_type)
  }

  return(list(preds = rf$preds[,2], pred_intervals = rf$preds[,c(1,3)], weights = featureWeights))
}

#' generate quantile RF
#'
#' This function is primarily meant to be used within the TungUbRF() function. All parameters are same as in TungUbRf().
#' @keywords random forest, internal
#' @noRd
#changes made to genRF; add to previous versions to maintain one function?
genRF <- function(formula = NULL, train_data = NULL, pred_data = NULL, num_trees = num_trees,
                  min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                  intervals = TRUE,
                  alpha = alpha, forest_type = "RF", importance = "none" , weights = NULL, num_threads = num_threads){

  #why does genRF need pred_data? it shouldn't if it is just generating a forest?
  #parse formula
  if (!is.null(formula)) {
    train_data <- ranger:::parse.formula(formula, data = train_data, env = parent.frame())
    #orders train_data by response
    #train_data <- train_data[order(train_data[,1]), ]
  } else {
    stop("Error: Please give formula!")
  }

  if (forest_type == "QRF"){
    quantreg <- TRUE
  } else {quantreg <- FALSE}

  #generate feature weights first
  rf <- ranger(formula, data = train_data, num.trees = num_trees,
               min.node.size = min_node_size, mtry = m_try,
               keep.inbag = keep_inbag, quantreg = quantreg,
               importance = importance, split.select.weights = weights,
               num.threads = num_threads)

}

#' generate weights for RF through feature bias reduction method outlined in Tung 2014.
#'
#' This function is primarily meant to be used within the TungUbRF() function. All parameters are same as in TungUbRf().
#' @keywords random forest, internal
#' @noRd
#call genRF in this function after sampling training data
genWeights <- function(formula = NULL, train_data = NULL, pred_data = NULL, feature_num_trees = feature_num_trees,
                       min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                       intervals = TRUE, alpha = alpha, forest_type = "RF", importance = "permutation",
                       R = R, num_threads = num_threads){

  artif_features <- list()

  #adjust this; currently includes an artifical response...
  vi <- matrix(0, nrow = R, ncol = dim(train_data)[2]*2-1)
  n <- dim(train_data)[1]

  #sample from training data within a feature; R replications
  vi <- foreach(i = 1:R, .combine = rbind, .export = c("genRF","ranger")) %dopar% {
    #need to remove dependent variable?
    artif_features[[i]] <- apply(train_data, FUN = sample, MARGIN = 2, size = n, replace = FALSE)
    colnames(artif_features[[i]]) <- paste0("a_", names(train_data))

    #combine artifical and true train data
    artif_train_data <- cbind(train_data, artif_features[[i]])

    #retrieve variable importance from each replicate
    vi[i, ] <- genRF(formula = formula, train_data = artif_train_data,
                     pred_data = artif_train_data, intervals = FALSE,
                     num_trees = feature_num_trees, min_node_size = min_node_size,
                     m_try = m_try, forest_type = "QRF", importance = "permutation",
                     num_threads = num_threads)$variable.importance
    #return the variable importance
    vi
  }

  vi
  #get VI_hat
  vi_hat_all <- apply(vi, FUN = mean, MARGIN = 2)
  vi_hat <- vi_hat_all[1:(dim(train_data)[2]-1)]

  #normalized weights; [0,1]
  weights <- (vi_hat - min(vi_hat))/ (max(vi_hat) - min(vi_hat))

  #do this normalizing method for now
  #weights <- vi_hat / sum(vi_hat)

  return(weights)
}

#' performs prediction debiasing from Tung 2014
#'
#' This function is primarily meant to be used within the TungUbRF() function. All parameters are same as in TungUbRf().
#' @keywords random forest, internal
#' @noRd
#prediction bias correction; two stage random forest; takes first stage rf object as input
#rf <- test$rf$stage1rf
predictionUbRF <- function(rf, formula = NULL, train_data = NULL, pred_data = NULL, num_trees = NULL,
                           min_node_size = NULL, m_try = NULL, keep_inbag = TRUE,
                           intervals = TRUE,
                           alpha = alpha, forest_type = "QRF", weights = NULL,
                           num_threads = num_threads, interval_type = NULL){


  #one sided intervals
  if(interval_type == "two-sided"){
    alpha <- alpha
  } else {
    alpha <- alpha*2
  }

  #calibrate unused currently; implemented for QRF so foresttype not neccessary
  #generalize for RF?
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
  #print(dep)

  #keeping inbag by default; get oob for each tree
  #getting oob index for each training data point
  oob <- unlist(rf$inbag.counts)
  dim(oob) <- c(num_trees, dim(train_data)[1])
  oob <- t(oob)
  oob <- oob == 0
  oob[oob == FALSE] <- NA

  #get oob predictions; average over tree where each data point is oob
  all_tree_preds <- predict(rf, data = train_data, predict.all = TRUE, num.threads = num_threads)
  oob_preds <- all_tree_preds$predictions*oob

  #median of all oob predictions
  q_oob_preds <- apply(oob_preds, FUN = median, MARGIN = 1, na.rm = TRUE)

  #get bias based on median of collection of oob residuals
  bias <- as.matrix(q_oob_preds - train_data[,dep])
  colnames(bias) <- "bias"

  #get new predictions; always get median...
  stage1_preds <- predict(rf, pred_data, type = "quantiles", quantiles = c(alpha/2, .5, 1-alpha/2), num.threads = num_threads)

  drop <- names(train_data) %in% dep
  aug_train_data <- cbind(bias, train_data[,!drop])

  #stage 2
  rf2 <- genRF(formula = bias ~. , train_data = aug_train_data,
              pred_data = pred_data, intervals = FALSE,
              num_trees = num_trees, min_node_size = min_node_size,
              m_try = m_try, forest_type = forest_type, weights = weights, importance = "none", alpha = alpha,
              num_threads = num_threads)

  #just need the median for this; median predicted bias...
  pred_bias <- predict(rf2, pred_data, type = "quantiles", quantiles = .5, num.threads = num_threads)
  dims <- dim(stage1_preds$predictions)

  #repeating predicted bias to apply to quantile predictions
  mat_pred_bias <- matrix(pred_bias$predictions, nrow = dims[1], ncol = dims[2], byrow = FALSE)

  #bis corrected quantile predictions
  bias_correct_preds <- stage1_preds$predictions - mat_pred_bias

  #print(bias_correct_preds)

  return(list(stage1rf = rf, stage2rf = rf2, bias = bias, preds = bias_correct_preds))

  #should we work on outputting details related to the bias corrected predictions?
  #in a similar format to ranger?

}

