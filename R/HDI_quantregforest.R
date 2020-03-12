# Title:  Highest Density Interval Regression Forest (HDI-forest)
# Author: Haozhe Zhang <haozhe.stat@gmail.com>
# Copyright 2019, Haozhe Zhang, All rights reserved.
# License: GPL-3
# Date First Created: 2019-09-13
# Reference:  Zhu, Lin, Jiaxin Lu, and Yihong Chen. "HDI-Forest: Highest Density Interval Regression Forest." arXiv preprint arXiv:1905.10101 (2019).

#' implements HDI RF prediction interval method in ...
#'
#' This function implements an HDI RF prediction interval method.
#' @param formula Object of class formula or character describing the model to fit. Interaction terms supported only for numerical variables.
#' @param train_data Training data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Matches ranger() requirements.
#' @param test_data Test data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Utilizes ranger::predict() to get prediction intervals for test data.
#' @param alpha Significance level for prediction intervals.
#' @param num_tree Number of trees.
#' @param min_node_size Minimum number of observations before split at a node.
#' @param mtry Number of variables to randomly select from at each split.
#' @param max_depth maximum depth of each tree in RF. ranger parameter.
#' @param replace Sample with replacement, or not. Utilized for the two different variants outlined in Ghosal, Hooker 2018. Currently variant 2 not implemented.
#' @param num_threads The number of threads to use in parallel. Default is the current number of cores.
#' @keywords prediction interval, random forest, boosting, BOP
#' @export
#' @examples
#' @noRd
HDI_quantregforest <- function(formula = NULL,
                               train_data = NULL,
                               test_data = NULL,
                               alpha = NULL,
                               num_tree = NULL,
                               mtry = NULL,
                               min_node_size = NULL,
                               max_depth = NULL,
                               replace = TRUE,
                               verbose = FALSE,
                               num_threads = NULL,
                               ...){

  ## sort the response of training data by ascending order
  if (!is.null(formula)) {
    train_data <- ranger:::parse.formula(formula, data = train_data, env = parent.frame())
    train_data <- train_data[order(train_data[,1]), ]
  } else {
    stop("Error: Please give formula!")
  }

  ## Some checks of inputs (inherited from quantregforest package)
  if(! class(train_data[,1]) %in% c("numeric","integer") )
    stop("The response must be numeric ")

  if(is.null(nrow(train_data[,-1])) || is.null(ncol(train_data[,-1])))
    stop("The training data contains no data ")

  if(is.null(nrow(test_data)) || is.null(ncol(test_data)))
    stop("The test data contains no data ")

  if (any(is.na(train_data))) stop("NA not permitted in the training data")
  if (any(is.na(test_data))) stop("NA not permitted in the test data")

  ## train a random forest via ranger
  rf <- ranger:::ranger(formula = formula,
                        data = train_data,
                        num.trees = num_tree,
                        mtry = mtry,
                        min.node.size = min_node_size,
                        max.depth = max_depth,
                        replace = replace,
                        sample.fraction = ifelse(replace, 1, 0.632),
                        keep.inbag = TRUE,
                        num.threads = num_threads)

  ## extract inbag and terminalnode matrices
  inbag <- rf$inbag.counts
  train_node <- predict(rf, train_data, type = "terminalNodes", num.threads = num_threads)$predictions
  test_node <- predict(rf, test_data, type = "terminalNodes", num.threads = num_threads)$predictions
  preds <- predict(rf, test_data)$predictions

  ## Compute random forest weights for test instances
  n_train <- nrow(train_data)
  n_test <- nrow(test_data)
  test_weights <- list()

  pred_int <- data.frame(matrix(nrow = n_test, ncol = 2))
  colnames(pred_int) <- c("lower", "upper")
  rownames(pred_int) <- rownames(test_data)

  for (ind in 1:n_test) {

    weight = rep(0, n_train)

    for (tree in 1:num_tree) {

      node_match <- which(train_node[, tree] == test_node[ind, tree])

      weight[node_match] <- weight[node_match] + inbag[[tree]][node_match] / sum(inbag[[tree]][node_match])
    }

    weight <- weight / num_tree

    test_weights[[ind]] <- weight

    lower_index <- 1
    upper_index <- 0
    prob <- 0
    lower_bound <- min(train_data[,1])
    upper_bound <- max(train_data[,1])
    interval_width <- upper_bound - lower_bound

    ## Find the narrowest interval under the coverage constraint
    for (lower_index in 1:n_train) {
      while (upper_index < n_train & prob < 1- alpha) {
        upper_index <- upper_index + 1
        prob <- prob + weight[upper_index]
      }
      if (prob >= 1 - alpha) {
        upper_index_optimal <- upper_index
        if (train_data[upper_index_optimal,1] - train_data[lower_index,1] < interval_width) {
          lower_bound <- train_data[lower_index,1]
          upper_bound <- train_data[upper_index_optimal,1]
          interval_width <- upper_bound - lower_bound
        }
      } else {
        break
      }
      prob <- prob - weight[lower_index]
    }

    pred_int$lower[ind] <- lower_bound
    pred_int$upper[ind] <- upper_bound

    if (verbose == TRUE) {
      print(pred_int[ind,])
    }

  }

  return(list(preds = preds, pred_intervals = pred_int, test_weights = test_weights))
}

## Test HDI_quantregforest function
#library(rfinterval)

#HDI_quantregforest(pm2.5~.,train_data = BeijingPM25[1:8000,],
#                   test_data = BeijingPM25[8001:8661,],
#                   alpha = 0.05,num_tree = 500, mtry = 8,
#                   min_node_size = 5, max_depth = 10,
#                   replace = TRUE,verbose = TRUE)
