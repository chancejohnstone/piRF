## ---------------------------
##
## Script name: rfweight.R
##
## Purpose of script:
##
## Author: Chancellor Johnstone
##
## Date Created: 2020-12-03
##
## Copyright (c) Chancellor Johnstone, 2020
## Email: cjohnsto@iastate.edu
##
## ---------------------------
##
## Notes:
## get rf weights for test data
##
## --------------------------
#' @keywords internal
rfweight <- function(res, train_data, test_data, bagged = TRUE){

  inbag <- res$inbag.counts
  num_tree <- length(inbag)

  #does not adjust for multiple copies of observation in bootstrapped sample
  if(!bagged){
    for(i in 1:num_tree){
      inbag[[i]] <- inbag[[i]] > 0
    }
  }

  n_test <- nrow(test_data)
  n_train <- length(inbag[[1]])
  w <- matrix(0, nrow = n_test, ncol = n_train)
  train_node <- predict(res, train_data, type = "terminalNodes")$predictions
  test_node <- predict(res, test_data, type = "terminalNodes")$predictions

  for (ind in 1:n_test) {

    weight = rep(0, n_train)

    for (tree in 1:num_tree) {

      node_match <- which(train_node[, tree] == test_node[ind, tree])
      weight[node_match] <- weight[node_match] + inbag[[tree]][node_match] / max(1,sum(inbag[[tree]][node_match]))

    }

    w[ind,] <- weight / num_tree
  }

  return(w)
}
