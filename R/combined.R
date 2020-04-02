## ---------------------------
##
## Script name: rfint()
##
## Purpose of script: incorporate all rf interval functions into single function
##
## Author: Chancellor Johnstone
##
## Date Created: 2020-01-13
##
## Copyright (c) Chancellor Johnstone, 2020
## Email: cjohnsto@iastate.edu
##
## ---------------------------
##
## Notes:
##
## allows for multiple methods to be chosen at once, but not multiple significance levels
## just does one set at a time, with predictions
##
## add in error messages; take HDI_quantregforest messages as start...
##
## Tung methodology very slow...
##
## better way to adjust options for different methods... options list?
## need to do this so we have more control over each function; calibrate() for example...
## add details section
## --------------------------

#source mean functions, var functions, data generating functions, all others...
#source("C:/Users/thechanceyman/Documents/RFIntervals/R/Simulation Functions.R")

#source the functions in each of the different sets of code
#source("C:/Users/thechanceyman/Documents/RFIntervals/R/Ghosal, Hooker 2018.R")
#source("C:/Users/thechanceyman/Documents/RFIntervals/R/HDI_quantregforest.R")
#source("C:/Users/thechanceyman/Documents/RFIntervals/R/Romano, Patterson, Candes 2018.R")
#source("C:/Users/thechanceyman/Documents/RFIntervals/R/Roy, Larocque 2019.R")
#source("C:/Users/thechanceyman/Documents/RFIntervals/R/Tung, Huang, Nyugen, Khan 2014.R")
#source("C:/Users/thechanceyman/Documents/RFIntervals/R/Ghosal, Hooker 2018.R")

###testing
#n <- 100
#ratio <- 2/3
#p <- 10
#full_n <- n*1/ratio
#x_data <- data_gen(full_n)
#colnames(x_data) <- paste0("X", 1:p)
#response <- linear_func(x_data) + rnorm(full_n, 0, 1)
#full_data <- cbind(response, x_data)

#train/test datasets
#subset <- sample(1:full_n, size = round(ratio*full_n))
#train <- as.data.frame(full_data[subset,])
#test <- as.data.frame(full_data[-subset,])
###

#' rfint()
#'
#' Single wrapper for seven different methods, and their variants, to create random forest prediction intervals. More methods to be added.
#' @param formula Object of class formula or character describing the model to fit. Interaction terms supported only for numerical variables.
#' @param train_data Training data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Matches ranger() requirements.
#' @param test_data Test data of class data.frame, matrix, dgCMatrix (Matrix) or gwaa.data (GenABEL). Utilizes ranger::predict() to get prediction intervals for test data.
#' @param method Choose what method to generate RF prediction intervals. ("Zhang", "quantile", "Romano", "Ghosal", "Roy", "Tung", "HDI). Defaults to "Zhang".
#' @param alpha Significance level for prediction intervals.
#' @param concise If TRUE, only predictions output. Defaults to FALSE
#' @param seed Seed for random number generation. Currently not utilized.
#' @param num_trees Number of trees used in the random forest.
#' @param min_node_size Minimum number of observations before split at a node.
#' @param m_try Number of variables to randomly select from at each split.
#' @param replace Sample with replacement, or not. Utilized for the two different variants outlined in Ghosal, Hooker 2018. Currently variant 2 not implemented.
#' @param prop Proportion of training data to sample for each tree. Currently variant 2 not implemented.
#' @param num_threads The number of threads to use in parallel. Default is the current number of cores.
#' @param symmetry True if constructing symmetric out-of-bag prediction intervals, False otherwise. Used only Zhang method. Defaults to TRUE.
#' @param calibrate If TRUE, intervals calibrated to achieve nominal coverage. Currently uses quantiles to calibrate. Only for method = "Roy".
#' @param variant Choose which variant to use. Only for method = "Ghosal". Currently variant 2 not implemented.
#' @param Ghosal_num_stages Number of total stages. Only for method = "Ghosal".
#' @param featureBias Remove feature bias. Only for method = "Tung".
#' @param predictionBias Remove prediction bias. Only for method = "Tung".
#' @param TungR Number of repetitions used in bias removal. Only for method = "Tung".
#' @param Tung_num_trees Number of trees used in bias removal. Only for method = "Tung".
#' @param interval_type Type of prediction interval to generate. Options = c("two-sided", "lower", "upper"). Default is "two-sided".
#' @return preds Predictions for test_data. Output with concise = FALSE. List format, organized by method name.
#' @return int Prediction intervals for test_data. List format, organized by method name.
#' @section Values
#' @keywords prediction interval, random forest, boosting, calibration
#' @importFrom parallel detectCores
#' @importFrom rfinterval rfinterval
#' @import ranger
#' @import foreach
#' @importFrom MASS ginv
#' @import doParallel
#' @export
#' @examples
#' library(piRF)
#'
#' getPILength <- function(x){
#' #average PI length across each set of predictions
#' l <- x[,2] - x[,1]
#' avg_l <- mean(l)
#' return(avg_l)
#' }
#'
#' getCoverage <- function(x, response){
#'   #output coverage for test data
#'   coverage <- sum((response >= x[,1]) * (response <= x[,2]))/length(response)
#'   return(coverage)
#' }
#'
#' #import airfoil self noise dataset
#' data(airfoil)
#'
#' #generate train and test data
#' ratio <- .9
#' nrow <- nrow(airfoil)
#' n <- floor(nrow*ratio)
#' samp <- sample(1:nrow, size = n)
#' train <- airfoil[samp,]
#' test <- airfoil[-samp,]
#'
#' #generate prediction intervals
#' res <- rfint(pressure ~ . , train_data = train, test_data = test,
#'              method = c("quantile", "Zhang", "Tung", "Romano", "Roy", "HDI", "Ghosal"))
#'
#' #empirical coverage, and average prediction interval length for each method
#' coverage <- sapply(res, FUN = getCoverage, response = test$pressure)
#' coverage
#' length <- sapply(res, FUN = getPILength)
#' length
rfint <- function(formula = formula, train_data = NULL, test_data = NULL, method = "Zhang", alpha = 0.1,
                  symmetry = TRUE, seed = NULL, m_try = 2, num_trees = 500,
                  min_node_size = 5, num_threads = detectCores(),
                  calibrate = FALSE, Roy_method = "quantile", featureBias = TRUE, predictionBias = TRUE,
                  Tung_R = 5, Tung_num_trees = 75, variant = 1, Ghosal_num_stages = 2, prop = .618,
                  concise = TRUE, interval_type = "two-sided", ...){

  #required libraries
  #require(parallel)
  #require(rfinterval)
  #require(ranger)
  #interval_type <- "upper"

  if(is.null(num_threads)){
    num_threads <- detectCores()
  } else {
    num_threads = num_threads
  }

  if(getDoParWorkers() != num_threads){
    clust <- makeCluster(num_threads)
    registerDoParallel(clust)
  }

  #list for all intervals; just gets all output for each interval
  res <- list()
  int <- list()
  pred <- list()

  #tracking list
  for(id in method){
    if(id == "Zhang"){
      #Zhang call; no one sided
      res[[id]] <- rfinterval(formula = formula, train_data = train, test_data = test,
                        method = "oob", alpha = alpha, symmetry = symmetry, seed = seed,
                        params_ranger = list(mtry = m_try, num.trees = num_trees, min.node.size = min_node_size,
                                             num.threads = num_threads))
      pred[[id]] <- res[[id]]$testPred
      int[[id]] <- res[[id]]$oob_interval
    } else if (id == "Roy"){
      #Roy, Larocque call
      res[[id]] <- RoyRF(formula = formula, train_data = train, pred_data = test, intervals = TRUE,
                   interval_method = Roy_method, alpha = alpha, num_trees = num_trees,
                   num_threads = num_threads, m_try = m_try, interval_type = interval_type)

      pred[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals
    } else if (id == "Ghosal"){
      #Ghosal, Hooker call
      res[[id]] <- GhosalBoostRF(formula = formula, train_data = train, pred_data = test, num_trees = num_trees,
                           min_node_size = min_node_size, m_try = m_try, keep_inbag = TRUE,
                           intervals = TRUE, alpha = alpha, prop = prop,
                           variant = variant, num_threads = num_threads, num_stages = Ghosal_num_stages,
                           interval_type = interval_type)

      pred[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals
    } else if (id == "Tung"){
      #Tung call
      res[[id]] <- TungUbRF(formula = formula, train_data = train, pred_data = test,
                      num_trees = num_trees, feature_num_trees = Tung_num_trees,
                      min_node_size = min_node_size, m_try = m_try, alpha = alpha, forest_type = "QRF",
                      featureBias = featureBias, predictionBias = predictionBias,
                      R = Tung_R, num_threads = num_threads, interval_type = interval_type)

      pred[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals
    } else if (id == "Romano"){
      #Romano, Candes call
      res[[id]] <- CQRF(formula = formula, train_data = train, pred_data = test,
                  num_trees = num_trees, min_node_size = min_node_size,
                  m_try = m_try, keep_inbag = TRUE, intervals = TRUE,
                  num_threads = num_threads, alpha = alpha, interval_type = interval_type)

      pred[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals
    } else if (id == "HDI"){
      #HDI forest call
      res[[id]] <- HDI_quantregforest(formula = formula,train_data = train, test_data = test,
                                alpha = alpha, num_tree = num_trees, mtry = m_try,
                                min_node_size = min_node_size, max_depth = 10,
                                replace = TRUE, verbose = FALSE, num_threads = num_threads)

      pred[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals
    } else if (id == "quantile"){
      #quantile RF call
      #intermediate step to get quantReg model
      res[[id]] <- ranger(formula = formula, data = train,
                    num.trees = num_trees, mtry = m_try, min.node.size = min_node_size,
                    quantreg = TRUE, num.threads = num_threads)

      #need to do this with predict for quantReg with ranger
      pred[[id]] <- predict(res[[id]], test, type = "quantiles", quantiles = c(alpha/2, 0.5, 1 - alpha/2))$predictions
      int[[id]] <- pred[[id]][,c(1,3)]
      #save only median
      pred[[id]] <- pred[[id]][,2]
    } else {
      stop(paste0(id, " is not a supported random forest prediction interval methodology"))
    }

    #one sided intervals; adjusts intervals based on on what sided interval is desired
    if(interval_type == "upper"){
      int[[id]][,1] <- -Inf
    } else if(interval_type == "lower"){
      int[[id]][,2] <- Inf
    }

    colnames(int[[id]]) <- c("lower", "upper")
  }

  #currently save list name as method name e.g. int$Zhang
  if(concise) {
    #output preds
    return(list(int = int))
  } else {
    #output oreds and intervals
    return(list(preds = pred, int = int))
  }

#stop cluster
#stopCluster(clust)
gc()
#getDoParWorkers()
}

#testing
#hold <- rfint(response ~ ., train_data = train, test_data = test, method = c("Roy"))



