## ---------------------------
##
## Script name: RF Interval Generation
##
## potential package names: prRFer, prRFect, RFPRInter, PrInterRF, piRF
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
## need to fix formula parsing formula issues...
##
## outputs everything from each method, not just the intervals
## lots of stuff needs to be cleaned/organized/consolidated/other adjective
##
## add in error messages; take HDI_quantregforest messages as start...
##
## need to adjust Simulation_Execution.R to use this function
##
## Tung methodology very slow...
##
## fix alphas
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
#' @param num_trees Number of trees.
#' @param min_node_size Minimum number of observations before split at a node.
#' @param m_try Number of variables to randomly select from at each split.
#' @param keep_inbag Saves matrix of observations and which tree(s) they occur in. Required to be true to generate variance estimates for Ghosal, Hooker 2018 method. *Should not be an option...
#' @param intervals Generate prediction intervals or not.
#' @param alpha Significance level for prediction intervals.
#' @param forest_type Determines what type of forest: regression forest vs. quantile regression forest. *Should not be an option...
#' @param replace Sample with replacement, or not. Utilized for the two different variants outlined in Ghosal, Hooker 2018. Currently variant 2 not implemented.
#' @param prop Proportion of training data to sample for each tree. Currently variant 2 not implemented.
#' @param variant Choose which variant to use. Currently variant 2 not implemented.
#' @param num_threads The number of threads to use in parallel. Default is the current number of cores.
#' @param concise Predictions for each method output in addition to intervals. Defaults to FALSE.
#' @keywords prediction interval, random forest, boosting
#' @export
#' @examples
#' library(piRF)
#'
#' getPILength <- function(PI){
#' #average PI length across each set of predictions
#' l <- PI[,2] - PI[,1]
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
                  Zhang_method = "oob", symmetry = TRUE, seed = NULL,
                  m_try = 2, num_trees = 500, min_node_size = 5, num_threads = detectCores(),
                  calibrate = FALSE, featureBias = TRUE, predictionBias = TRUE,
                  Tung_R = 5, Tung_num_trees = 75, Ghosal_variant = 1, concise = FALSE, prop = .632){

  #required libraries
  require(parallel)
  require(rfinterval)
  require(ranger)

  #list for all intervals; just gets all output for each interval
  int <- list()

  #tracking list
  for(m in method){
    id <- m
    print(paste0("Generating prediction intervals using ", m, " methodology..."))

    if(m == "Zhang"){
      #Zhang call
      int[[id]] <- rfinterval(formula = formula, train_data = train, test_data = test,
                        method = Zhang_method, alpha = alpha, symmetry = symmetry, seed = seed,
                        params_ranger = list(mtry = m_try, num.trees = num_trees, min.node.size = min_node_size,
                                             num.threads = num_threads))$oob_interval
    } else if (m == "Roy"){
      #Roy, Larocque call
      int[[id]] <- RoyRF(formula = formula, train_data = train, pred_data = test, intervals = TRUE,
                   interval_type = "quantile", alpha = alpha, num_trees = num_trees,
                   num_threads = num_threads, m_try = m_try)$pred_intervals
    } else if (m == "Ghosal"){
      #Ghosal, Hooker call
      int[[id]] <- GhosalBoostRF(formula = formula, train_data = train, pred_data = test, num_trees = num_trees,
                           min_node_size = min_node_size, m_try = m_try, keep_inbag = TRUE,
                           intervals = TRUE, alpha = alpha, prop = prop,
                           variant = Ghosal_variant, num_threads = num_threads)$pred_intervals
    } else if (m == "Tung"){
      #Tung call
      int[[id]] <- TungUbRF(formula = formula, train_data = train, pred_data = test,
                      num_trees = num_trees, feature_num_trees = Tung_num_trees,
                      min_node_size = min_node_size, m_try = m_try, alpha = alpha, forest_type = "QRF",
                      featureBias = featureBias, predictionBias = predictionBias,
                      R = Tung_R, num_threads = num_threads)$pred_intervals
    } else if (m == "Romano"){
      #Romano, Candes call
      int[[id]] <- CQRF(formula = formula, train_data = train, pred_data = test,
                  num_trees = num_trees, min_node_size = min_node_size,
                  m_try = m_try, keep_inbag = TRUE, intervals = TRUE,
                  num_threads = num_threads, alpha = alpha)$pred_intervals
    } else if (m == "HDI"){
      #HDI forest call
      int[[id]] <- HDI_quantregforest(formula = formula,train_data = train, test_data = test,
                                alpha = alpha, num_tree = num_trees, mtry = m_try,
                                min_node_size = min_node_size, max_depth = 10,
                                replace = TRUE, verbose = FALSE, num_threads = num_threads)$pred_intervals
    } else if (m == "quantile"){
      #quantile RF call
      #intermediate step to get quantReg model
      res <- ranger(formula = formula, data = train,
                    num.trees = num_trees, mtry = m_try, min.node.size = min_node_size,
                    quantreg = TRUE, num.threads = num_threads)

      #need to do this with predict for quantReg with ranger
      int[[id]] <- predict(res, test, type = "quantiles", quantiles = c(alpha/2, 1 - alpha/2))$predictions
    } else {
      print(paste0("Error: ", m, " is not a supported random forest prediction interval methodology"))
    }

    colnames(int[[id]]) <- c("lower", "upper")

    if(!concise) {
      #output preds
      #not implemented
    }
  }

#currently save list name as method name e.g. int$Zhang
return(int)
}

#testing
#hold <- rfint(response~., train_data = train, test_data = test, method = c("Tung", "Romano", "HDI", "Ghosal", "quantile"))


