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
## --------------------------

#' rfint()
#'
#' Implements seven different random forest prediction interval methods.
#'
#' The seven methods implemented are cited in the References section.
#' Additional information can be found within those references.
#' Each of these methods are implemented by utilizing the ranger package.
#' For \code{method = "Zhang"}, prediction intervals are generated using out-of-bag residuals.
#' \code{method = "Romano"} utilizes a split-conformal approach.
#' \code{method = "Roy"} uses a bag-of-predictors approach.
#' \code{method = "Ghosal"} performs boosting to reduce bias in the random forest, and estimates variance.
#' The authors provide multiple variants to their methodology.
#' \code{method = "Tung"} debiases feature selection and prediction. Prediction intervals are generated using quantile regression forests.
#' \code{method = "HDI"} delivers prediction intervals through highest-density interval regression forests.
#' \code{method = "quantile"} utilizes quantile regression forests.
#'
#' @author Chancellor Johnstone
#' @author Haozhe Zhang
#' @param formula Object of class formula or character describing the model to fit. Interaction terms supported only for numerical variables.
#' @param train_data Training data of class data.frame.
#' @param test_data Test data of class data.frame. Utilizes ranger::predict() to produce prediction intervals for test data.
#' @param method Choose what method to generate RF prediction intervals. Options are \code{method = c("Zhang", "quantile", "Romano", "Ghosal", "Roy", "Tung", "HDI")}. Defaults to \code{method = "Zhang"}.
#' @param alpha Significance level for prediction intervals. Defaults to \code{alpha = 0.1}.
#' @param num_trees Number of trees used in the random forest.
#' @param min_node_size Minimum number of observations before split at a node.
#' @param m_try Number of variables to randomly select from at each split.
#' @param concise If  concise = TRUE, only predictions output. Defaults to \code{concise = FALSE}.
#' @param seed Seed for random number generation. Currently not utilized.
#' @param prop Proportion of training data to sample for each tree. Only for \code{method = "Ghosal"}.
#' @param num_threads The number of threads to use in parallel. Default is the current number of cores.
#' @param symmetry True if constructing symmetric out-of-bag prediction intervals, False otherwise. Used only \code{method = "Zhang"}. Defaults to \code{symmetry = TRUE}.
#' @param calibrate If \code{calibrate = TRUE}, intervals are calibrated to achieve nominal coverage. Currently uses quantiles to calibrate. Only for \code{method = "Roy"}.
#' @param Roy_method Interval method for \code{method = "Roy"}.
#' Options are \code{Roy_method = c("quantile", "HDI", "CHDI")}.
#' @param variant Choose which variant to use. Options are \code{method = c("1", "2")}. Only for \code{method = "Ghosal"}.
#' @param Ghosal_num_stages Number of total stages. Only for \code{method = "Ghosal"}.
#' @param featureBias Remove feature bias. Only for \code{method = "Tung"}.
#' @param predictionBias Remove prediction bias. Only for \code{method = "Tung"}.
#' @param Tung_R Number of repetitions used in bias removal. Only for \code{method = "Tung"}.
#' @param Tung_num_trees Number of trees used in bias removal. Only for \code{method = "Tung"}.
#' @param interval_type Type of prediction interval to generate.
#' Options are \code{method = c("two-sided", "lower", "upper")}. Default is  \code{method = "two-sided"}.
#' @export
#' @importFrom Rdpack reprompt
#' @seealso \link[ranger]{ranger}
#' @seealso \link[rfinterval]{rfinterval}
#'
#' @return
#'    \item{\code{int}}{Default output. Includes prediction intervals for all methods in \code{methods}.}
#'    \item{\code{preds}}{Predictions for test data for all methods in \code{methods}. Output when \code{concise = FALSE}.}
#' @examples
#' \donttest{
#' library(piRF)
#'
#' #functions to get average length and average coverage of output
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
#' method_vec <- c("quantile", "Zhang", "Tung", "Romano", "Roy", "HDI", "Ghosal")
#' #generate train and test data
#' ratio <- .975
#' nrow <- nrow(airfoil)
#' n <- floor(nrow*ratio)
#' samp <- sample(1:nrow, size = n)
#' train <- airfoil[samp,]
#' test <- airfoil[-samp,]
#'
#' #generate prediction intervals
#' res <- rfint(pressure ~ . , train_data = train, test_data = test,
#'              method = method_vec,
#'              concise= FALSE,
#'              num_threads = 1)
#'
#' #empirical coverage, and average prediction interval length for each method
#' coverage <- sapply(res$int, FUN = getCoverage, response = test$pressure)
#' coverage
#' length <- sapply(res$int, FUN = getPILength)
#' length
#'
#' #plotting intervals and predictions
#' par(mfrow = c(2,2))
#' for(i in 1:7){
#'    col <- ((test$pressure >= res$int[[i]][,1]) *
#'    (test$pressure <= res$int[[i]][,2])-1)*(-1)+1
#'    plot(x = res$preds[[i]], y = test$pressure, pch = 20,
#'       col = "black", ylab = "true", xlab = "predicted", main = method_vec[i])
#'    abline(a = 0, b = 1)
#'    segments(x0 = res$int[[i]][,1], x1 = res$int[[i]][,2],
#'       y1 = test$pressure, y0 = test$pressure, lwd = 1, col = col)
#' }
#' }
#' @references
#' \insertRef{breiman2001random}{piRF}
#'
#' \insertRef{ghosal2018boosting}{piRF}
#'
#' \insertRef{meinshausen2006quantile}{piRF}
#'
#' \insertRef{romano2019conformalized}{piRF}
#'
#' \insertRef{roy2019prediction}{piRF}
#'
#' \insertRef{tung2014bias}{piRF}
#'
#' \insertRef{zhang2019random}{piRF}
#'
#' \insertRef{zhu2019hdi}{piRF}

rfint <- function(formula = formula,
                  train_data = NULL,
                  test_data = NULL,
                  method = "Zhang",
                  alpha = 0.1,
                  symmetry = TRUE,
                  seed = NULL,
                  m_try = 2,
                  num_trees = 500,
                  min_node_size = 5,
                  num_threads = parallel::detectCores(),
                  calibrate = FALSE,
                  Roy_method = "quantile",
                  featureBias = TRUE,
                  predictionBias = TRUE,
                  Tung_R = 5,
                  Tung_num_trees = 75,
                  variant = 1,
                  Ghosal_num_stages = 2,
                  prop = .618,
                  concise = TRUE,
                  interval_type = "two-sided"){

  #check for currently running parallel processes
  if(is.null(num_threads)){
    num_threads <- parallel::detectCores()
  } else {
    num_threads = num_threads
  }

  if(foreach::getDoParWorkers() != num_threads){
    clust <- parallel::makeCluster(num_threads)
    doParallel::registerDoParallel(clust)
  } else {

  }

  #list for all intervals; just gets all output for each interval
  res <- list()
  int <- list()
  pred <- list()

  #tracking list
  for(id in method){
    if(id == "Zhang"){
      #Zhang call; no one sided
      res[[id]] <- rfinterval::rfinterval(formula = formula, train_data = train_data, test_data = test_data,
                        method = "oob", alpha = alpha, symmetry = symmetry, seed = seed,
                        params_ranger = list(mtry = m_try, num.trees = num_trees, min.node.size = min_node_size,
                                             num.threads = num_threads))
      pred[[id]] <- res[[id]]$testPred
      int[[id]] <- res[[id]]$oob_interval
    } else if (id == "Roy"){
      #Roy, Larocque call
      res[[id]] <- RoyRF(formula = formula, train_data = train_data, pred_data = test_data, intervals = TRUE,
                   interval_method = Roy_method, alpha = alpha, num_trees = num_trees,
                   num_threads = num_threads, m_try = m_try, interval_type = interval_type)

      pred[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals
    } else if (id == "Ghosal"){
      #Ghosal, Hooker call
      res[[id]] <- GhosalBoostRF(formula = formula, train_data = train_data, pred_data = test_data, num_trees = num_trees,
                           min_node_size = min_node_size, m_try = m_try, keep_inbag = TRUE,
                           intervals = TRUE, alpha = alpha, prop = prop,
                           variant = variant, num_threads = num_threads, num_stages = Ghosal_num_stages,
                           interval_type = interval_type)

      pred[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals
    } else if (id == "Tung"){
      #Tung call
      res[[id]] <- TungUbRF(formula = formula, train_data = train_data, pred_data = test_data,
                      num_trees = num_trees, feature_num_trees = Tung_num_trees,
                      min_node_size = min_node_size, m_try = m_try, alpha = alpha, forest_type = "QRF",
                      featureBias = TRUE, predictionBias = TRUE,
                      R = Tung_R, num_threads = num_threads, interval_type = interval_type)

      pred[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals
    } else if (id == "Romano"){
      #Romano, Candes call
      res[[id]] <- CQRF(formula = formula, train_data = train_data, pred_data = test_data,
                  num_trees = num_trees, min_node_size = min_node_size,
                  m_try = m_try, keep_inbag = TRUE, intervals = TRUE,
                  num_threads = num_threads, alpha = alpha, interval_type = interval_type)

      pred[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals
    } else if (id == "HDI"){
      #HDI forest call
      res[[id]] <- HDI_quantregforest(formula = formula,train_data = train_data, test_data = test_data,
                                alpha = alpha, num_tree = num_trees, mtry = m_try,
                                min_node_size = min_node_size, max_depth = 10,
                                replace = TRUE, verbose = FALSE, num_threads = num_threads)

      pred[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals
    } else if (id == "quantile"){
      #quantile RF call
      #intermediate step to get quantReg model
      res[[id]] <- ranger::ranger(formula = formula, data = train_data,
                    num.trees = num_trees, mtry = m_try, min.node.size = min_node_size,
                    quantreg = TRUE, num.threads = num_threads)

      #need to do this with predict for quantReg with ranger
      pred[[id]] <- predict(res[[id]], test_data, type = "quantiles", quantiles = c(alpha/2, 0.5, 1 - alpha/2))$predictions
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

foreach::registerDoSEQ()
}



