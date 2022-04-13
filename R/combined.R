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
#' For \code{method = "oob"}, prediction intervals are generated using out-of-bag residuals.
#' \code{method = "cqrf"} utilizes a split-conformal approach.
#' \code{method = "bop"} uses a bag-of-predictors approach.
#' \code{method = "brf"} performs boosting to reduce bias in the random forest, and estimates variance.
#' The authors provide multiple variants to their methodology.
#' \code{method = "bcqrf"} debiases feature selection and prediction. Prediction intervals are generated using quantile regression forests.
#' \code{method = "hdi"} delivers prediction intervals through highest-density interval regression forests.
#' \code{method = "quantile"} utilizes quantile regression forests.
#'
#' @author Chancellor Johnstone
#' @author Haozhe Zhang
#' @param formula Object of class formula or character describing the model to fit. Interaction terms supported only for numerical variables.
#' @param train_data Training data of class data.frame.
#' @param test_data Test data of class data.frame. Utilizes ranger::predict() to produce prediction intervals for test data.
#' @param method Choose what method to generate RF prediction intervals. Options are \code{method = c("oob", "quantile", "cqrf", "brf", "bop", "bcqrf", "hdi")}. Defaults to \code{method = "oob"}.
#' @param alpha Significance level for prediction intervals. Defaults to \code{alpha = 0.1}.
#' @param num_trees Number of trees used in the random forest.
#' @param min_node_size Minimum number of observations before split at a node.
#' @param m_try Number of variables to randomly select from at each split.
#' @param concise If  concise = TRUE, only predictions output. Defaults to \code{concise = FALSE}.
#' @param seed Seed for random number generation. Currently not utilized.
#' @param prop Proportion of training data to sample for each tree. Only for \code{method = "brf"}.
#' @param num_threads The number of threads to use in parallel. Default is the current number of cores.
#' @param symmetry True if constructing symmetric out-of-bag prediction intervals, False otherwise. Used only \code{method = "oob"}. Defaults to \code{symmetry = TRUE}.
#' @param calibrate If \code{calibrate = TRUE}, intervals are calibrated to achieve nominal coverage. Currently uses quantiles to calibrate. Only for \code{method = "bop"}.
#' @param Roy_method Interval method for \code{method = "bop"}.
#' Options are \code{Roy_method = c("quantile", "hdi", "CHDI")}.
#' @param variant Choose which variant to use. Options are \code{method = c("1", "2")}. Only for \code{method = "brf"}.
#' @param Ghosal_num_stages Number of total stages. Only for \code{method = "brf"}.
#' @param featureBias Remove feature bias. Only for \code{method = "bcqrf"}.
#' @param predictionBias Remove prediction bias. Only for \code{method = "bcqrf"}.
#' @param Tung_R Number of repetitions used in bias removal. Only for \code{method = "bcqrf"}.
#' @param Tung_num_trees Number of trees used in bias removal. Only for \code{method = "bcqrf"}.
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
#' method_vec <- c("quantile", "oob", "bcqrf", "cqrf", "bop", "hdi", "brf")
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
#' #get current mfrow setting
#' opar <- par(mfrow = c(2,2))
#'
#' #plotting intervals and predictions
#' for(i in 1:7){
#'    col <- ((test$pressure >= res$int[[i]][,1]) *
#'    (test$pressure <= res$int[[i]][,2])-1)*(-1)+1
#'    plot(x = res$preds[[i]], y = test$pressure, pch = 20,
#'       col = "black", ylab = "true", xlab = "predicted", main = method_vec[i])
#'    abline(a = 0, b = 1)
#'    segments(x0 = res$int[[i]][,1], x1 = res$int[[i]][,2],
#'       y1 = test$pressure, y0 = test$pressure, lwd = 1, col = col)
#' }
#' par(opar)
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

rfint <- function(formula = NULL,
                  train_data = NULL,
                  test_data = NULL,
                  method = "oob",
                  alpha = 0.1,
                  symmetry = TRUE,
                  seed = NULL,
                  m_try = NULL,
                  num_trees = 500,
                  min_node_size = NULL,
                  num_threads = parallel::detectCores(),
                  calibrate = FALSE,
                  Roy_method = "quantile",
                  featureBias = FALSE,
                  predictionBias = TRUE,
                  Tung_R = 5,
                  Tung_num_trees = 75,
                  variant = 1,
                  Ghosal_num_stages = 2,
                  prop = .618,
                  concise = TRUE,
                  interval_type = "two-sided"){

  .Deprecated("pirf", msg = "rfint() is deprecated. Utilize pirf() to build your forest model(s) and predict.pirf() to get test predictions and prediction intervals.")

}

