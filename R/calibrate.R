## ---------------------------
##
## Script name: calibrate()
##
## Purpose of script: general calibration procedure based off of obb collections
##
## Author: Chancellor Johnstone
##
## Date Created: 2020-01-21
##
## Copyright (c) Chancellor Johnstone, 2020
## Email: cjohnsto@iastate.edu
##
## ---------------------------
##
## Notes:
## Currently specific to Roy methodology; attempting to generalize...
## might be better to output all oob stuff for each function...
## does not currently emphasize overcoverage over undercoverage...
## other methods might not have potential for calibration, if they already are using oob to generate intervals...
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
#m_try <- 10
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
#response_data <- train$response
#r <- RoyRF(response~., train_data = train, pred_data = test, intervals = TRUE,
#           interval_type = "quantile", alpha = .1, num_trees = 500,
#           num_threads = 1, m_try = 5)
#rf <- r$rf
#oob  <- r$rf$oobBOP
#alpha <- .1
#dep <- "response"
#test <- as.data.frame(full_data[-subset,])
#tolerance <- .01
###

#' calibrate()
#'
#' This function outputs a calibrated significance level based on coverage of prediction intervals generated using oob collections. Primarily for use in RoyRF(). Attempting to see which other methods could utilize this procedure.
#' @param oob collection of oob predictions for training data (in list form).
#' @param alpha nominal significance level. Defaults to 0.01.
#' @param response_data response data of class data.frame. Must have names() attribute.
#' @param tolerance tolerance allowed around nominal alpha. Default is 0.25.
#' @param step_ratio ratio absolute difference between empirical oob coverage and nominal coverage to adjust when calibrating. Defaults to 0.618.
#' @param under Allow undercoverage. Defaults to TRUE. Not currently implemented.
#' @param method Method to calibrate prediction intervals with. Defaults to "quantile"). Current only "quantile" implemented.
#' @param max_iter Maximum number of iterations. Defaults to 10.
#' @keywords random forest, calibration
#' @export
#' @examples
#' calibrate <- function(oob, alpha = alpha, response_data, dep, tolerance = .025)
calibrate <- function(oob, alpha = .1, response_data, tolerance = .025,
                      step_percent = .618, under = TRUE, method = "quantile",
                      max_iter = 10) {

  #find calibrated alpha
  good_cov <- FALSE
  adjust <- 0
  cov_list <- c()
  alpha_cal <- alpha
  iter <- 0

  #loop until coverage is within tolerance of nominal
  while(good_cov == FALSE) {
    #track iterations
    iter <- iter + 1

    #how should we determine step size?
    alpha_cal <- adjust + alpha

    if(iter > max_iter){
      stop("The maximum number of iterations has been exceeded. Increase tolerance or decrease step_percent.")
    }

    if(!as.logical((alpha_cal < 1)*(alpha_cal > 0))){
      stop("The calibrated alpha is not within the interval (0,1). Adjust parameters.")
    }

    #change to more than just quantile for calibration; fix...
    oob_quant <- unlist(lapply(oob, FUN = quantile, na.rm = TRUE, probs = c(alpha_cal/2, 1 - alpha_cal/2)))
    dim(oob_quant) <- c(2, length(oob))
    oob_quant <- t(oob_quant)

    #coverage with current iteration of calibrated alpha
    cov <- sum((response_data >= oob_quant[,1]) * (response_data <= oob_quant[,2]))/dim(oob_quant)[1]

    #step by ratio of distance
    step <- abs(cov - (1-alpha))*step_percent

    #change based on if under or over
    if(cov < alpha) {adjust <- adjust - step} else {adjust <- adjust + step}

    #must adjust tolerance based on number of samples...
    good_cov <- as.logical((cov <= (1 - alpha) + tolerance) * (cov >= (1 - alpha) - tolerance))
    #print(step)
  }

  #returns calibrate alpha value to achieve nominal coverage, within the proposed tolerance
  return(alpha_cal)
}
