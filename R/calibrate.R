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

#' calibrate()
#'
#' This function outputs a calibrated significance level based on coverage of prediction intervals generated using oob collections. Primarily for use in RoyRF(). Attempting to see which other methods could utilize this procedure.
#' @param oob collection of oob predictions for training data (in list form).
#' @param alpha nominal significance level. Defaults to 0.01.
#' @param response_data response data of class data.frame. Must have names() attribute.
#' @param tolerance tolerance allowed around nominal alpha. Default is 0.25.
#' @param step_ratio ratio absolute difference between empirical oob coverage and nominal coverage to adjust when calibrating. Defaults to 0.618.
#' @param undercoverage Allow undercoverage. Defaults to TRUE. Not currently implemented.
#' @param method Method to calibrate prediction intervals with. Defaults to "quantile"). Current only "quantile" implemented.
#' @param max_iter Maximum number of iterations. Defaults to 10.
#' @keywords random forest, calibration
#' @export
#' @examples
#' calibrate <- function(oob, alpha = alpha, response_data, dep, tolerance = .025)
#' @noRd
calibrate <- function(oob, alpha = .1, response_data, tolerance = .025,
                      step_percent = .618, undercoverage = FALSE, method = "quantile",
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

    #make it so that alpha_cal > 0
    if(!as.logical((alpha_cal < 1)*(alpha_cal > 0))){
      stop("The calibrated alpha is not within the interval (0,1). Adjust parameters.")
    }

    #currently calibration is based on quantile prediction intervals
    #which quantile function should we use?
    oob_quant <- unlist(lapply(oob, FUN = quantile, na.rm = TRUE, probs = c(alpha_cal/2, 1 - alpha_cal/2)))
    dim(oob_quant) <- c(2, length(oob))
    oob_quant <- t(oob_quant)

    #coverage with current iteration of calibrated alpha
    cov <- sum((response_data >= oob_quant[,1]) * (response_data <= oob_quant[,2]))/dim(oob_quant)[1]

    #step by ratio of distance
    step <- abs(cov - (1-alpha))*step_percent

    #change based on if undercoverage or overcoverage
    if(cov < alpha) {adjust <- adjust - step} else {adjust <- adjust + step}

    #must adjust tolerance based on number of samples...
    good_cov <- as.logical((cov <= (1 - alpha) + tolerance) * (cov >= (1 - alpha) - tolerance))
  }

  #returns calibrate alpha value to achieve nominal coverage, within the proposed tolerance
  return(alpha_cal)
}
