#' pirf()
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
#' @seealso \link[ranger]{ranger}
#' @seealso \link[rfinterval]{rfinterval}
#' @seealso \link[piRF]{predict.pirf}
pirf <- function(formula = NULL,
                 train_data = NULL,
                 method = "oob",
                 alpha = 0.1,
                 symmetry = TRUE,
                 seed = NULL,
                 m_try = 2,
                 num_trees = 500,
                 min_node_size = NULL,
                 num_threads = parallel::detectCores(),
                 calibrate = FALSE,
                 split = .5,
                 Roy_method = "quantile",
                 featureBias = FALSE,
                 predictionBias = TRUE,
                 Tung_R = 5,
                 Tung_num_trees = 75,
                 variant = 1,
                 Ghosal_num_stages = 2,
                 prop = .618,
                 concise = TRUE,
                 interval_type = "two-sided",
                 ...){

  #check for currently running parallel processes
  if(is.null(num_threads)){
    num_threads <- parallel::detectCores()
  } else {
    num_threads = num_threads
  }

  #fix this
  if(foreach::getDoParWorkers() != num_threads){
    clust <- parallel::makeCluster(num_threads)
    doParallel::registerDoParallel(clust)
  } else {

  }

  #list for all intervals; just gets all output for each interval
  res <- list()
  int <- list()
  preds <- list()
  newpreds <- list()

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

  #browser()
  #tracking list
  for(id in method){
    if(id == "oob"){
      #oob call; no one sided

      #res[[id]] <- ranger::ranger(formula = formula, data = train_data,
      #                            num.trees = num_trees, mtry = m_try, min.node.size = min_node_size,
      #                            num.threads = num_threads)

      res[[id]] <- rfinterval::rfinterval(formula = formula, train_data = train_data, test_data = train_data,
                                          method = "oob", alpha = alpha, symmetry = symmetry, seed = seed,
                                          params_ranger = list(mtry = m_try, num.trees = num_trees, min.node.size = min_node_size,
                                                               num.threads = num_threads))

      #predictions for training data
      preds[[id]] <- res[[id]]$testPred
      int[[id]] <- res[[id]]$oob_interval
    } else if (id == "bop"){
      #bop

      #new fit function
      res[[id]] <- fit_bop(formula = formula,
                           train_data = train_data,
                           num_trees = num_trees,
                           num_threads = num_threads,
                           m_try = m_try)

      #new predict function
      newpreds[[id]] <- predict_bop(res[[id]],
                                    pred_data = train_data,
                                    interval_method = Roy_method,
                                    calibrate = calibrate,
                                    alpha = alpha,
                                    tolerance = tolerance,
                                    interval_type = interval_type,
                                    num_threads = num_threads)

      #predictions for training data
      preds[[id]] <- newpreds[[id]]$preds
      int[[id]] <- newpreds[[id]]$pred_intervals

    } else if (id == "brf"){
      #brf

      #new fit function
      res[[id]] <- fit_brf(formula,
                           train_data = train_data,
                           num_trees = num_trees,
                           min_node_size = min_node_size,
                           m_try = m_try,
                           keep_inbag = TRUE,
                           prop = prop,
                           variant = variant,
                           num_threads = num_threads,
                           num_stages = Ghosal_num_stages)

      #new predict function
      newpreds[[id]] <- predict_brf(res[[id]],
                                    pred_data = train_data,
                                    alpha = alpha,
                                    interval_type = interval_type,
                                    variant = variant)

      res[[id]]$mse_est <- newpreds[[id]]$mse_est

      #predictions for training data
      preds[[id]] <- newpreds[[id]]$preds
      int[[id]] <- newpreds[[id]]$pred_intervals

    } else if (id == "bcqrf"){
      #bcqrf

      #new fit function
      res[[id]] <- fit_bcqrf(formula = formula,
                             train_data = train_data,
                             num_trees = num_trees,
                             min_node_size = min_node_size,
                             m_try = m_try,
                             keep_inbag = TRUE,
                             feature_num_trees = Tung_num_trees,
                             featureBias = TRUE,
                             predictionBias = TRUE,
                             R = Tung_R,
                             num_threads = num_threads)

      #new predict function
      newpreds[[id]] <- predict_bcqrf(res[[id]],
                                      train_data = train_data,
                                      pred_data = train_data,
                                      intervals = TRUE,
                                      alpha = alpha,
                                      num_threads = num_threads,
                                      interval_type = interval_type)

      #predictions for training data
      preds[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals

    } else if (id == "cqrf"){
      #cqrf

      #new fit function
      res[[id]] <- fit_cqrf(formula = formula,
                            train_data = train_data,
                            split = split,
                            num_trees = num_trees,
                            min_node_size = min_node_size,
                            m_try = m_try,
                            keep_inbag = TRUE,
                            forest_type = "RF",
                            num_threads = num_threads)

      #new predict function
      newpreds[[id]] <- predict_cqrf(res[[id]],
                                     pred_data = train_data,
                                     alpha = alpha,
                                     num_threads = num_threads,
                                     intervals = TRUE,
                                     interval_type = "two-sided")

      #predictions for training data
      preds[[id]] <- res[[id]]$preds
      int[[id]] <- res[[id]]$pred_intervals

    } else if (id == "hdi"){
      #hdi

      #new fit function
      res[[id]] <- fit_hdi(formula,
                           train_data = train_data,
                           num_tree = num_trees,
                           mtry = m_try,
                           min_node_size = min_node_size,
                           max_depth = 10,
                           replace = TRUE,
                           verbose = FALSE,
                           num_threads = num_threads)

      #new predict function
      newpreds[[id]] <- predict_hdi(res[[id]],
                                    train_data,
                                    alpha = alpha,
                                    num_threads = num_threads)

      #predictions for training data
      preds[[id]] <- newpreds[[id]]$preds
      int[[id]] <- newpreds[[id]]$pred_intervals

    } else if (id == "quantile"){
      #qrf call
      #intermediate step to get quantReg model
      res[[id]] <- ranger::ranger(formula = formula, data = train_data,
                                  num.trees = num_trees, mtry = m_try, min.node.size = min_node_size,
                                  quantreg = TRUE, num.threads = num_threads)

      #need to do this with predict for quantReg with ranger
      preds[[id]] <- predict(res[[id]], train_data, type = "quantiles", quantiles = c(alpha1, 0.5, alpha2))$predictions
      int[[id]] <- preds[[id]][,c(1,3)]

      #save only median
      preds[[id]] <- preds[[id]][,2]
    } else {
      stop(paste0(id, " is not a supported random forest prediction interval methodology"))
    }

    #one sided intervals; adjusts intervals based on on what sided interval is desired
    #if(interval_type == "upper"){
    #  int[[id]][,1] <- -Inf
    #} else if(interval_type == "lower"){
    #  int[[id]][,2] <- Inf
    #}

    #colnames(int[[id]]) <- c("lower", "upper")
  }

  res$call <- formula
  res$train_data <- train_data

  foreach::registerDoSEQ()

  class(res) = "pirf"
  return(res)
}

#' predict.pirf()
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
predict.pirf <- function(res,
                         test_data = NULL,
                         num_threads = parallel::detectCores(),
                         concise = TRUE,
                         interval_type = "two-sided",
                         alpha = 0.1,
                         symmetry = TRUE,
                         seed = NULL,
                         m_try = 2,
                         num_trees = 500,
                         min_node_size = NULL,
                         calibrate = FALSE,
                         tolerance = .01,
                         split = .5,
                         Roy_method = "quantile",
                         featureBias = FALSE,
                         predictionBias = TRUE,
                         Tung_R = 5,
                         Tung_num_trees = 75,
                         variant = 1,
                         Ghosal_num_stages = 2,
                         prop = .618,
                         ...){

  #browser()

  #list of models fit
  method <- names(res)[!(names(res) %in% c("call", "train_data"))]

  #check for currently running parallel processes
  if(is.null(num_threads)){
    num_threads <- parallel::detectCores()
  } else {
    num_threads = num_threads
  }

  #fix this
  if(foreach::getDoParWorkers() != num_threads){
    clust <- parallel::makeCluster(num_threads)
    doParallel::registerDoParallel(clust)
  } else {

  }

  #list for all intervals; just gets all output for each interval
  #res <- list()
  int <- list()
  preds <- list()
  newpreds <- list()

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

  #tracking list
  train_data <- res$train_data
  formula <- as.formula(res$call)
  for(id in method){
    if(id == "oob"){
      #oob

      #fix this; retrains a new model for every new set of test predictions
      res[[id]] <- rfinterval::rfinterval(formula = formula, train_data = train_data, test_data = test_data,
                                          method = "oob", alpha = alpha, symmetry = symmetry, seed = seed,
                                          params_ranger = list(num.threads = num_threads))

      #predictions for test data
      preds[[id]] <- res[[id]]$testPred
      int[[id]] <- res[[id]]$oob_interval

    } else if (id == "bop"){
      #bop

      #new predict function
      newpreds[[id]] <- predict_bop(res[[id]],
                                    pred_data = test_data,
                                    interval_method = Roy_method,
                                    calibrate = calibrate,
                                    alpha = alpha,
                                    tolerance = tolerance,
                                    interval_type = interval_type,
                                    num_threads = num_threads)

      preds[[id]] <- newpreds[[id]]$preds
      int[[id]] <- newpreds[[id]]$pred_intervals

    } else if (id == "brf"){
      #brf

      #new predict function
      newpreds[[id]] <- predict_brf(res[[id]],
                                    pred_data = test_data,
                                    alpha = alpha,
                                    interval_type = interval_type,
                                    variant = variant)

      preds[[id]] <- newpreds[[id]]$preds
      int[[id]] <- newpreds[[id]]$pred_intervals

    } else if (id == "bcqrf"){
      #bcqrf

      #browser()

      #new predict function
      newpreds[[id]] <- predict_bcqrf(res[[id]],
                                      train_data = train_data,
                                      pred_data = test_data,
                                      intervals = TRUE,
                                      alpha = alpha,
                                      num_threads = num_threads,
                                      interval_type = interval_type)

      preds[[id]] <- newpreds[[id]]$preds
      int[[id]] <- newpreds[[id]]$pred_intervals

    } else if (id == "cqrf"){
      #cqrf

      #new predict function
      newpreds[[id]] <- predict_cqrf(res[[id]],
                                     pred_data = test_data,
                                     alpha = alpha,
                                     num_threads = num_threads,
                                     intervals = TRUE,
                                     interval_type = "two-sided")

      preds[[id]] <- newpreds[[id]]$preds
      int[[id]] <- newpreds[[id]]$pred_intervals

    } else if (id == "hdi"){
      #hdi

      #new predict function
      newpreds[[id]] <- predict_hdi(res[[id]],
                                    test_data,
                                    alpha = alpha,
                                    num_threads = num_threads)

      preds[[id]] <- newpreds[[id]]$preds
      int[[id]] <- newpreds[[id]]$pred_intervals

    } else if (id == "quantile"){
      #qrf call

      #need to do this with predict for quantReg with ranger
      preds[[id]] <- predict(res[[id]], test_data, type = "quantiles", quantiles = c(alpha1, 0.5, alpha2))$predictions
      int[[id]] <- preds[[id]][,c(1,3)]

      #save only median
      preds[[id]] <- preds[[id]][,2]
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

  return(list(method = method, preds = preds, int = int))
  foreach::registerDoSEQ()

}
