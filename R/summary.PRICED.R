#' Summarizing PRICED fits
#'
#' summary method for class "PRICED"
#'
#' @param object an object of class "PRICED"
#' @param CI_level level for credible intervals
#' @param estimator character, either "mean" or "median" for respective posterior point estimate. Ignored if method = "normal"
#' @param CI_method character, either "HPD" (highest posterior density) or "central".  Ignored if method = "normal"
#' @param print logical
#'
#' @examples
#' PRICED_data =
#' simulate_PRICED(seed = 2023,N_subj = 500, N_time = 10)
#'
#' PRICED_fit1 =
#'   PRICED(formula_prevalence = ~ x1,
#'          formula_incidence = ~ x1 + x3,
#'          formula_clearance = p_observed ~ x1 + x2 + (time | subject),
#'          PRICED_data$data[[1]],
#'          prior_prevalence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          prior_incidence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          prior_clearance = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          method = "normal",
#'          seed = 1,
#'          sensitivity = 0.75,
#'          specificity = 0.98,
#'          verbose = TRUE)
#' summary(PRICED_fit1)
#'
#' set.seed(123)
#' PRICED_fit2 =
#'   PRICED(formula_prevalence = ~ x1,
#'          formula_incidence = ~ x1 + x3,
#'          formula_clearance = p_observed ~ x1 + x2 + (time | subject),
#'          PRICED_data$data[[1]],
#'          prior_prevalence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          prior_incidence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          prior_clearance = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          method = "adaptMCMC",
#'          n_draws = 5e3,
#'          seed = 1,
#'          sensitivity = 0.75,
#'          specificity = 0.98,
#'          verbose = TRUE,
#'          acc.rate = 0.234)
#' summary(PRICED_fit2)
#' @export

summary.PRICED = function(object,
                          CI_level = 0.95,
                          estimator = c("mean","median"),
                          CI_method = c("HPD","central")[1],
                          print = TRUE,
                          burnin = 0){

  a = 1 - CI_level

  summ =
    tibble(Parameter = names(coef(object)))

  if(is.null(object$samples)){

    summ %<>%
      mutate(Estimate = coef(object)) %>%
      mutate(lower = Estimate - qnorm(1 - a/2) * sqrt(diag(object$covariance)),
             upper = Estimate + qnorm(1 - a/2) * sqrt(diag(object$covariance))) %>%
      rename(`Posterior mode` = Estimate)

  }else{

    if("mean" %in% estimator){
      summ %<>%
        mutate(`Posterior mean` =
                 colMeans(object$samples[(burnin+1):nrow(object$samples),]))
    }
    if("median" %in% estimator){
      summ %<>%
        mutate(`Posterior median` =
                 apply(object$samples[(burnin+1):nrow(object$samples),],2,median))
    }

    if(CI_method == "central"){
      temp =
        rbind(
          t(apply(object$samples[(burnin+1):nrow(object$samples),],2,quantile,probs = c(a/2,1-a/2)))
        )
      summ %<>%
        mutate(lower = temp[,1],
               upper = temp[,2])
    }else{
      temp =
        coda::HPDinterval(as.mcmc(object$samples[(burnin+1):nrow(object$samples),]),CI_level)
      summ %<>%
        mutate(lower = temp[,1],
               upper = temp[,2])
    }
  }


  colnames(summ)[which(colnames(summ) == "lower")] =
    paste0(a / 2 * 100,"%")
  colnames(summ)[which(colnames(summ) == "upper")] =
    paste0((1 - a / 2) * 100,"%")

  if(print) print(summ,n = Inf)

  invisible(summ)
}
