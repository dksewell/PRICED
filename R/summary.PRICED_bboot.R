#' Summarizing PRICED fits
#'
#' summary method for class "PRICED"
#'
#' @param object an object of class "PRICED_boot"
#' @param CI_level level for credible intervals
#' @param estimator character, either "mean", "median", or "mode" for respective posterior point estimate.
#' @param CI_method character, either "HPD" (highest posterior density) or "central".
#' @param print logical
#'
#' @examples
#' \dontrun{
#' PRICED_data =
#'   simulate_PRICED(seed = 1,N_subj = 500, N_time = 10)
#' cl = makeCluster(detectCores())
#' PRICED_fit_parallel = PRICED_bboot(formula_prevalence = ~ x1 + x2 + x3,
#'                                    formula_incidence = ~ x1 + x2 + x3,
#'                                    formula_clearance = p_observed ~ x1 + x2 + x3 + (time | subject),
#'                                    data = PRICED_data$data[[1]],
#'                                    prior_prevalence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'                                    prior_incidence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'                                    prior_clearance = list(location = 0, scale = 2.5, autoscale = TRUE),
#'                                    sensitivity = 0.75,
#'                                    specificity = 0.99,
#'                                    n_bootstraps = detectCores() * 20,
#'                                    cl,
#'                                    verbose = TRUE)
#' stopCluster(cl)
#' summary(PRICED_fit_parallel)
#' }
#' 
#' @export

summary.PRICED_bboot = function(object,
                                CI_level = 0.95,
                                estimator = c("mean","median","mode"),
                                CI_method = c("HPD","central")[1],
                                print = TRUE,
                                burnin = 0){

  a = 1 - CI_level

  summ =
    tibble(Parameter = colnames(object$t))
  
  
  if("mean" %in% estimator){
    summ %<>%
      mutate(`Posterior mean` =
               colMeans(object$t))
  }
  if("median" %in% estimator){
    summ %<>%
      mutate(`Posterior median` =
               apply(object$t,2,median))
  }
  if("mode" %in% estimator){
    summ %<>%
      mutate(`Posterior mode` = object$t0)
  }

  if(CI_method == "central"){
    temp =
      rbind(
        t(apply(object$t,2,quantile,probs = c(a/2,1-a/2)))
      )
    summ %<>%
      mutate(lower = temp[,1],
             upper = temp[,2])
  }else{
    temp =
      coda::HPDinterval(as.mcmc(object$t),CI_level)
    summ %<>%
      mutate(lower = temp[,1],
             upper = temp[,2])
  }
  
  
  colnames(summ)[which(colnames(summ) == "lower")] =
    paste0(a / 2 * 100,"%")
  colnames(summ)[which(colnames(summ) == "upper")] =
    paste0((1 - a / 2) * 100,"%")

  if(print) print(summ,n = Inf)

  invisible(summ)
}
