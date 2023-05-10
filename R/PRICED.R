#' Perform inference for the PRICED model
#'
#' Obtain posterior samples based on data fit to
#' the PRICED model, a set of regression models for
#' CLearance, Incidence, and Prevalence, accounting for
#' diagnostic errors.
#'
#' @param formula_clearance Formula of the form y ~ x_1 + x_2 + ... + (time_variable | subject_ID)
#' @param formula_incidence Right handed formula of the form  ~ x_1 + x_2 + ...
#' @param formula_prevalence Right handed formula of the form  ~ x_1 + x_2 + ...
#' @param data Data frame or tibble
#' @param prior_clearance,prior_incidence,prior_beta_pr a list giving the normal hyperparameters location and scale (i.e., sd).
#' Additionally, autoscale can be set to TRUE, in which case the scale will be
#' divided by the standard deviation of the covariates for each regression coefficient.
#' @param sensitivity Diagnostic test sensitivity
#' @param specificity Diagnostic test specificity
#' @param method character, either "normal" for asymptotic normal approximation of the posterior (default), or
#' "adaptMCMC" to use adaptive MCMC sampling using the adaptMCMC package.
#' @param n_draws integer, number of posterior draws to obtain (not used if method = "normal")
#' @param seed a single value, interpreted as an integer.  Used for replication purposes.
#' @param verbose logical.  Whether to print messages as to where the model fitting algorithm is.
#' @param ... Further arguments passed to adaptMCMC::MCMC() (if method = "adaptMCMC", otherwise ignored)
#' For example, you'll probably want to pass in acc.rate = 0.234.
#' @returns  Object of class 'PRICED' which has the following elements:
#' * beta_C A matrix where each row is a posterior sample and each column is a regression
#' coefficient in the clearance rate model
#' * beta_I A matrix where each row is a posterior sample and each column is a regression
#' coefficient in the incidence rate model
#' * beta_pr  A matrix where each row is a posterior sample and each column is a regression
#' coefficient in the prevalence model
#' * data Original data frame/tibble with the following two columns added: (1) tau,
#' the time between consecutive tests (within subject), and (2) presence_estimate,
#' the posterior mean that the pathogen was present
#' * priors List of prior hyperparameters used in the model fitting
#' * acceptance_rate Acceptance rate for the adaptive Metropolis-Hastings
#'
#' @examples
#' PRICED_data =
#'   simulate_PRICED(seed = 2023,N_subj = 500, N_time = 10)
#'
#' formula_clearance =  p_observed ~ x1 + x2 + (time | subject)
#' formula_clearance_oracle =  p ~ x1 + x2 + (time | subject)
#' formula_incidence =  ~ x1 + x3
#' formula_prevalence =  ~ x1
#'
#' # Fit the oracle model (no diagnostic error introduced)
#' PRICED_fit0 =
#'   PRICED(formula_prevalence,
#'          formula_incidence,
#'          formula_clearance_oracle,
#'          PRICED_data$data[[1]],
#'          prior_prevalence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          prior_incidence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          prior_clearance = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          method = "normal",
#'          seed = 1,
#'          sensitivity = 1,
#'          specificity = 1,
#'          verbose = TRUE)
#'
#' # Fit the PRICED model, accounting for diagnostic error
#' PRICED_fit1 =
#'   PRICED(formula_prevalence,
#'          formula_incidence,
#'          formula_clearance,
#'          PRICED_data$data[[1]],
#'          prior_prevalence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          prior_incidence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          prior_clearance = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          method = "normal",
#'          seed = 1,
#'          sensitivity = 0.75,
#'          specificity = 0.98,
#'          verbose = TRUE)
#'
#' # Fit the model, ignoring diagnostic error
#' PRICED_fit2 =
#'   PRICED(formula_prevalence,
#'          formula_incidence,
#'          formula_clearance,
#'          PRICED_data$data[[1]],
#'          prior_prevalence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          prior_incidence = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          prior_clearance = list(location = 0, scale = 2.5, autoscale = TRUE),
#'          method = "normal",
#'          seed = 1,
#'          sensitivity = 1,
#'          specificity = 1,
#'          verbose = TRUE)
#'
#' # Collate the results
#' library(tibble)
#' collated_results =
#'   tibble(parameter = names(unlist(PRICED_data$parameters[[1]][c(3:1)])),
#'          truth = unlist(PRICED_data$parameters[[1]][3:1]),
#'          `oracle estimate` = coef(PRICED_fit0),
#'          `PRICED estimate` = coef(PRICED_fit1),
#'          `naive estimate` = coef(PRICED_fit2))
#' print(collated_results)
#'
#' #Compute MSE:
#' mean((collated_results$`oracle estimate` - collated_results$truth)^2)
#' mean((collated_results$`PRICED estimate` - collated_results$truth)^2)
#' mean((collated_results$`naive estimate` - collated_results$truth)^2)
#'
#'
#' @export
#' @import magrittr
#' @import dplyr
#' @import adaptMCMC
#' @import methods
#' @import coda
#' @import tibble
#' @importFrom numDeriv hessian
#' @exportClass PRICED


PRICED = function(formula_prevalence,
                  formula_incidence,
                  formula_clearance,
                  data,
                  prior_prevalence = list(location = 0, scale = 2.5, autoscale = TRUE),
                  prior_incidence = list(location = 0, scale = 2.5, autoscale = TRUE),
                  prior_clearance = list(location = 0, scale = 2.5, autoscale = TRUE),
                  sensitivity = 0.75,
                  specificity = 0.98,
                  method = c("normal","adaptMCMC")[1],
                  n_draws = 5e3,
                  seed = NULL,
                  verbose = TRUE,
                  ...){

  if(!is.null(seed)){
    set.seed(seed)
  }else{
    warning("Make sure to set your seed!")
  }


  # Curating data -----------------------------------------------------------

  if(verbose) cat("\nCurating data and determining which responses/covariates are estimable\n")
  # Drop NAs
  data_aug =
    data %>%
    dplyr::select(c(all.vars(formula_prevalence),
                    all.vars(formula_incidence),
                    all.vars(formula_clearance))) %>%
    na.omit()

  # get names of ...
  form_char = as.character(formula_clearance)[-1]
  ## pathogen columns
  pathogen_var = form_char[1]
  ## covariates for clearance
  X_vars = list()
  X_vars$clearance =
    gsub("\ \\+","",
         substr(form_char[2],
                start = 1,
                stop = gregexpr("\\(",form_char[2])[[1]][1] - 3)
    ) %>%
    strsplit(" ")
  X_vars$clearance = X_vars$clearance[[1]]

  ## subject variable
  subject_var =
    gsub('[[:punct:]]',"",
         substr(form_char[2],
                start = gregexpr("\\|",form_char[2])[[1]][1] + 2,
                stop = 1e4)
    )
  ## time variable
  time_var =
    gsub('[[:punct:]]',"",
         substr(form_char[2],
                start = gregexpr("\\(",form_char[2])[[1]][1] + 1,
                stop = gregexpr("\\|",form_char[2])[[1]][1] - 2)
    )

  form_char = as.character(formula_incidence)[-1]
  X_vars$incidence =
    gsub("\ \\+","", form_char) %>%
    strsplit(" ")
  X_vars$incidence = X_vars$incidence[[1]]

  form_char = as.character(formula_prevalence)[-1]
  X_vars$prevalence =
    gsub("\ \\+","",form_char) %>%
    strsplit(" ")
  X_vars$prevalence = X_vars$prevalence[[1]]

  # Store basic quantities
  NT = nrow(data)
  N =
    data %>%
    dplyr::select(all_of(subject_var)) %>%
    unlist() %>%
    unique() %>%
    length()
  Q = lapply(X_vars, function(x) length(x) + 1) # + 1 for the intercept

  # arrange data by subject and by time
  data_aug %<>%
    arrange(across(starts_with(subject_var)),
            across(starts_with(time_var)))

  # Compute time between observations (tau)
  data_aug %<>%
    mutate(time0 = c(TRUE,data_aug[[subject_var]][-1] !=
                       data_aug[[subject_var]][-NT]),
           t_counter = 0L) %>%
    mutate(tau =
             ifelse(!time0,
                    diff(c(data_aug[[time_var]][1],data_aug[[time_var]])),
                    0.0)
    )
  for(i in 1:nrow(data_aug)){
    data_aug$t_counter[i] =
      ifelse(data_aug$time0[i] == TRUE,
             0L,
             data_aug$t_counter[i - 1] + 1L)
  }
  T_max = max(data_aug$t_counter)

  # Get time indices as a list
  t_indices =
    lapply(1 + 0:T_max,function(tt) which(data_aug$t_counter == tt - 1))

  # Get design matrices
  X = list()
  X$clearance =
    model.matrix(as.formula(paste0("~ ",
                                   paste(X_vars$clearance,collapse = "+"))),
                 data = data_aug)
  X$incidence =
    model.matrix(as.formula(paste0("~ ",
                                   paste(X_vars$incidence,collapse = "+"))),
                 data = data_aug)
  X$prevalence =
    model.matrix(as.formula(paste0("~ ",
                                   paste(X_vars$prevalence,collapse = "+"))),
                 data =
                   data_aug %>%
                   filter(t_counter == 0))

  # Get quantities that only depend on data and Se,Sp
  Se_y = sensitivity^data_aug[[pathogen_var]] * (1 - sensitivity)^(1 - data_aug[[pathogen_var]])
  Sp_y = (1 - specificity)^data_aug[[pathogen_var]] * specificity^(1 - data_aug[[pathogen_var]])



  # Set prior hyperparameter values -----------------------------------------

  sd_x = lapply(X,function(x) apply(as.matrix(x[,-1]),2,sd))
  names(sd_x) = names(X)
  prior_list = list(beta_I = list(mean = prior_incidence$location,
                                  sd =
                                    prior_incidence$scale /
                                    c(1,sd_x$incidence)^prior_incidence$autoscale),
                    beta_C = list(mean = prior_clearance$location,
                                  sd =
                                    prior_clearance$scale /
                                    c(1,sd_x$clearance)^prior_clearance$autoscale),
                    beta_pr = list(mean = prior_prevalence$location,
                                   sd =
                                     prior_prevalence$scale /
                                     c(1,sd_x$prevalence)^prior_prevalence$autoscale))
  rm(sd_x)



  # Create likelihood and posterior functions -------------------------------

  loglik_function = function(x){
    # Convert x vector into parameters
    beta_P = x[1:Q$prevalence]
    beta_I = x[Q$prevalence + 1:Q$incidence]
    beta_C = x[Q$prevalence + Q$incidence + 1:Q$clearance]

    # Compute key quantities
    lambda_tau =
      data_aug$tau * exp(drop(X$incidence %*% beta_I))
    eta_tau =
      data_aug$tau * exp(drop(X$incidence %*% beta_C))
    rho =
      1 / (1 + exp(-drop(X$prevalence %*% beta_P)))

    # Set up matrix for recursive relations
    pqr = matrix(0.0,nrow(data_aug),3,
                 dimnames = list(NULL,letters[16:18]))

    # 1. Set $r_{i0}=\rho_i$
    pqr[t_indices[[1]],"r"] = rho

    # 2. Set $q_{i0} = \frac{r_{i0}S_e^{y_{i0}}(1-S_e)^{1-y_{i0}}}{r_{i0}S_e^{y_{i0}}(1-S_e)^{1-y_{i0}} + (1-r_{i0})S_p^{1 - y_{i0}}(1-S_p)^{y_{i0}}}$
    pqr[t_indices[[1]],"q"] =
      pqr[t_indices[[1]],"r"] * Se_y[t_indices[[1]]]
    pqr[t_indices[[1]],"q"] =
      pqr[t_indices[[1]],"q"] /
      ( pqr[t_indices[[1]],"q"] + (1 - pqr[t_indices[[1]],"r"]) * Sp_y[t_indices[[1]]])

    # 3. For $t = 1,\ldots,T_i$,
    for(tt in 1:T_max){
      # 3a.  Set $r_{it} = q_{i(t-1)}e^{-\eta_{it}\tau_{it}} + (1-q_{i(t-1)})\left(1 - e^{-\lambda_{it}\tau_{it}}\right))$
      pqr[t_indices[[1 + tt]],"r"] =
        pqr[t_indices[[tt]],"q"] * exp(- eta_tau[t_indices[[1 + tt]]]) +
        (1 - pqr[t_indices[[tt]],"q"]) * (1 - exp(-lambda_tau[t_indices[[1 + tt]]]))

      # 3b. Set $q_{it} = \frac{r_{it}S_e^{y_{it}}(1-S_e)^{1-y_{it}}}{r_{it}S_e^{y_{it}}(1-S_e)^{1-y_{it}} + (1-r_{it})S_p^{1 - y_{it}}(1-S_p)^{y_{it}}}$
      if(tt < T_max){
        pqr[t_indices[[1 + tt]],"q"] =
          pqr[t_indices[[1 + tt]],"r"] * Se_y[t_indices[[1 + tt]]]
        pqr[t_indices[[1 + tt]],"q"] =
          pqr[t_indices[[1 + tt]],"q"] /
          ( pqr[t_indices[[1 + tt]],"q"] + (1 - pqr[t_indices[[1 + tt]],"r"]) * Sp_y[t_indices[[1 + tt]]])
      }
    }

    # 4. Finally, for $t=0,1,\ldots,T_i$, set $p_{it} =r_{it}S_e + (1 - r_{it})(1-S_p)$
    pqr[,"p"] = pqr[,"r"] * sensitivity + (1 - pqr[,"r"]) * (1 - specificity)


    # return log likelihood
    return( sum( dbinom(data_aug[[pathogen_var]],size = 1, prob = pqr[,"r"], log = TRUE) ) )

  }

  logpost = function(x){
    # Convert x vector into parameters
    beta_P = x[1:Q$prevalence]
    beta_I = x[Q$prevalence + 1:Q$incidence]
    beta_C = x[Q$prevalence + Q$incidence + 1:Q$clearance]

    return(
      loglik_function(x) +
        sum(dnorm(c(beta_P),
                  mean = c(prior_list$beta_pr$mean),
                  sd = c(prior_list$beta_pr$sd),
                  log = TRUE)) +
        sum(dnorm(c(beta_I),
                  mean = c(prior_list$beta_I$mean),
                  sd = c(prior_list$beta_I$sd),
                  log = TRUE)) +
        sum(dnorm(c(beta_C),
                  mean = c(prior_list$beta_C$mean),
                  sd = c(prior_list$beta_C$sd),
                  log = TRUE))
    )
  }

  get_filtered_z = function(x){
    # Convert x vector into parameters
    beta_P = x[1:Q$prevalence]
    beta_I = x[Q$prevalence + 1:Q$incidence]
    beta_C = x[Q$prevalence + Q$incidence + 1:Q$clearance]

    # Compute key quantities
    lambda_tau =
      data_aug$tau * exp(drop(X$incidence %*% beta_I))
    eta_tau =
      data_aug$tau * exp(drop(X$incidence %*% beta_C))
    rho =
      1 / (1 + exp(-drop(X$prevalence %*% beta_P)))

    # Set up matrix for recursive relations
    pqr = matrix(0.0,nrow(data_aug),3,
                 dimnames = list(NULL,letters[16:18]))

    # 1. Set $r_{i0}=\rho_i$
    pqr[t_indices[[1]],"r"] = rho

    # 2. Set $q_{i0} = \frac{r_{i0}S_e^{y_{i0}}(1-S_e)^{1-y_{i0}}}{r_{i0}S_e^{y_{i0}}(1-S_e)^{1-y_{i0}} + (1-r_{i0})S_p^{1 - y_{i0}}(1-S_p)^{y_{i0}}}$
    pqr[t_indices[[1]],"q"] =
      pqr[t_indices[[1]],"r"] * Se_y[t_indices[[1]]]
    pqr[t_indices[[1]],"q"] =
      pqr[t_indices[[1]],"q"] /
      ( pqr[t_indices[[1]],"q"] + (1 - pqr[t_indices[[1]],"r"]) * Sp_y[t_indices[[1]]])

    # 3. For $t = 1,\ldots,T_i$,
    for(tt in 1:T_max){
      # 3a.  Set $r_{it} = q_{i(t-1)}e^{-\eta_{it}\tau_{it}} + (1-q_{i(t-1)})\left(1 - e^{-\lambda_{it}\tau_{it}}\right))$
      pqr[t_indices[[1 + tt]],"r"] =
        pqr[t_indices[[tt]],"q"] * exp(- eta_tau[t_indices[[1 + tt]]]) +
        (1 - pqr[t_indices[[tt]],"q"]) * (1 - exp(-lambda_tau[t_indices[[1 + tt]]]))

      # 3b. Set $q_{it} = \frac{r_{it}S_e^{y_{it}}(1-S_e)^{1-y_{it}}}{r_{it}S_e^{y_{it}}(1-S_e)^{1-y_{it}} + (1-r_{it})S_p^{1 - y_{it}}(1-S_p)^{y_{it}}}$
      if(tt < T_max){
        pqr[t_indices[[1 + tt]],"q"] =
          pqr[t_indices[[1 + tt]],"r"] * Se_y[t_indices[[1 + tt]]]
        pqr[t_indices[[1 + tt]],"q"] =
          pqr[t_indices[[1 + tt]],"q"] /
          ( pqr[t_indices[[1 + tt]],"q"] + (1 - pqr[t_indices[[1 + tt]],"r"]) * Sp_y[t_indices[[1 + tt]]])
      }
    }

    # 4. Finally, for $t=0,1,\ldots,T_i$, set $p_{it} =r_{it}S_e + (1 - r_{it})(1-S_p)$
    pqr[,"p"] = pqr[,"r"] * sensitivity + (1 - pqr[,"r"]) * (1 - specificity)


    # return log likelihood
    return( pqr[,"q"] )

  }


  # Make inference ----------------------------------------------------------

  return_object =
    list(coefficients = NULL,
         covariance = NULL,
         samples = NULL,
         acceptance_rate = NULL,
         formula_prevalence = formula_prevalence,
         formula_incidence = formula_incidence,
         formula_clearance = formula_clearance,
         priors = prior_list)


  if(method == "normal"){

    # Find posterior mode
    if(verbose) cat("\nFinding posterior mode\n")
    opt =
      optim(par = rnorm(sum(unlist(Q))),
            fn = logpost,
            method = "Nelder-Mead",
            control = list(fnscale = -1,
                           maxit = 5e4))

    return_object$coefficients =
      opt$par

    # Find hessian of the log posterior
    if(verbose) cat("\nFinding hessian for covariance matrix\n")
    opt_hess =
      numDeriv::hessian(logpost,opt$par)

    # Flip to get posterior covariance matrix
    return_object$covariance =
      qr.solve(-opt_hess)
  }

  if(method == "adaptMCMC"){

    if(verbose) cat("\nPerforming adaptive MCMC\n")
    samp =
      MCMC(p = logpost,
           n = n_draws + round(n_draws * 0.1),
           init = rnorm(sum(unlist(Q))),
           ...)

    return_object$samples = samp$samples[-c(1:round(n_draws * 0.1)),]
    return_object$coefficients = colMeans(return_object$samples)
    return_object$acceptance_rate = samp$acceptance.rate
  }

  names(return_object$coefficients) =
    c(paste("Prevalence:",colnames(X$prevalence)),
      paste("Incidence:",colnames(X$incidence)),
      paste("Clearance:",colnames(X$clearance)))


  # Get final filtered estimate of Z
  return_object$data =
    data_aug %>%
    dplyr::select(all_of(c(subject_var,time_var,"tau",pathogen_var,
                           unique(c(X_vars$prevalence,X_vars$clearance,X_vars$incidence))))) %>%
    mutate(estimated_pathogen_presence = get_filtered_z(return_object$coefficients))


  # Return results
  class(return_object) = "PRICED"
  return(return_object)
}
