#' Perform inference for the PRICED model using the Bayesian bootstrap for inference
#'
#' Obtain posterior samples based on data fit to
#' the PRICED model, a set of regression models for
#' CLearance, Incidence, and Prevalence, accounting for
#' diagnostic errors.  Bayesian bootstrapping is used to get
#' posterior samples of the regression coefficients.
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
#' @param n_draws integer, number of posterior draws to obtain (not used if method = "normal")
#' @param seed a single value, interpreted as an integer.  Used for replication purposes.
#' @param verbose logical.  Whether to print messages as to where the model fitting algorithm is.
#' @returns  Object of class 'PRICED_boot' which has the following elements:
#' * t0 The regression coefficients using the original sample
#' * t A matrix where each row is a posterior sample and each column is a regression
#' coefficient
#' * R  The number of bootstrap samples
#' * formula_prevalence, formula_incidence, formula_clearance 
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
#'
#' @export
#' @import magrittr
#' @import dplyr
#' @import adaptMCMC
#' @import methods
#' @import coda
#' @import tibble
#' @import parallel
#' @importFrom numDeriv hessian
#' @importFrom bayesboot rudirichlet
#' @exportClass PRICED

PRICED_bboot = function(formula_prevalence,
                        formula_incidence,
                        formula_clearance,
                        data,
                        prior_prevalence = list(location = 0, scale = 2.5, autoscale = TRUE),
                        prior_incidence = list(location = 0, scale = 2.5, autoscale = TRUE),
                        prior_clearance = list(location = 0, scale = 2.5, autoscale = TRUE),
                        sensitivity = 0.75,
                        specificity = 0.98,
                        n_bootstraps = 250,
                        cl,
                        verbose = TRUE){
  
  
  
  PRICED_no_hessian = function(formula_prevalence,
                               formula_incidence,
                               formula_clearance,
                               data,
                               prior_prevalence,
                               prior_incidence,
                               prior_clearance,
                               sensitivity,
                               specificity,
                               verbose
                               ){
    
    # Curating data -----------------------------------------------------------
    
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
    
    # Make inference ----------------------------------------------------------
    
    return_object =
      list(coefficients = NULL)
    
    
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
    
    names(return_object$coefficients) =
      c(paste("Prevalence:",colnames(X$prevalence)),
        paste("Incidence:",colnames(X$incidence)),
        paste("Clearance:",colnames(X$clearance)))
    
    # Return results
    class(return_object) = "PRICED"
    return(return_object)
  }
  
  if(verbose) cat("\n-----------------------------------------------\n")
  if(verbose) cat("Obtaining PRICED estimates from original sample\n")
  if(verbose) cat("-----------------------------------------------\n\n")
  
  orig_fit = 
    PRICED_no_hessian(formula_prevalence = formula_prevalence,
                      formula_incidence = formula_incidence,
                      formula_clearance = formula_clearance,
                      data = data,
                      prior_prevalence = prior_prevalence,
                      prior_incidence = prior_incidence,
                      prior_clearance = prior_clearance,
                      sensitivity = sensitivity,
                      specificity = specificity,
                      verbose = verbose)
  p = length(orig_fit[[1]])
  
  form_char = as.character(formula_clearance)[-1]
  pathogen_var = form_char[1]
  subject_var =
    gsub('[[:punct:]]',"",
         substr(form_char[2],
                start = gregexpr("\\|",form_char[2])[[1]][1] + 2,
                stop = 1e4)
    )
  
  
  if(verbose) cat("\n---------------------------------\n")
  if(verbose) cat("Performing Bayesian bootstrapping\n")
  if(verbose) cat("---------------------------------\n\n")
  
  IDs = unique(data[[subject_var]])
  N = NROW(IDs)
  subject_indices = 
    lapply(IDs,function(id) which(data[[subject_var]] == id))
  T_i = sapply(subject_indices,length)
  
  bs_weights = 
    bayesboot::rudirichlet(n_bootstraps,N)
  
  
  if(missing(cl)){
    
    
    coef_matrix = 
      matrix(NA,n_bootstraps + 1,p)
    coef_matrix[1,] = coef(orig_fit)
    colnames(coef_matrix) = names(coef(orig_fit))
    
    
    if(verbose) pb = txtProgressBar(0,n_bootstraps,style=3)
    for(i in 1:n_bootstraps){
      try({
        resamp_subjects = sample(N,N,T,prob = bs_weights[i,])
        bs_data = 
          data[unlist(subject_indices[resamp_subjects]),] %>% 
          mutate(subject = rep(1:N,T_i[resamp_subjects]))
        bs_fit = 
          PRICED_no_hessian(formula_prevalence = formula_prevalence,
                            formula_incidence = formula_incidence,
                            formula_clearance = formula_clearance,
                            data = bs_data,
                            prior_prevalence = prior_prevalence,
                            prior_incidence = prior_incidence,
                            prior_clearance = prior_clearance,
                            sensitivity = sensitivity,
                            specificity = specificity,
                            verbose = FALSE)
        coef_matrix[1 + i,] = bs_fit$coefficients
      },silent=T)
      
      setTxtProgressBar(pb,i)
    }
    
    
  }else{
    
    wrapper_fun = function(i){
      coefs = rep(NA,p)
      
      try({
        resamp_subjects = sample(N,N,T,prob = bs_weights[i,])
        bs_data = 
          data[unlist(subject_indices[resamp_subjects]),] %>% 
          mutate(subject = rep(1:N,T_i[resamp_subjects]))
        bs_fit = 
          PRICED_no_hessian(formula_prevalence = formula_prevalence,
                            formula_incidence = formula_incidence,
                            formula_clearance = formula_clearance,
                            data = bs_data,
                            prior_prevalence = prior_prevalence,
                            prior_incidence = prior_incidence,
                            prior_clearance = prior_clearance,
                            sensitivity = sensitivity,
                            specificity = specificity,
                            verbose = FALSE)
        coefs = bs_fit$coefficients
      },silent=T)
      
      return(coefs)
    }
    
    clusterEvalQ(cl,{library(PRICED);library(magrittr);library(dplyr)})
    clusterExport(cl,
                  c("formula_prevalence",
                    "formula_incidence",
                    "formula_clearance",
                    "data",
                    "prior_prevalence",
                    "prior_incidence",
                    "prior_clearance",
                    "sensitivity",
                    "specificity",
                    "N",
                    "subject_indices",
                    "T_i",
                    "bs_weights",
                    "PRICED_no_hessian",
                    "p"),
                  envir = environment())
    coef_matrix = 
      rbind(orig_fit$coefficients,
            t(parSapply(cl,1:n_bootstraps,wrapper_fun))
      )
  }
  
  
  
  # Return bootstraps
  
  return_object =
    list(t0 = coef_matrix[1,],
         t = coef_matrix[-1,],
         R = n_bootstraps,
         formula_prevalence = formula_prevalence,
         formula_incidence = formula_incidence,
         formula_clearance = formula_clearance)
  
  class(return_object) = "PRICED_bboot"
  return(return_object)
}
