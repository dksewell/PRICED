#' Simulate data from the PRICED model
#'
#' @param nsim integer. The number of simulated data sets to create.
#' @param seed a single value, interpreted as an integer.  Used for replication purposes.
#' @param N_subj positive integer. Number of subjects in data sets
#' @param N_time Either positive integer or vector of integers of length N_subj.
#' Number of time points per subject.
#' @param bl_prevalence positive numeric, giving the
#' baseline prevalences for each pathogen type.
#' @param bl_clearance_rate positive numeric, giving the overall rate of clearance
#' @param sensitivity numeric, giving diagnostic test sensitivity
#' @param specificity numeric, giving diagnostic test specificity
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

simulate_PRICED = function(nsim = 1,
                           seed = NULL,
                           N_subj = 250,
                           N_time = 5,
                           bl_prevalence = 0.05,
                           bl_clearance_rate = 1/3,
                           sensitivity = 0.75,
                           specificity = 0.98,
                           beta_C_slopes = c(0,1),
                           beta_I_slopes = c(-1.5,0.5),
                           beta_pr_slope = 1.5){

  # Set seed for reproducibility
  if(!is.null(seed)){
    set.seed(seed)
  }else{
    warning("Make sure to set your seed!")
  }

  # Set the number of time points for each subject, if not done already
  if(length(N_time) == 1) N_time = rep(N_time,N_subj)

  # Compute incidence rate
  bl_incidence_rate = bl_prevalence / (1 - bl_prevalence) / bl_clearance_rate


  # Set up objects to be returned
  parameters = sic_data = list()


  # Run nsim simulations
  for(it in 1:nsim){

    # Set up tibble with subject and time information
    sic_data[[it]] =
      tibble(subject = rep(1:N_subj,N_time),
             time = rgamma(sum(N_time),
                           shape = 2,
                           rate = 2)) %>%
      arrange(subject,
              time)

    # Compute time between observations
    sic_data[[it]] %<>%
      mutate(tau =
               ifelse(
                 c(TRUE,sic_data[[it]]$subject[-1] == sic_data[[it]]$subject[-nrow(sic_data[[it]])]),
                 diff(c(sic_data[[it]]$time[1],sic_data[[it]]$time)),
                 0.0)
      )

    # Generate one binary and two continuous covariates
    sic_data[[it]] %<>%
      mutate(x1 = rbinom(nrow(sic_data[[it]]),
                         1,
                         0.25),
             x2 = 1 + rnorm(nrow(sic_data[[it]])),
             x3 = rnorm(nrow(sic_data[[it]])))

    # Set up parameter values
    beta_C = numeric(3)
    beta_I = numeric(3)
    beta_pr = numeric(2)

    beta_C[1] = log(bl_clearance_rate) - sum(beta_C_slopes)
    beta_C[-1] = beta_C_slopes
    beta_I[1] = log(bl_incidence_rate)
    beta_I[-1] = beta_I_slopes
    beta_pr[1] = log(bl_prevalence / (1 - bl_prevalence))
    beta_pr[2] = beta_pr_slope

    parameters[[it]] =
      list(beta_C = beta_C,
           beta_I = beta_I,
           beta_pr = beta_pr,
           sensitivity = sensitivity,
           specificity = specificity)

    # Create design matrix
    X = list()
    X$clearance =
      sic_data[[it]] %>%
      select(x1,x2) %>%
      as.matrix()
    X$incidence =
      sic_data[[it]] %>%
      select(x1,x3) %>%
      as.matrix()
    X$prevalence =
      sic_data[[it]] %>%
      select(x1) %>%
      as.matrix()
    X %<>% lapply(function(x) cbind(`(Intercept)` = 1,x))



    # Initialize y using baseline prevalences
    ind = which(sic_data[[it]]$tau == 0)
    sic_data[[it]] %<>%
      mutate(p = 0L)
    sic_data[[it]]$p[ind] =
      rbinom(length(ind),
             1,
             1 / (1 + exp(-X$prevalence[ind,] %*% beta_pr)))


    # Simulate data for each subject
    for(i in 1:N_subj){
      # Pull out appropriate rows of sic_data[[it]]
      i_index = which(sic_data[[it]]$subject == i)

      # Find out what their initial infection statuses are
      y_tm1 = sic_data[[it]]$p[i_index[1]]

      # Simulate over each time point
      for(tt in 2:length(i_index)){

        # Create rate
        lambda_it =
          (1 - y_tm1) *
          (
            X$incidence[i_index[tt],] %*% beta_I
          ) +
          y_tm1 *
          (
            X$clearance[i_index[tt],] %*% beta_C
          )
        lambda_it = exp(lambda_it)

        # Compute the probability of detecting each pathogen.  Will
        # vary based on the pathogen presence previously
        prob1 =
          ifelse(y_tm1,
                 exp(-sic_data[[it]]$tau[i_index[tt]] *
                       lambda_it),
                 1 - exp(-sic_data[[it]]$tau[i_index[tt]] *
                           lambda_it))

        # Simulate new data at time tt
        y_tm1 = # Yes, this is actually y_t, but it will serve next as y_{t-1}
          rbinom(1,1,prob1)
        sic_data[[it]]$p[i_index[tt]] = y_tm1
      }

    }


    # Include diagnostic testing error
    sic_data[[it]]$p_observed =
      ifelse(sic_data[[it]]$p == 1,
             rbinom(nrow(sic_data[[it]]),1,sensitivity),
             rbinom(nrow(sic_data[[it]]),1,1 - specificity))

  }

  # Return data and stochastically generated parameters
  return(list(data = sic_data,
              parameters = parameters))
}
