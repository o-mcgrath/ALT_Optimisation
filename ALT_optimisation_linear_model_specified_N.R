#########################################################################
### Test plan optimisation for Weibull linear AFT model (Specified N) ###
#########################################################################
#The user specifies initial model parameters (following a preliminary experiment and/or SME input)
#The user specifies the required number of stresses (N)
#The optimisation algorithm returns the optimal stress levels and associated number of test units

###################################################
### Load DEoptim package and required functions ###
###################################################

### Install DEoptim package and load DEoptim library ###
#install.packages("DEoptim") 
library(DEoptim) 

### Function 1. Weibull log likelihood function (linear AFT model) ###
loglikfnc <- function(theta, data) 
{
  b0 <- theta[1]
  b1 <- theta[2]
  a0 <- theta[3] 
  S <- data$stress
  n <- nrow(data)
  ti <- data$lifetime
  di <- data$cens
  mu <- b0 + b1*S
  lam <- exp(mu)
  gam <- exp(a0) 
  
  epsilon <- .Machine$double.eps
  lam[lam == 0] <- epsilon
  ti[ti == 0] <- epsilon
  gam[gam == 0] <- epsilon 
  
  if(any(!is.finite((ti/lam)^gam))) 
  {
    return(10^10) 
  }
  
  loglik <- sum(di*log(gam) + di*(gam - 1)*log(ti) - di*gam*log(lam) - (ti/lam)^gam)
  
  if(!is.finite(loglik)) 
  {
    return(10^10) 
  }
  
  dl_db0 <- sum(-di*gam + gam*((ti/lam)^gam)) 
  dl_db1 <- sum(-di*gam*S + gam*S*(ti/lam)^gam)
  dl_da0 <- sum(di + (log(ti) - log(lam))*di*gam - log(ti/lam)*gam*(ti/lam)^gam)
  grad_vec <- -c(dl_db0, dl_db1, dl_da0) 
  
  if(any(!is.finite(grad_vec))) 
  {
    return(10^10) 
  }
  
  attr(loglik, "gradient") <- grad_vec
  return(-1*loglik) 
}

### Function 2. Test plan objective function ###
objective_fnc <- function(x, n_sim = 1000, test_units = 100, test_duration = 1*24*365, min_units = 1, 
                          stress_distance = 10^(-5), max_stress = 0.9, design_stress = 0.05,
                          p_cens_limit = 0.8, filter_limit = 0.95, 
                          pquant = 0.5, model_params = c(13, -19, 2), dp = 5) 
{
  DEparams <- as.integer(length(x))
  design_pts <- as.integer(DEparams/2) 
  s1_d <- round(x[1:design_pts], dp) 
  s_UB <- seq(from = max_stress - (stress_distance*(design_pts-2)), to = max_stress, by = stress_distance)
  
  s <- numeric(design_pts)
  s[1] <- s1_d[1]
  for (i in 2:design_pts) 
  {
    s[i] <- min(s_UB[i-1], s[i-1] + s1_d[i])
  }
  
  zi <- x[(design_pts + 1):(2*design_pts)]
  p <- exp(zi)/sum(exp(zi))
  n_free <- test_units - (min_units*design_pts) 
  n <- round(n_free*p) + min_units
  diff <- test_units - sum(n)
  n[which.max(p)] <- n[which.max(p)] + diff
  
  B0 <- model_params[1]
  B1 <- model_params[2]
  gamma <- model_params[3]
  stress <- rep(s, times = n)
  lambda <- exp(B0 + B1*stress)
  true_quantile <- exp(B0 + B1*design_stress)*(-log(1 - pquant))^(1/gamma)
  
  est_quantile <- numeric(n_sim)
  lifetime <- cens <- numeric(test_units)
  i <- 1
  while (i <= n_sim) 
  {
    lifetime[] <- lambda*(-log(1 - runif(n = test_units)))^(1/gamma) 
    cens[] <- ifelse(lifetime >= test_duration, 0, 1)
    lifetime[] <- pmin(lifetime, test_duration)
    dataset <- data.frame(stress = stress, lifetime = lifetime, cens = cens)
    
    if (mean(cens == 0) > p_cens_limit) 
    {
      return(10^10) 
    }
    
    epsilon <- .Machine$double.eps
    no_params <- length(model_params)
    lam_exp <- -log(sum(dataset$cens + epsilon)/sum(dataset$lifetime + epsilon))
    exp_mle <- c(lam_exp, rep(0.01, times = (no_params - 1)))
    exp_small_jitter <- exp_mle + runif(no_params, min = -0.2, max = 0.2)
    exp_large_jitter <- exp_mle + runif(no_params, min = -1, max = 1)
    init_values <- list(exp_mle, exp_small_jitter, exp_large_jitter)
    
    est_quantile[i] <- NA
    
    for (j in 1:length(init_values)) 
    {
      result <- try(nlm(loglikfnc, p = init_values[[j]], data = dataset, iterlim = 1000, check.analyticals = FALSE), silent = TRUE)
      
      if (!inherits(result, "try-error") && is.finite(result$minimum)) 
      {
        mle_b0 <- result$estimate[1]
        mle_b1 <- result$estimate[2]
        mle_gamma <- exp(result$estimate[3])
        if (is.finite(mle_gamma)) 
        {
          est_quantile[i] <- exp(mle_b0 + mle_b1*design_stress)*(-log(1 - pquant))^(1/mle_gamma)
          break
        }
      }
    }
    
    i <- i + 1
  }
  
  est_quantile <- est_quantile[is.finite(est_quantile) & est_quantile > 0]
  Q1 <- quantile(est_quantile, probs = 0.25)
  Q3 <- quantile(est_quantile, probs = 0.75)
  IQR <- Q3 - Q1
  est_quantile <- est_quantile[est_quantile >= Q1 - 6*IQR & est_quantile <= Q3 + 6*IQR]
  s2 <- var(est_quantile)
  bias <- mean(est_quantile) - true_quantile
  
  if (length(est_quantile) <= n_sim*filter_limit || !is.finite(s2 + bias^2)) 
  {
    return(10^10) 
  }else 
    {
      return(sqrt(s2 + bias^2))  
    }
}

### Function 3. Optimisation results function ###
results_fnc <- function(x, test_units = 100, min_units = 1, stress_distance = 10^(-5), max_stress = 0.9, dp = 5, rmse) 
{
  DEparams <- as.integer(length(x))
  design_pts <- as.integer(DEparams/2)
  s1_d <- round(x[1:design_pts], dp)
  s_UB <- seq(from = max_stress - (stress_distance*(design_pts-2)), to = max_stress, by = stress_distance)
  s <- numeric(design_pts)
  s[1] <- s1_d[1]
  for (i in 2:design_pts) 
  {
    s[i] <- min(s_UB[i-1], s[i-1] + s1_d[i])
  }
  
  zi <- x[(design_pts + 1):(2*design_pts)]
  p <- exp(zi)/sum(exp(zi))
  n_free <- test_units - (min_units*design_pts) 
  n <- round(n_free*p) + min_units
  diff <- test_units - sum(n)
  n[which.max(p)] <- n[which.max(p)] + diff
  
  optimal_values <- c(design_pts, s, n, rmse)
  names(optimal_values) <- c("N", paste("s", 1:design_pts, sep = ""), paste("n", 1:design_pts, sep = ""), "Min. RMSE")
  return(optimal_values)
}

### Function 4. Experiment configuration function ###
experiment_config_fnc <- function(n_sim = 1000, test_units = 100, test_duration = 1*24*365, min_units = 1, 
                                  stress_levels = 2, stress_distance = 10^(-5), max_stress = 0.9, min_stress = 0.1, design_stress = 0.05,   
                                  p_cens_limit = 0.8, filter_limit = 0.95,  
                                  pquant = 0.5, no_model_params = 3, dp = 5,
                                  DE_iter = 200, DE_step = 200, DE_tol = 0.05)
{
  config_check <- c(
  check1 = test_units > 0,
  check2 = all(c(design_stress, min_stress, max_stress, stress_distance) >= 0),
  check3 = test_duration > 0,
  check4 = min_stress < max_stress,
  check5 = dp >= 1,
  check6 = dp == round(dp),
  check7 = min_units >= 1,
  check8 = min_units == round(min_units),
  check9 = pquant >= 0 & pquant <= 1,
  check10 = no_model_params == 3,
  check11 = p_cens_limit >= 0 & p_cens_limit <= 1,
  check12 = filter_limit >= 0 & filter_limit <= 1,
  check13 = n_sim >= 1,
  check14 = n_sim == round(n_sim),
  check15 = stress_levels >= 2,
  check16 = all(c(max_stress, stress_distance) > 0),
  check17 = min_units <= test_units/stress_levels,
  check18 = stress_distance*stress_levels <= max_stress,
  check19 = all(c(DE_iter, DE_step, DE_tol) > 0))
  
  error_message <- c(
  check1 = "The number of test units must be greater than zero.",
  check2 = "The stress values cannot be negative.",
  check3 = "The test duration must be greater than zero.",
  check4 = "The minimum stress must be less than the maximum stress.",
  check5 = "The number of decimal places of accuracy must be at least one.",
  check6 = "The number of decimal places must be an integer.",
  check7 = "The number of test units at each stress level must be at least one.",
  check8 = "The minimum number of test units must be an integer.",
  check9 = "The specified quantile must be in the range [0,1].",
  check10 = "The number of initial model parameters must be three.",
  check11 = "The censoring proportion limit must be in the range [0,1].",
  check12 = "The outlier proportion limit must be in the range [0,1].",
  check13 = "The number of simulation replicates must be at least one.",
  check14 = "The number of simulation replicates must be an integer.",
  check15 = "The minimum allowable number of stress levels is two.",
  check16 = "The maximum stress and stress distance must be greater than zero.",
  check17 = "The desired number of test units at each stress level exceeds the total test units.",
  check18 = "(number of stress levels * distance between each stress) exceeds the maximum stress.",
  check19 = "DEoptim control/stopping criteria values must be greater than zero.")
  
  if(any(config_check == FALSE)) 
  {
    stop(paste("\n", error_message[config_check == FALSE], "\n"))
  }
}

### Function 5. DEoptim optimisation function ###
DEoptim_fnc <- function()
{
  experiment_config_check <- try(experiment_config_fnc(n_sim = n_sim, test_units = test_units, test_duration = test_duration, min_units = min_units, 
                                                       stress_levels = stress_levels, stress_distance = stress_distance, max_stress = max_stress, min_stress = min_stress, 
                                                       design_stress = design_stress, p_cens_limit = p_cens_limit, filter_limit = filter_limit, 
                                                       pquant = pquant, no_model_params = length(model_params),  dp = dp, 
                                                       DE_iter = DE_iter, DE_step = DE_step, DE_tol = DE_tol))                     
  
  if(inherits(experiment_config_check, "try-error"))
  {
    stop("Review and adjust experiment configuration settings before proceeding.")
  }else
    {
      #DEoptim parameter bounds
      s1LB <- min_stress
      s1UB <- max_stress - (stress_distance*(stress_levels - 1))
      dLB <- stress_distance
      dUB <- (s1UB + stress_distance) - s1LB
      DELB <- c(s1LB, rep(dLB, times = (stress_levels - 1)), rep(-5, times = stress_levels))
      DEUB <- c(s1UB, rep(dUB, times = (stress_levels - 1)), rep(5, times = stress_levels))
    
      DEoptim_results <- DEoptim(objective_fnc, lower = DELB, upper = DEUB, control = DEoptim.control(reltol = DE_tol, steptol = DE_step, itermax = DE_iter),
                               n_sim = n_sim, test_units = test_units, test_duration = test_duration, min_units = min_units,
                               stress_distance = stress_distance, max_stress = max_stress, design_stress = design_stress,  
                               p_cens_limit = p_cens_limit, filter_limit = filter_limit,  
                               pquant = pquant, model_params = model_params, dp = dp)
    
      return(DEoptim_results) 
    }
}

### Function 6. Print results function ###
print_results_fnc <- function()
{
  experiment_config_check <- try(experiment_config_fnc(n_sim = n_sim, test_units = test_units, test_duration = test_duration, min_units = min_units, 
                                                       stress_levels = stress_levels, stress_distance = stress_distance, max_stress = max_stress, min_stress = min_stress, 
                                                       design_stress = design_stress, p_cens_limit = p_cens_limit, filter_limit = filter_limit, 
                                                       pquant = pquant, no_model_params = length(model_params),  dp = dp, 
                                                       DE_iter = DE_iter, DE_step = DE_step, DE_tol = DE_tol)) 
  
  if(inherits(experiment_config_check, "try-error"))
  {
    stop("Review and adjust experiment configuration settings before proceeding.")
  }else
    {
      print_results <- results_fnc(x  = optimisation_results$optim$bestmem, test_units = test_units, min_units = min_units, 
                                 stress_distance = stress_distance, max_stress = max_stress, 
                                 dp = dp, rmse = optimisation_results$optim$bestval)
    
      return(print_results)
    }
}

##################################################
### User input - Test experiment configuration ###
##################################################

#Initial model parameters (b0, b1, gamma) (specified from preliminary experiment and/or SME input)
init_b0 <- 12.54 
init_b1 <- -19.48
init_gamma <- 1.98
model_params <- c(init_b0, init_b1, init_gamma) 

#Test experiment settings
test_units <- 100 #Number of components to be tested
test_duration <- 1*24*365 #Test experiment duration (hours)
min_units <- 1 #Minimum number of test components at each stress level (must be >= 1)
p_cens_limit <- 0.8 #Maximum allowable proportion of censoring in data set
filter_limit <- 0.95 #Minimum allowable proportion of observations remaining after removing outliers 
pquant <- 0.5 #Quantile at which to estimate design-stress lifetime
dp <- 5 #Stress level resolution, i.e., number of decimal places of accuracy (must >= 1)
n_sim <- 2 #Number of simulation replicates (fitted models) for a given test plan
stress_distance <- 10^(-5) #Minimum allowable distance between subsequent stress levels (must be > 0)
min_stress <- 0.1 #Minimum test stress
max_stress <- 0.9 #Maximum test stress
design_stress <- 0.05 #Nominal/operating stress
stress_levels <- 6 #Number of required test plan stress levels (design points) (must be >= 2)

#Optimisation stopping criteria 
DE_iter <- 2 #Number of DEoptim iterations (generations)
DE_step <- 2 #DEoptim step tolerance
DE_tol <- 0.05 #DEoptim relative convergence tolerance

#Note: User can alter the number of iterations and/or specify their own stopping criteria by choosing different 'DE_step' and 'DE_tol' values
#See ?DEoptim.control for further details

##################################
### Run test plan optimisation ###
##################################

optimisation_results <- DEoptim_fnc()
optimisation_results

###############################
### Print optimal test plan ###
###############################

optimal_test_plan <- print_results_fnc()
optimal_test_plan