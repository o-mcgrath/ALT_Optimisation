############################################################################
### Test plan optimisation for Weibull quadratic AFT model (Specified N) ###
############################################################################

#The user specifies initial model parameters (following a preliminary experiment
#and/or consultation with subject matter expert(s)).
#The user specifies the required number of stresses (N).
#The optimisation algorithm returns the optimal test plan, i.e.,
#the optimal stress levels and the associated number of test units.

###################################################
### Load DEoptim package and required functions ###
###################################################

### Install DEoptim package and load DEoptim library ###
#install.packages("DEoptim") 
library(DEoptim) 

### Function 1. Weibull log-likelihood function (quadratic AFT model) ###

#Function inputs are the model parameters and (simulated) lifetime data.
#Function output is the negative log likelihood.

loglikfnc <- function(theta, data) 
{
  b0 <- theta[1]
  b1 <- theta[2]
  b2 <- theta[3]
  a0 <- theta[4] 
  S <- data$stress
  n <- nrow(data)
  ti <- data$lifetime
  di <- data$cens
  mu <- b0 + b1*S + b2*(S^2)
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
  dl_db2 <- sum(-di*gam*(S^2) + gam*(S^2)*(ti/lam)^gam)
  dl_da0 <- sum(di + (log(ti) - log(lam))*di*gam - log(ti/lam)*gam*(ti/lam)^gam)
  grad_vec <- -c(dl_db0, dl_db1, dl_db2, dl_da0) 
  
  if(any(!is.finite(grad_vec))) 
  {
    return(10^10) 
  }
  
  attr(loglik, "gradient") <- grad_vec
  return(-1*loglik) 
}

### Function 2. Test plan objective function ###

#This is the optimisation objective function. The function inputs are the test
#plan parameters. For a given test plan configuration, an AFT model is fitted,
#and the design-stress lifetime (at a specified quantile) is estimated. 
#This is repeated 'n_sim' times and the function returns the RMSE of estimates.

objective_function <- function(x, model_params = c(13, -38, 18, 2), test_units = 100,
                               test_duration = 1*24*365, min_units = 1, 
                               p_cens_limit = 0.8, filter_limit = 0.95,
                               pquant = 0.5, dp = 5, stress_distance = 10^(-5),
                               max_stress = 0.9, design_stress = 0.05,  
                               n_sim = 1000) 
{
  DEparams <- as.integer(length(x))
  design_pts <- as.integer(DEparams/2) 
  s1_d <- round(x[1:design_pts], dp) 
  s_UB <- seq(from = max_stress - (stress_distance*(design_pts-2)), 
              to = max_stress, by = stress_distance)
  
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
  B2 <- model_params[3]
  gamma <- model_params[4]
  stress <- rep(s, times = n)
  lambda <- exp(B0 + B1*stress + B2*(stress^2))
  true_quantile <- exp(B0 + B1*design_stress + B2*(design_stress^2))*
    (-log(1 - pquant))^(1/gamma)
  
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
      result <- try(nlm(loglikfnc, p = init_values[[j]], data = dataset, 
                        iterlim = 1000, check.analyticals = FALSE), silent = TRUE)
      
      if (!inherits(result, "try-error") && is.finite(result$minimum)) 
      {
        mle_b0 <- result$estimate[1]
        mle_b1 <- result$estimate[2]
        mle_b2 <- result$estimate[3]
        mle_gamma <- exp(result$estimate[4])
        if (is.finite(mle_gamma)) 
        {
          est_quantile[i] <- exp(mle_b0 + mle_b1*design_stress + 
                      mle_b2*(design_stress^2))*(-log(1 - pquant))^(1/mle_gamma)
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
  est_quantile <- est_quantile[est_quantile >= Q1 - 6*IQR 
                               & est_quantile <= Q3 + 6*IQR]
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

### Function 3. ALT optimisation function ###

#This is the ALT test plan optimisation function. The experiment configuration is 
#checked before the DE optimisation algorithm (DEoptim) is called. 
#For the specified test plan configuration, DEoptim minimises the objective 
#function and returns the optimisation output and a printout of the optimal test 
#plan settings. The printout includes the number of stresses (N), the stress 
#levels (s1 s2 ... sj), the test units (n1 n2 ... nj) and the RMSE.

ALT_optimise <- function(model_params = c(13, -38, 18, 2), test_units = 100,  
                         test_duration = 1*24*365, min_units = 1,
                         p_cens_limit = 0.8, filter_limit = 0.95,
                         pquant = 0.5, dp = 5, stress_distance = 10^(-5),
                         min_stress = 0.1, max_stress = 0.9, design_stress = 0.05,
                         stress_levels = 3, DE_tol = 0.05, n_sim = 1000, 
                         DE_iter = 200, DE_step = 200)
{
  #Check test experiment configuration
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
  check10 = length(model_params) == 4,
  check11 = p_cens_limit >= 0 & p_cens_limit <= 1,
  check12 = filter_limit >= 0 & filter_limit <= 1,
  check13 = n_sim >= 1,
  check14 = n_sim == round(n_sim),
  check15 = stress_levels >= 3,
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
  check10 = "The number of initial model parameters must be four.",
  check11 = "The censoring proportion limit must be in the range [0,1].",
  check12 = "The outlier proportion limit must be in the range [0,1].",
  check13 = "The number of simulation replicates must be at least one.",
  check14 = "The number of simulation replicates must be an integer.",
  check15 = "The minimum allowable number of stress levels is three.",
  check16 = "The maximum stress and stress distance must be greater than zero.",
  check17 = "The desired number of test units at each stress level exceeds the total test units.",
  check18 = "(number of stress levels * distance between each stress) exceeds the maximum stress.",
  check19 = "DEoptim control/stopping criteria values must be greater than zero.")
  
  if(any(config_check == FALSE)) 
  {
    stop(paste("\n", error_message[config_check == FALSE], "\n"))
  }
  
  #Run DE optimisation algorithm (via DEoptim)
  s1LB <- min_stress
  s1UB <- max_stress - (stress_distance*(stress_levels - 1))
  dLB <- stress_distance
  dUB <- (s1UB + stress_distance) - s1LB
  DELB <- c(s1LB, rep(dLB, times = (stress_levels - 1)), 
            rep(-5, times = stress_levels))
  DEUB <- c(s1UB, rep(dUB, times = (stress_levels - 1)), 
            rep(5, times = stress_levels))
  
  DEoptim_results <- DEoptim(objective_function, lower = DELB, upper = DEUB, 
                             control = DEoptim.control(reltol = DE_tol, 
                             steptol = DE_step, itermax = DE_iter),
                             model_params = model_params, test_units = test_units, 
                             test_duration = test_duration, min_units = min_units,
                             p_cens_limit = p_cens_limit, filter_limit = filter_limit, 
                             pquant = pquant, dp = dp, stress_distance = stress_distance, 
                             max_stress = max_stress, design_stress = design_stress,  
                             n_sim = n_sim)
  
  #Print optimal test plan
  x <- DEoptim_results$optim$bestmem
  optim_rmse <- DEoptim_results$optim$bestval
  DEparams <- as.integer(length(x))
  design_pts <- as.integer(DEparams/2)
  s1_d <- round(x[1:design_pts], dp)
  s_UB <- seq(from = max_stress - (stress_distance*(design_pts-2)), 
              to = max_stress, by = stress_distance)
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
  
  optimal_test_plan <- c(design_pts, s, n, optim_rmse)
  names(optimal_test_plan) <- c("N", paste("s", 1:design_pts, sep = ""), 
                                paste("n", 1:design_pts, sep = ""), "Min. RMSE")
  
  optimisation_results <- list(DEoptim_results, optimal_test_plan)
  names(optimisation_results) <- c("DEoptim output", "Optimal test plan")
  
  #Print details of the optimal test plan
  return(optimisation_results) 
  
}

##################################################
### User input - Test experiment configuration ###
##################################################

#The user specifies the test experiment settings in this section.
#User inputs:
#1. Initial/preliminary model parameters.
#2. Specific test experiment configuration/constraints.
#3. Optimisation stopping criteria.

### 1. Initial model parameters (b0, b1, b2, gamma) ###

#Specified from preliminary experiment and/or subject matter expert(s).
init_b0 <- 13.41 
init_b1 <- -37.89
init_b2 <- 17.72
init_gamma <- 2.07
model_params <- c(init_b0, init_b1, init_b2, init_gamma) 

### 2. Test experiment configuration ###

#Number of components to be tested.
test_units <- 100 
#Test experiment duration (hours).
test_duration <- 1*24*365 
#Minimum number of test components at each stress level (must be >= 1).
min_units <- 1 
#Maximum allowable proportion of censoring in data set.
p_cens_limit <- 0.8
#Minimum allowable proportion of observations remaining after removing outliers.
filter_limit <- 0.95
#Quantile at which to estimate design-stress lifetime.
pquant <- 0.5 
#Stress level resolution, i.e., number of decimal places of accuracy (must >= 1).
dp <- 5 
#Minimum allowable distance between subsequent stress levels (must be > 0).
stress_distance <- 10^(-5) 
#Minimum test stress.
min_stress <- 0.1
#Maximum test stress.
max_stress <- 0.9 
#Nominal/operating stress.
design_stress <- 0.05
#Number of simulation replicates (fitted models) for a given test plan.
n_sim <- 1000 
#Number of required test plan stress levels (design points) (must be >= 3).
stress_levels <- 6

### 3. Optimisation stopping criteria ###

DE_iter <- 200 #Number of DEoptim iterations (generations).
DE_step <- 200 #DEoptim step tolerance.
DE_tol <- 0.05 #DEoptim relative convergence tolerance.

#User can specify their own stopping criteria by adjusting the above values.
#See ?DEoptim.control for further details.

##############################################
### Test plan optimisation (Demonstration) ###
##############################################

#**The following simulation is for demonstration purposes only**
#The settings are chosen to demonstrate the optimisation process/output.
#**These settings should not be used in practice.**
#Execution time ~ minute

set.seed(202508)
demo_optimisation_results <- ALT_optimise(model_params = model_params, 
                                          test_units = test_units, 
                                          test_duration = test_duration, 
                                          min_units = min_units,
                                          p_cens_limit = p_cens_limit, 
                                          filter_limit = filter_limit,
                                          pquant = pquant, dp = dp, 
                                          stress_distance = stress_distance, 
                                          min_stress = min_stress, 
                                          max_stress = max_stress,  
                                          design_stress = design_stress,   
                                          stress_levels = stress_levels,
                                          DE_tol = DE_tol,
                                          n_sim = 5, DE_iter = 5, DE_step = 5)
demo_optimisation_results

#####################################################
### Test plan optimisation (exploratory analysis) ###
#####################################################

#The following optimisation settings may be suitable for exploratory analysis.
#Execution time ~ hour

set.seed(202508)
exploratory_optimisation_results <- ALT_optimise(model_params = model_params, 
                                                 test_units = test_units, 
                                                 test_duration = test_duration, 
                                                 min_units = min_units,
                                                 p_cens_limit = p_cens_limit, 
                                                 filter_limit = filter_limit,
                                                 pquant = pquant, dp = dp, 
                                                 stress_distance = stress_distance, 
                                                 min_stress = min_stress, 
                                                 max_stress = max_stress,  
                                                 design_stress = design_stress,   
                                                 stress_levels = stress_levels,
                                                 DE_tol = DE_tol,
                                                 n_sim = 250, DE_iter = 25, 
                                                 DE_step = 25)
exploratory_optimisation_results

##################################################
### Test plan optimisation (detailed analysis) ###
##################################################

#The following optimisation settings may be suitable for detailed/final analysis.
#Values of n_sim = 1000, DE_iter =  200, and DE_step = 200 were used in the paper.
#Execution time ~ day

set.seed(202508)
final_optimisation_results <- ALT_optimise(model_params = model_params, 
                                           test_units = test_units, 
                                           test_duration = test_duration, 
                                           min_units = min_units,
                                           p_cens_limit = p_cens_limit, 
                                           filter_limit = filter_limit,
                                           pquant = pquant, dp = dp, 
                                           stress_distance = stress_distance, 
                                           min_stress = min_stress, 
                                           max_stress = max_stress,  
                                           design_stress = design_stress,   
                                           stress_levels = stress_levels,
                                           DE_tol = DE_tol,
                                           n_sim = 1000, DE_iter = 200, 
                                           DE_step = 200)
final_optimisation_results

