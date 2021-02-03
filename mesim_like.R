# Pure measurement error functions for simulation:
#   getMSE
#   me_like

# Load the broom package for the tidy() function
# Load tidyverse for tibble()
# (only needed if the calling code doesn't load these)
pacman::p_load(broom, tidyverse)
#
#
# getMSE function:
getMSE <- function(obs,pred) {
    # obs is the outcome variable
    # pred is the prediction from a model
    obs_avg <- mean(obs)
    MSE_obs <- mean((obs - obs_avg) ^ 2)
    MSE_pred <- mean((obs - pred) ^ 2)
    result <- c(sqrt(MSE_pred),
                max(1 - MSE_pred / MSE_obs, 0))
    names(result) <-  c("RMSE", "R2(MSE-based)")
    
    return(result)
}

create_dataset <- function(n, sd_s, sd_eta, alpha_0, alpha) {
    tibble(
        s_1 = rnorm(n, sd = sd_s[1]),
        s_2 = rnorm(n, sd = sd_s[2]),
        s_3 = rnorm(n, sd = sd_s[3]),
        x = alpha_0 + alpha[1] * s_1 + alpha[2] * s_2 + alpha[3] * s_3 +
            rnorm(n, sd = sd_eta)
    )
}

# me_like function
me_like <- function(n_subj = 10000, n_samp = 100, s3_sd1 = 1, s3_sd2 = 0.3) {
# definition of terms:
#   n for sample size (n_subj and n_samp), referring to subject and
#       samples/monitors, respectively
#   s for exposure model covariates
#   y for outcome
# the following terms are passed from outside if not using defaults:
#   n_subj=10,000 (default) is the no. subjects
#   n_samp=100 (default) is the no. monitors/monitor locations/exposure samples
#   s3_sd is the SD for the s3 variable in the sample 
#       and allowed to be 1 or 0.3  (All other SDs for the s's are 1)
#       s3_sd1=1   (default)
#       s3_sd2=0.3 (default)
# the following terms are set inside the program:
#   sd_eta=4 is the SD for the error in the exposure model
#   sd_eps=25 is the SD for the error in the disease model
#   alpha1=4  is the parameter for s1  in the exposure model
#   alpha2=4  is the parameter for s2  in the exposure model
#   alpha3=4  is the parameter for s3  in the exposure model
#   beta0=1   is the intercept in the disease model
#   beta1=2   is the slope in the disease model
# define the exposure and health model parameters for fixed effects
    alpha_0 <- 0
    alpha <- c(4, 4, 4)
    beta  <- c(1, 2)
    
# define the SDs for all the components for the subjects and samples
# Note:  for sd1, the true and sample1 dataset will have the same SDs
    sd1_s <- c(1, 1, s3_sd1)
    sd2_s <- c(1, 1, s3_sd2)
    sd_eta <- 4
    sd_e <- c(4, 4, 4)
    sd_eps <- 25
    
# first create the subject dataset, using n_subj as supplied in the
#   function's parameter list:
    
    subj <- create_dataset(n_subj, sd1_s, sd_eta, alpha_0, alpha)
    subj <- subj %>% 
        mutate(y = beta[1] + beta[2] * x + rnorm(n_subj, sd = sd_eps))

# now create the 1st monitoring dataset drawn from the same distribution as the
# subject data

    samp1 <- create_dataset(n_samp, sd1_s, sd_eta, alpha_0, alpha)
    
# now create the 2nd monitoring dataset drawn from a different distribution not
# the same as the subject data distribution
    
    samp2 <- create_dataset(n_samp, sd2_s, sd_eta, alpha_0, alpha)
    
 # Now develop exposure predictions using the first monitor datset and predict
 #  on the subjects
    # Fully specified exposure model
    lm_1full <- lm(x ~ s_1 + s_2 + s_3, data = samp1)
    #xhat_1full <- predict(lm_1full, data = samp1)
    a3hat_1full <- tidy(lm_1full)$estimate[4]
    a3var_1full <- tidy(lm_1full)$std.error[4]^2
    r2_1full <- glance(lm_1full)$r.squared
    
    subj <- subj %>%
            modelr::add_predictions(lm_1full,"xhat_1full")
    
    # Reduced exposure model
    lm_1red <- lm(x ~ s_1 + s_2, data = samp1)
    r2_1red <- glance(lm_1red)$r.squared
    
    subj <- subj %>%
        modelr::add_predictions(lm_1red,"xhat_1red")
    
# Now develop exposure predictions using the second monitor datset and predict
# on the subjects
    # Fully specified exposure model
    lm_2full <- lm(x ~ s_1 + s_2 + s_3, data = samp2)
    a3hat_2full <- tidy(lm_2full)$estimate[4]
    a3var_2full <- tidy(lm_2full)$std.error[4]^2
    r2_2full <- glance(lm_2full)$r.squared
    
    subj <- subj %>%
        modelr::add_predictions(lm_2full,"xhat_2full")
    
    # Reduced exposure model
    lm_2red <- lm(x ~ s_1 + s_2, data = samp2)
    r2_2red <- glance(lm_2red)$r.squared
    
    subj <- subj %>%
        modelr::add_predictions(lm_2red,"xhat_2red")
    
# The next step is commented out because it is subsumed in the following 
#   lapply() step
# Do health analyses:
    # True exposure
#    lm_xtrue <- lm(y ~ x, data = subj)
    # Predicted expsure, monitor set 1 fully specified exposure model
#    lm_x1full <- lm(y ~ xhat_1full, data = subj)
    # Predicted expsure, monitor set 1 reduced exposure model
#    lm_x1red <- lm(y ~ xhat_1red, data = subj)
    # Predicted expsure, monitor set 2 fully specified exposure model
#    lm_x2full <- lm(y ~ xhat_2full, data = subj)
    # Predicted expsure, monitor set 2 reduced exposure model
#    lm_x2red <- lm(y ~ xhat_2red, data = subj)
     
 # collect b1, Vb1, R2 for disease model, variance of xhat, 
    predictors <- c("x",
                    "xhat_1full",
                    "xhat_1red",
                    "xhat_2full",
                    "xhat_2red")
    
    # Define a list for 4 parameters from disease model fits x 5 exposures
    ret.list <- lapply(predictors, function(i) {
        lmfit <- lm(subj$y ~ subj[[i]])
        list(tidy(lmfit)$estimate[2], 
             tidy(lmfit)$std.error[2],
             as.numeric(getMSE(subj$y, lmfit$fitted.values)[2]),
             var(subj[[i]]))
    })
    
    # Set names for list items and the items contained in each item
    names(ret.list) <- predictors
    for (i in 1:length(ret.list)) {
        names(ret.list[[i]]) <- c('b1', 'seb1', 'R2_y', 'exp_var')
    }
    
    # Add additional parameters (a3hat, a3var, r2) for (1, 2) and (full, red)
    ret.list$x[c('a3hat', 'a3var', 'r2')] <- NA
    for (i in 1:2) {
        for (j in c('full', 'red')) {
            for (k in c('a3hat', 'a3var', 'r2')) {
                varname <- paste0(k, '_', i, j)
                ret.list[[paste0('xhat_', i, j)]][[k]] <- 
                    ifelse(exists(varname), get(varname), NA)
            }
        }
    }
    
    # Return the list
    return(ret.list)
}


