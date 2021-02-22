# Pure measurement error functions for simulation:
#   getMSE
#   me_pure_data
#   me_pure

# Load the broom package for the tidy() function
pacman::p_load(broom, stringr)


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

# make dataframe to simulate pure measurement error
me_pure_data <- function(n_subj = 10000) {
    # definition of terms:
    #   n for sample size (n_subj), referring to subject 
    #   s for exposure model covariates
    #   y for outcome
    # the following terms are passed from outside if not using defaults:
    #   n_subj=10000 (default) is the # subjects
    # the following terms are set inside the program: (could edit the function to
    # pass these in from outside)
    #   sd_eta=4  is the SD for the error in the exposure model
    #   sd_e1=4   is the SD for the extra classical error of exposure
    #   sd_e2=4   is the SD for the extra classical error of exposure
    #   sd_e3=4   is the SD for the extra classical error of exposure
    #   sd_eps=25 is the SD for the error in the disease model
    #   alpha_0=0 is the intercept parameter in the exposure model
    #   alpha[1]=4is the parameter for s1 in the exposure model
    #   alpha[2]=4is the parameter for s2 in the exposure model
    #   alpha[3]=4is the parameter for s3 in the exposure model
    #   beta[1]=1   is the intercept in the disease model
    #   beta[2]=2   is the slope in the disease model (called \beta_x in the lab)
    # define the exposure and health model parameters for fixed effects
    # define the coefficients alpha and beta (description given above)
    alpha_0 <- 0
    alpha <- c(4, 4, 4)
    beta  <- c(1, 2)
    
    # define the SDs for all the components for the subjects and samples 
    sd_s <- c(1, 1, 1)
    sd_eta <- 4
    sd_e <- c(4, 4, 4)
    sd_eps <- 25
    
    
    # create the subject dataset, using n_subj as supplied in the
    # function's parameter list:
    tibble( 
        s_1 = rnorm(n_subj, sd = sd_s[1]),
        s_2 = rnorm(n_subj, sd = sd_s[2]),
        s_3 = rnorm(n_subj, sd = sd_s[3]),
        
        x = alpha_0 + alpha[1] * s_1 + alpha[2] * s_2 + alpha[3] * s_3 +
            rnorm(n_subj, sd = sd_eta),
        y = beta[1] + beta[2] * x + rnorm(n_subj, sd = sd_eps),
        
        Berk_1 = alpha_0 + alpha[1] * s_1,
        Berk_2 = alpha_0 + alpha[1] * s_1 + alpha[2] * s_2,
        Berk_3 = alpha_0 + alpha[1] * s_1 + alpha[2] * s_2 + alpha[3] * s_3,
        
        class_1 = x       + rnorm(n_subj, sd = sd_e[1]),
        class_2 = class_1 + rnorm(n_subj, sd = sd_e[2]),
        class_3 = class_2 + rnorm(n_subj, sd = sd_e[3])
    )
    
}


# me_pure function
me_pure <- function(n_subj = 10000){
    
    # create dataframe with simulated measurement error
    d <- me_pure_data(n_subj)
    
    # list predictors in d, looking for x, berkson and classical variable names
    predictors <- str_subset(names(d), "x|Berk_|class_")
    
    # Define a tibble of 4 parameters x 7 models
    ret <- lapply(predictors, function(i) {
        
        # specify formula
        frmla <- as.formula(paste("y ~", i))
        
        # fit linear model
        lmfit <- lm(frmla, data = d)
        
        #compile parameters of interest
        tibble(b1 = tidy(lmfit)$estimate[2],
               seb1 = tidy(lmfit)$std.error[2],
               R2_W_reg = cor(d[[i]],d$x)^2, 
               R2_W_MSE = as.numeric(getMSE(d$x, d[[i]])[2]),
               exp_var = var(d[[i]])
        )
    }) %>% 
        
        # Set names for list items
        setNames(predictors) %>% 
        
        # bind list rows together
        bind_rows(.id = "predictor")
    
    # Return the dataframe
    return(ret)
}

