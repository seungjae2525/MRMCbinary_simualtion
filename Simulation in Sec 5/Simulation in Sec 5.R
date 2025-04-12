################################################################################
################################################################################
rm(list=ls())

library(survival)
library(geepack)
library(lme4)
library(MASS)
library(dplyr)
library(glmmTMB)
library(R.utils)                 # For "withTimeout"

library(parallel)                # For "makeCluster"
library(doSNOW)                  # For "registerDoSNOW" and ".options.snow"
library(foreach)                 # For "%dopar%" and "foreach"
library(doParallel)              # For "registerDoParallel"
library(doRNG)                   # For "registerDoRNG"
library(progress)                # For "progress_bar"

################################################################################
## A set of parameters
B <- 1000

## Number of modality compared
n_K <- 3

## Combination of parameters
par.comb <- expand.grid(I = c(100, 500),
                        n_J = c(3, 5, 10),
                        beta_0 = c(-2.5, 0, 2.5),
                        beta_2k = c(-1.0, -0.5, 0.5, 1.0),
                        beta_3k = c(-1.0, -0.5, 0.5, 1.0))
par.comb <- par.comb[order(par.comb$beta_0, par.comb$beta_2k, par.comb$beta_3k, 
                           par.comb$I, par.comb$n_J), ]
rownames(par.comb) <- NULL

################################################################################
## Set a master seed number for reproducibility
set.seed(250103)

## Generate a unique set of seeds for each parameter combination in `par.comb`
seeds <- sample(1e8, nrow(par.comb))

## Setting for parallel processing
total_core <- parallel::detectCores()
cores <- c(1, 2, 5, 10, 20, 50, 100)[findInterval(total_core, c(0, 2, 5, 10, 20, 50, 100))]
cores

## Make progress bar
name.format <- paste0("(:spin) Progress: [:bar] :percent (:current/", nrow(par.comb),
                      ") [Elapsed time: :elapsedfull || Estimated time remaining: :eta]")
pb <- progress_bar$new(format = name.format,
                       total = nrow(par.comb),
                       complete = "=",     # Completion bar character
                       incomplete = "-",   # Incomplete bar character
                       current = ">",      # Current bar character
                       clear = FALSE,      # If TRUE, clears the bar when finish
                       width = 100)        # Width of the progress bar

################################################################################
### Make storage list
results.sim <- list()

## Start simulation
start <- Sys.time()
system.time(
  for(kkk in 1:nrow(par.comb)){
    if(kkk == 1){
      cat(paste0('[',Sys.time(),'] Start simulation! \n'))
    }
    
    ## Updates the current state
    pb$tick()
    
    ## split data by ourselves
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl, cores=cores)
    
    ## Get a unique seed for the current parameter combination
    seed <- seeds[kkk]
    ## Register a reproducible RNG stream for each iteration
    doRNG::registerDoRNG(seed)
    
    ## Run the simulation in parallel using foreach
    results_kkk <- foreach::foreach(iii = 1:B,
                                    .packages = c("parallel", "survival", "geepack",
                                                  "glmmTMB", "MASS", "dplyr", 
                                                  "R.utils","lme4"),
                                    .combine = 'c') %dopar%
      {
        ########################################################################
        ############################# Make dataset #############################
        n_I <- par.comb[kkk, 1]
        n_J <- par.comb[kkk, 2]
        
        ## Beta coefficients
        beta_0 <- par.comb[kkk, 3]
        beta_2k <- par.comb[kkk, 4]
        beta_3k <- par.comb[kkk, 5]
        
        ## Reader effects
        r <- rnorm(n = n_J, mean = 0, sd = 1)
        
        ## New beta0 and R_coef
        beta0_control <- r[1]
        r_coef <- r - beta0_control
        beta0_coef <- beta_0 - beta0_control
        
        ## Case effects
        v <- rnorm(n = n_I, mean = 0, sd = 1)
        
        ## Interaction effects
        u <- matrix(rnorm(n = n_I * n_J, mean = 0, sd = 1), nrow = n_I, ncol = n_J)
        
        ## generate index
        i_vals <- rep(x = 1:n_I, each = n_J * n_K)
        j_vals <- rep(x = rep(x = 1:n_J, each = n_K), times = n_I)
        k_vals <- rep(x = 1:n_K, times = n_I * n_J)
        
        ## Z Indicator
        z_ij2 <- ifelse(k_vals == 2, 1, 0)
        z_ij3 <- ifelse(k_vals == 3, 1, 0)
        
        ## binary outcome
        pp <- plogis(beta_0 + beta_2k * z_ij2 + beta_3k * z_ij3 + 
                       r[j_vals] + v[i_vals] + u[cbind(i_vals, j_vals)])
        y_ijk <- rbinom(n = length(pp), size = 1, prob = pp)
        
        ## Generate data
        dat <- data.frame(ID = i_vals, Reader = j_vals, Modality = k_vals, Y = y_ijk)
        dat$ID <- factor(dat$ID)
        dat$Reader <- factor(dat$Reader)
        dat$Modality <- factor(dat$Modality)
        
        ########################################################################
        rm(n_I, n_J, r, beta0_control, beta0_coef, v, u, 
           i_vals, j_vals, k_vals, z_ij2, z_ij3, pp, y_ijk)
        
        ########################################################################
        ### Fitting Clogit
        ## CLR
        fit_CLR_message <- "No warnings or errors."
        fit_CLR <- tryCatch({
          R.utils::withTimeout({
            withCallingHandlers({
              model <- survival::clogit(formula = Y ~ Modality + strata(ID, Reader), 
                                        data = dat,
                                        control = coxph.control(iter.max = 1e3))
              model
            }, warning = function(w) {
              fit_CLR_message <<- paste0("WWWW ", conditionMessage(w))
              invokeRestart("muffleWarning")
            })
          }, timeout = 120)  # Set timeout to 120 seconds
        }, error = function(e) {
          fit_CLR_message <<- paste0("EEEE ", conditionMessage(e))
          NA
        })
        if (is.na(fit_CLR_message)) {
          fit_CLR_message <- "No warnings or errors."
        }; rm(model)
        if (class(fit_CLR)[1] == "clogit") {
          fit_CLR <- summary(fit_CLR)$coefficients
        } else {
          fit_CLR <- NA
        }
        
        
        ### Fitting GLM
        ## GLM
        fit_GLM_message <- "No warnings or errors."
        fit_GLM <- tryCatch({
          R.utils::withTimeout({
            withCallingHandlers({
              model <- stats::glm(formula = Y ~ Modality, 
                                  family = binomial(link = "logit"), 
                                  data = dat, control = glm.control(maxit = 1e3))
              model
            }, warning = function(w) {
              fit_GLM_message <<- paste0("WWWW ", conditionMessage(w))
              invokeRestart("muffleWarning")
            })
          }, timeout = 120)  # Set timeout to 120 seconds
        }, error = function(e) {
          fit_GLM_message <<- paste0("EEEE ", conditionMessage(e))
          NA
        })
        if (is.na(fit_GLM_message)) {
          fit_GLM_message <- "No warnings or errors."
        }; rm(model)
        if (class(fit_GLM)[1] == "glm") {
          fit_GLM <- summary(fit_GLM)$coefficients
        } else {
          fit_GLM <- NA
        }
        
        
        ### Fitting GEE
        ## GEE1
        fit_GEE1_message <- "No warnings or errors."
        fit_GEE1 <- tryCatch({
          R.utils::withTimeout({
            withCallingHandlers({
              model <- geepack::geeglm(formula = Y ~ Modality, id = ID,
                                       family = binomial(link = "logit"), 
                                       corstr = "independence", data = dat, 
                                       control = geese.control(maxit = 1e3))
              model
            }, warning = function(w) {
              fit_GEE1_message <<- paste0("WWWW ", conditionMessage(w))
              invokeRestart("muffleWarning")
            })
          }, timeout = 120)  # Set timeout to 120 seconds
        }, error = function(e) {
          fit_GEE1_message <<- paste0("EEEE ", conditionMessage(e))
          NA
        })
        if (is.na(fit_GEE1_message)) {
          fit_GEE1_message <- "No warnings or errors."
        }; rm(model)
        fit_GEE1_e_which <- sum(diag(vcov(fit_GEE1)) < 0) != 0
        if (class(fit_GEE1)[1] == "geeglm") {
          fit_GEE1 <- summary(fit_GEE1)$coefficients
        } else {
          fit_GEE1 <- NA
        }
        
        ## GEE2
        fit_GEE2_message <- "No warnings or errors."
        fit_GEE2 <- tryCatch({
          R.utils::withTimeout({
            withCallingHandlers({
              model <- geepack::geeglm(formula = Y ~ Modality, id = interaction(ID, Reader),
                                       family = binomial(link = "logit"), 
                                       corstr = "independence", data = dat, 
                                       control = geese.control(maxit = 1e3))
              model
            }, warning = function(w) {
              fit_GEE2_message <<- paste0("WWWW ", conditionMessage(w))
              invokeRestart("muffleWarning")
            })
          }, timeout = 120)  # Set timeout to 120 seconds
        }, error = function(e) {
          fit_GEE2_message <<- paste0("EEEE ", conditionMessage(e))
          NA
        })
        if (is.na(fit_GEE2_message)) {
          fit_GEE2_message <- "No warnings or errors."
        }; rm(model)
        fit_GEE2_e_which <- sum(diag(vcov(fit_GEE2)) < 0) != 0
        if (class(fit_GEE2)[1] == "geeglm") {
          fit_GEE2 <- summary(fit_GEE2)$coefficients
        } else {
          fit_GEE2 <- NA
        }
        
        
        ### Fitting GLMM using glmmTMB
        ## glmmTMB1
        fit_glmmTMB1_message <- "No warnings or errors."
        fit_glmmTMB1 <- tryCatch({
          R.utils::withTimeout({
            withCallingHandlers({
              model <- glmmTMB::glmmTMB(formula = Y ~ Modality + (1 | ID),
                                        family = binomial(link = "logit"), 
                                        data = dat, REML = FALSE,
                                        control = glmmTMBControl(optCtrl = list(iter.max = 1e3, 
                                                                                eval.max = 1e3)))
              model
            }, warning = function(w) {
              fit_glmmTMB1_message <<- paste0("WWWW ", conditionMessage(w))
              invokeRestart("muffleWarning")
            })
          }, timeout = 120)  # Set timeout to 120 seconds
        }, error = function(e) {
          fit_glmmTMB1_message <<- paste0("EEEE ", conditionMessage(e))
          NA
        })
        if (is.na(fit_glmmTMB1_message)) {
          fit_glmmTMB1_message <- "No warnings or errors."
        }; rm(model)
        if (class(fit_glmmTMB1)[1] == "glmmTMB") {
          fit_glmmTMB1 <- summary(fit_glmmTMB1)$coefficients$cond
        } else {
          fit_glmmTMB1 <- NA
        }
        
        ## glmmTMB2
        fit_glmmTMB2_message <- "No warnings or errors."
        fit_glmmTMB2 <- tryCatch({
          R.utils::withTimeout({
            withCallingHandlers({
              model <- glmmTMB::glmmTMB(formula = Y ~ Modality + (1 | ID) + (1 | Reader) + (1 | ID:Reader),
                                        family = binomial(link = "logit"), 
                                        data = dat, REML = FALSE,
                                        control = glmmTMBControl(optCtrl = list(iter.max = 1e3, 
                                                                                eval.max = 1e3)))
              model
            }, warning = function(w) {
              fit_glmmTMB2_message <<- paste0("WWWW ", conditionMessage(w))
              invokeRestart("muffleWarning")
            })
          }, timeout = 120)  # Set timeout to 120 seconds
        }, error = function(e) {
          fit_glmmTMB2_message <<- paste0("EEEE ", conditionMessage(e))
          NA
        })
        if (is.na(fit_glmmTMB2_message)) {
          fit_glmmTMB2_message <- "No warnings or errors."
        }; rm(model)
        if (class(fit_glmmTMB2)[1] == "glmmTMB") {
          fit_glmmTMB2 <- summary(fit_glmmTMB2)$coefficients$cond
        } else {
          fit_glmmTMB2 <- NA
        }
        
        
        ### Fitting GLMM using glmer
        ## glmer1
        fit_glmer1_message <- "No warnings or errors."
        fit_glmer1 <- tryCatch({
          R.utils::withTimeout({
            withCallingHandlers({
              model <- lme4::glmer(formula = Y ~ Modality + (1 | ID),
                                   family = binomial(link = "logit"), data = dat)
              model
            }, warning = function(w) {
              fit_glmer1_message <<- paste0("WWWW ", conditionMessage(w))
              invokeRestart("muffleWarning")
            })
          }, timeout = 120)  # Set timeout to 120 seconds
        }, error = function(e) {
          fit_glmer1_message <<- paste0("EEEE ", conditionMessage(e))
          NA
        })
        if (is.na(fit_glmer1_message)) {
          fit_glmer1_message <- "No warnings or errors."
        }; rm(model)
        fit_glmer1_w_which <- ifelse(class(fit_glmer1)[1] == "logical", TRUE, isSingular(fit_glmer1))
        if (class(fit_glmer1)[1] == "glmerMod") {
          fit_glmer1 <- summary(fit_glmer1)$coefficients
        } else {
          fit_glmer1 <- NA
        }
        
        ## glmer2
        fit_glmer2_message <- "No warnings or errors."
        fit_glmer2 <- tryCatch({
          R.utils::withTimeout({
            withCallingHandlers({
              model <- lme4::glmer(formula = Y ~ Modality + (1 | ID) + (1 | Reader) + (1 | ID:Reader),
                                   family = binomial(link = "logit"), data = dat)
              model
            }, warning = function(w) {
              fit_glmer2_message <<- paste0("WWWW ", conditionMessage(w))
              invokeRestart("muffleWarning")
            })
          }, timeout = 120)  # Set timeout to 120 seconds
        }, error = function(e) {
          fit_glmer2_message <<- paste0("EEEE ", conditionMessage(e))
          NA
        })
        if (is.na(fit_glmer2_message)) {
          fit_glmer2_message <- "No warnings or errors."
        }; rm(model)
        fit_glmer2_w_which <- ifelse(class(fit_glmer2)[1] == "logical", TRUE, isSingular(fit_glmer2))
        if (class(fit_glmer2)[1] == "glmerMod") {
          fit_glmer2 <- summary(fit_glmer2)$coefficients
        } else {
          fit_glmer2 <- NA
        }
        
        ########################################################################
        list(dplyr::lst(iii, r_coef,
                        fit_CLR, fit_CLR_message,
                        
                        fit_GLM, fit_GLM_message, 
                        
                        fit_GEE1, fit_GEE1_message, fit_GEE1_e_which,
                        fit_GEE2, fit_GEE2_message, fit_GEE2_e_which,
                        
                        fit_glmmTMB1, fit_glmmTMB1_message,
                        fit_glmmTMB2, fit_glmmTMB2_message,
                        
                        fit_glmer1, fit_glmer1_message, fit_glmer1_w_which,
                        fit_glmer2, fit_glmer2_message, fit_glmer2_w_which))
      }
    ############################################################################
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    
    ## Save global results
    results.sim[[kkk]] <- results_kkk
    rm(results_kkk)
  })
(end <- Sys.time() - start)
cat(paste0('[',Sys.time(),'] End simulation! \n'))

################################################################################
# setwd("~/")
save.image(file="Simulation in Sec 5.RData")

