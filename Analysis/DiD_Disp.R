############################################
#
# This file runs the DiDs
#
#############################################


# Preamble ----------------------------------------------------------------

# Clearing the workspace
rm(list = ls())

# Loading the packages that will be needed. 
require(caret)
require(labelled)
require(vtreat)
require(doParallel)
require(tidyverse)
require(gtools)
library(MatrixModels)
library(matrixStats)
library(Matrix)
library(lfe)
library(magrittr)


# Load in functions and options file.
source(file = "./Code/Functions_and_Options.R")

# Setting a seed 
set.seed(42)

# Number of bootstraps
nboot <- 100


# Load Data ---------------------------------------------------------------
load(file = "./Data/ucr_data_stack_cleaned_subset.RData")

# Subset to local police and sheriffs
ucr_dat_stack <- 
  ucr_dat_stack %>% 
  filter(agency_type %in% c("local police department", "sheriffs office"))

# Drop early legalizations
ucr_dat_stack <- 
  ucr_dat_stack %>% 
  filter(disp_eligible == 1 & 
           year < 2020)



# DiD Function ----------------------------------------------------------
sparse_impute_did <- function(x, y, stage2, data, treat_var, weights = NULL){
  
  treat_bin <- ceiling(treat_var) # Get 0,1 treatment variable
  
  l0 <- which(treat_bin == 0) # get subset that is untreated
  
  # Fit OLS model on control units
  mod1 <- MatrixModels:::lm.fit.sparse(y = y[l0],
                                       x = x[l0,],
                                       w = weights[l0])
  
  # Get estimated treatment effect
  treat_effect <- (y - x%*%mod1)
  data$treat_effect <- treat_effect[,1]
  
  # Run model of treatment effects on
  mod1 <- felm(stage2,
               weights = weights,
               data = data,
               exactDOF = TRUE)
  
  return(mod1)
  
  
}


# Triple Diffs -------------------------------------------------------------
sparse_dat_stack <- sparse.model.matrix(as.formula(paste("", "~ type_ori + type_year + year_ori")), data = ucr_dat_stack)




# Bootstrapping -----------------------------------------------------------
#### These SEs will be correct ####

# Initialize weights
ucr_dat_stack$boot_weights <- 1

# Initialize coefficients
did3_rob_coefs <- matrix(0, nrow = nboot, ncol = 2)
did3_ass_coefs <- matrix(0, nrow = nboot, ncol = 2)
did3_lnrob_coefs <- matrix(0, nrow = nboot, ncol = 2)
did3_lnass_coefs <- matrix(0, nrow = nboot, ncol = 2)
did3_lnrob_coefse <- matrix(0, nrow = nboot, ncol = 2)
did3_lnass_coefse <- matrix(0, nrow = nboot, ncol = 2)
did3_rob_coefs_e <- matrix(0, nrow = nboot, ncol = 2)
did3_ass_coefs_e <- matrix(0, nrow = nboot, ncol = 2)
did3_lnrob_coefs_e <- matrix(0, nrow = nboot, ncol = 2)
did3_lnass_coefs_e <- matrix(0, nrow = nboot, ncol = 2)

# This makes sure there is no population weighting
ucr_dat_stack$population <- 1

# for(boot in 1:nboot){
  # 
  # Testing
  boot <- 1
  
  # Progress
  print(boot/nboot)
  
  # Get weights
  # This is a clever alternative to resampling
  if(boot > 1){
    bweights_block <- c(rdirichlet(1, rep(1, length(unique(ucr_dat_stack$state)))))
    
    for(s in 1:length(unique(ucr_dat_stack$state))){
      
      k <- which(ucr_dat_stack$state == unique(ucr_dat_stack$state)[s])
      ucr_dat_stack$boot_weights[k] <- bweights_block[s]
    }
  }
  
  
  
  # 3d Robbery
  mod_3d_rob <- sparse_impute_did(x = sparse_dat_stack,
                                  y = ucr_dat_stack$robbery,
                                  stage2 = treat_effect ~ trip_treat:I(1-disp25) + trip_treat:disp25-1|0|0|0, 
                                  data = ucr_dat_stack, 
                                  treat_var = ucr_dat_stack$trip_treat,
                                  weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)
  
  did3_rob_coefs[boot,1:2] <- mod_3d_rob$coefficients[1:2]
  
  
  
  # 3d Assault
  mod_3d_ass <- sparse_impute_did(x = sparse_dat_stack,
                                  y = ucr_dat_stack$assault,
                                  stage2 = treat_effect ~ trip_treat:I(1-disp25) + trip_treat:disp25-1|0|0|0, 
                                  data = ucr_dat_stack, 
                                  treat_var = ucr_dat_stack$trip_treat,
                                  weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)
  
  did3_ass_coefs[boot,1:2 ] <- mod_3d_ass$coefficients[1:2]
  
  
  # 3d ln Robbery
  mod_3d_lnrob <- sparse_impute_did(x = sparse_dat_stack,
                                  y = ucr_dat_stack$ln_rob,
                                  stage2 = treat_effect ~ trip_treat:I(1-disp25) + trip_treat:disp25-1|0|0|0, 
                                  data = ucr_dat_stack, 
                                  treat_var = ucr_dat_stack$trip_treat,
                                  weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)
  
  did3_lnrob_coefs[boot,1:2] <- mod_3d_lnrob$coefficients[1:2]
  
  
  # 3d ln Assault
  mod_3d_lnass <- sparse_impute_did(x = sparse_dat_stack,
                                  y = ucr_dat_stack$ln_ass,
                                  stage2 = treat_effect ~ trip_treat:I(1-disp25) + trip_treat:disp25-1|0|0|0, 
                                  data = ucr_dat_stack, 
                                  treat_var = ucr_dat_stack$trip_treat,
                                  weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)
  
  did3_lnass_coefs[boot,1:2] <- mod_3d_lnass$coefficients[1:2]
  
  
# Ever Dispensary ---------------------------------------------------------

  
  
  # 3d Robbery
  mod_3d_rob <- sparse_impute_did(x = sparse_dat_stack,
                                  y = ucr_dat_stack$robbery,
                                  stage2 = treat_effect ~ trip_treat + trip_treat:everDisp-1|0|0|0, 
                                  data = ucr_dat_stack, 
                                  treat_var = ucr_dat_stack$trip_treat,
                                  weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)
  
  did3_rob_coefs_e[boot,1:2] <- mod_3d_rob$coefficients[1:2]
  
  
  
  # 3d Assault
  mod_3d_ass <- sparse_impute_did(x = sparse_dat_stack,
                                  y = ucr_dat_stack$assault,
                                  stage2 = treat_effect ~ trip_treat + trip_treat:everDisp-1|0|0|0, 
                                  data = ucr_dat_stack, 
                                  treat_var = ucr_dat_stack$trip_treat,
                                  weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)
  
  did3_ass_coefs_e[boot,1:2 ] <- mod_3d_ass$coefficients[1:2]
  
  
  # 3d ln Robbery
  mod_3d_lnrob <- sparse_impute_did(x = sparse_dat_stack,
                                    y = ucr_dat_stack$ln_rob,
                                    stage2 = treat_effect ~ trip_treat + trip_treat:everDisp-1|0|0|0, 
                                    data = ucr_dat_stack, 
                                    treat_var = ucr_dat_stack$trip_treat,
                                    weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)
  
  did3_lnrob_coefs_e[boot,1:2] <- mod_3d_lnrob$coefficients[1:2]
  
  
  # 3d ln Assault
  mod_3d_lnass <- sparse_impute_did(x = sparse_dat_stack,
                                    y = ucr_dat_stack$ln_ass,
                                    stage2 = treat_effect ~ trip_treat + trip_treat:everDisp-1|0|0|0, 
                                    data = ucr_dat_stack, 
                                    treat_var = ucr_dat_stack$trip_treat,
                                    weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)
  
  did3_lnass_coefs_e[boot,1:2] <- mod_3d_lnass$coefficients[1:2]  
}

did3_rob_coefs %>% colQuantiles(probs = c(.05, .95))
did3_ass_coefs %>% colQuantiles(probs = c(.05, .95))
did3_lnrob_coefs %>% colQuantiles(probs = c(.05, .95))
did3_lnass_coefs %>% colQuantiles(probs = c(.05, .95))

did3_rob_coefs_e %>% colQuantiles(probs = c(.05, .95))
did3_ass_coefs_e %>% colQuantiles(probs = c(.05, .95))
did3_lnrob_coefs_e %>% colQuantiles(probs = c(.05, .95))
did3_lnass_coefs_e %>% colQuantiles(probs = c(.05, .95))

did_disp_results <- list(did3_rob_coefs = did3_rob_coefs,
                         did3_ass_coefs = did3_ass_coefs,
                         did3_lnrob_coefs = did3_lnrob_coefs,
                         did3_lnass_coefs = did3_lnass_coefs,
                         did3_rob_coefs_e = did3_rob_coefs_e,
                         did3_ass_coefs_e = did3_ass_coefs_e,
                         did3_lnrob_coefs_e = did3_lnrob_coefs_e,
                         did3_lnass_coefs_e = did3_lnass_coefs_e)

save(did_disp_results, file = "./Data/did_disp_results.RData")








