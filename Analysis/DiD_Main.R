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

# Split of unstacked data
ucr_dat <- ucr_dat_stack %>% filter(type == "gun")

# Clean Data --------------------------------------------------------------

# Make sure year is a factor
ucr_dat$year <- as.character(ucr_dat$year)
ucr_dat$year_num <- as.numeric(as.character(ucr_dat$year))
ucr_dat_stack$year <- as.character(ucr_dat_stack$year)





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

# Base Diffs --------------------------------------------------------------
sparse_dat <- sparse.model.matrix(as.formula(paste("", "~ ori + year + 
                                                   stand_your_ground +
                                                   ccw_shall_issue +
                                                   ccw_prohibit +
                                                   backgrd_check3 +
                                                   backgrd_check2 +
                                                   backgrd_check1")), data = ucr_dat)


# Triple Diffs -------------------------------------------------------------
sparse_dat_stack <- sparse.model.matrix(as.formula(paste("", "~ type_ori + type_year + year_ori")), data = ucr_dat_stack)




# Bootstrapping -----------------------------------------------------------
#### These SEs will be correct ####

# Initialize weights
ucr_dat$boot_weights <- 1
ucr_dat_stack$boot_weights <- 1

# Initialize coefficients
did_rob_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_ass_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_lnrob_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_lnass_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_rob_K_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_ass_K_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_lnrob_K_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_lnass_K_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_rob_u_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_ass_u_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_lnrob_u_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_lnass_u_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_rob_o_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_ass_o_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_lnrob_o_coefs <- matrix(0, nrow = nboot, ncol = 11)
did_lnass_o_coefs <- matrix(0, nrow = nboot, ncol = 11)
did3_rob_coefs <- matrix(0, nrow = nboot, ncol = 11)
did3_ass_coefs <- matrix(0, nrow = nboot, ncol = 11)
did3_lnrob_coefs <- matrix(0, nrow = nboot, ncol =11)
did3_lnass_coefs <- matrix(0, nrow = nboot, ncol = 11)

# This makes sure there is no population weighting
ucr_dat$population <- 1
ucr_dat_stack$population <- 1

for(boot in 1:nboot){
# 
# Testing
# boot <- 2

# Progress
print(boot/nboot)

# Get weights
# This is a clever alternative to resampling
if(boot > 1){
  bweights_block <- c(rdirichlet(1, rep(1, length(unique(ucr_dat$state)))))
  
  for(s in 1:length(unique(ucr_dat$state))){
    
    l <- which(ucr_dat$state == unique(ucr_dat$state)[s])
    ucr_dat$boot_weights[l] <- bweights_block[s]
    
    k <- which(ucr_dat_stack$state == unique(ucr_dat$state)[s])
    ucr_dat_stack$boot_weights[k] <- bweights_block[s]
  }
}


# DID - GUN ---------------------------------------------------------------


# Robbery
mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_robbery_with_a_gun,
                                 stage2 = treat_effect ~ treat - 1|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_rob_coefs[boot,1] <- mod_did_rob$coefficients[1]

mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_robbery_with_a_gun,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_rob_coefs[boot,2:11] <- mod_did_rob$coefficients[2:11]


# Assault
mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_assault_with_a_gun,
                                 stage2 = treat_effect ~ treat - 1|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_ass_coefs[boot,1] <- mod_did_ass$coefficients[1]

mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_assault_with_a_gun,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_ass_coefs[boot,2:11] <- mod_did_ass$coefficients[2:11]

# ln Robbery
mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_rob_gun,
                                 stage2 = treat_effect ~ treat - 1|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_lnrob_coefs[boot, 1 ] <- mod_did_rob$coefficients[1]

mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_rob_gun,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_lnrob_coefs[boot,2:11] <- mod_did_rob$coefficients[2:11]


# ln Assault
mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_ass_gun,
                                 stage2 = treat_effect ~ treat - 1|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_lnass_coefs[boot,1 ] <- mod_did_ass$coefficients[1]

mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_ass_gun,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_lnass_coefs[boot,2:11] <- mod_did_ass$coefficients[2:11]


# DID - Knife -------------------------------------------------------------


# Robbery
mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_robbery_with_a_knife,
                                 stage2 = treat_effect ~ treat - 1|0|0|0,  
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_rob_K_coefs[boot,1] <- mod_did_rob$coefficients[1]

mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_robbery_with_a_knife,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_rob_K_coefs[boot,2:11] <- mod_did_rob$coefficients[2:11]


# Assault
mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_assault_with_a_knife,
                                 stage2 = treat_effect ~ treat - 1|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_ass_K_coefs[boot,1] <- mod_did_ass$coefficients[1]

mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_assault_with_a_knife,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_ass_K_coefs[boot,2:11] <- mod_did_ass$coefficients[2:11]

# Logs 
# Robbery
mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_rob_knife,
                                 stage2 = treat_effect ~ treat - 1|0|0|0,  
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_lnrob_K_coefs[boot,1] <- mod_did_rob$coefficients[1]

mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_rob_knife,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_lnrob_K_coefs[boot,2:11] <- mod_did_rob$coefficients[2:11]


# Assault
mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_ass_knife,
                                 stage2 = treat_effect ~ treat - 1|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_lnass_K_coefs[boot,1] <- mod_did_ass$coefficients[1]

mod_lndid_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_ass_knife,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_lnass_K_coefs[boot,2:11] <- mod_lndid_ass$coefficients[2:11]


# DID - Unarmed -------------------------------------------------------------


# Robbery
mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_robbery_unarmed,
                                 stage2 = treat_effect ~ treat - 1|0|0|0,  
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_rob_u_coefs[boot,1] <- mod_did_rob$coefficients[1]

mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_robbery_unarmed,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_rob_u_coefs[boot,2:11] <- mod_did_rob$coefficients[2:11]


# Assault
mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_assault_unarmed,
                                 stage2 = treat_effect ~ treat - 1|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_ass_u_coefs[boot,1] <- mod_did_ass$coefficients[1]

mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_assault_unarmed,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_ass_u_coefs[boot,2:11] <- mod_did_ass$coefficients[2:11]

# Logs 
# Robbery
mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_rob_unarmed,
                                 stage2 = treat_effect ~ treat - 1|0|0|0,  
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_lnrob_u_coefs[boot,1] <- mod_did_rob$coefficients[1]

mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_rob_unarmed,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_lnrob_u_coefs[boot,2:11] <- mod_did_rob$coefficients[2:11]


# Assault
mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_ass_unarmed,
                                 stage2 = treat_effect ~ treat - 1|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_lnass_u_coefs[boot,1] <- mod_did_ass$coefficients[1]

mod_lndid_ass <- sparse_impute_did(x = sparse_dat,
                                   y = ucr_dat$ln_ass_unarmed,
                                   stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                   data = ucr_dat, 
                                   treat_var = ucr_dat$treat,
                                   weights = ucr_dat$boot_weights * ucr_dat$population)


did_lnass_u_coefs[boot,2:11] <- mod_lndid_ass$coefficients[2:11]


# DID - Other Weapon -------------------------------------------------------------


# Robbery
mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_robbery_other_weapon,
                                 stage2 = treat_effect ~ treat - 1|0|0|0,  
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_rob_o_coefs[boot,1] <- mod_did_rob$coefficients[1]

mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_robbery_other_weapon,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_rob_o_coefs[boot,2:11] <- mod_did_rob$coefficients[2:11]


# Assault
mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_assault_other_weapon,
                                 stage2 = treat_effect ~ treat - 1|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_ass_o_coefs[boot,1] <- mod_did_ass$coefficients[1]

mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$actual_assault_other_weapon,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_ass_o_coefs[boot,2:11] <- mod_did_ass$coefficients[2:11]

# Logs 
# Robbery
mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_rob_other_weapon,
                                 stage2 = treat_effect ~ treat - 1|0|0|0,  
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_lnrob_o_coefs[boot,1] <- mod_did_rob$coefficients[1]

mod_did_rob <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_rob_other_weapon,
                                 stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)

did_lnrob_o_coefs[boot,2:11] <- mod_did_rob$coefficients[2:11]


# Assault
mod_did_ass <- sparse_impute_did(x = sparse_dat,
                                 y = ucr_dat$ln_ass_other_weapon,
                                 stage2 = treat_effect ~ treat - 1|0|0|0, 
                                 data = ucr_dat, 
                                 treat_var = ucr_dat$treat,
                                 weights = ucr_dat$boot_weights * ucr_dat$population)


did_lnass_o_coefs[boot,1] <- mod_did_ass$coefficients[1]

mod_lndid_ass <- sparse_impute_did(x = sparse_dat,
                                   y = ucr_dat$ln_ass_other_weapon,
                                   stage2 = treat_effect ~ year_since_treat|0|0|0, 
                                   data = ucr_dat, 
                                   treat_var = ucr_dat$treat,
                                   weights = ucr_dat$boot_weights * ucr_dat$population)


did_lnass_o_coefs[boot,2:11] <- mod_lndid_ass$coefficients[2:11]
# Trip Diffs --------------------------------------------------------------



# 3d Robbery
mod_3d_rob <- sparse_impute_did(x = sparse_dat_stack,
                                y = ucr_dat_stack$robbery,
                                stage2 = treat_effect ~ trip_treat - 1|0|0|0, 
                                data = ucr_dat_stack, 
                                treat_var = ucr_dat_stack$trip_treat,
                                weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)

did3_rob_coefs[boot,1] <- mod_3d_rob$coefficients[1]

mod_3d_rob <- sparse_impute_did(x = sparse_dat_stack,
                                y = ucr_dat_stack$robbery,
                                stage2 = treat_effect ~ trip_year_since_treat|0|0|0, 
                                data = ucr_dat_stack, 
                                treat_var = ucr_dat_stack$trip_treat,
                                weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)

did3_rob_coefs[boot,2:11] <- mod_3d_rob$coefficients[2:11]


# 3d Assault
mod_3d_ass <- sparse_impute_did(x = sparse_dat_stack,
                                y = ucr_dat_stack$assault,
                                stage2 = treat_effect ~ trip_treat - 1|0|0|0, 
                                data = ucr_dat_stack, 
                                treat_var = ucr_dat_stack$trip_treat,
                                weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)

did3_ass_coefs[boot,1] <- mod_3d_ass$coefficients[1]

mod_3d_ass <- sparse_impute_did(x = sparse_dat_stack,
                                y = ucr_dat_stack$assault,
                                stage2 = treat_effect ~ trip_year_since_treat|0|0|0, 
                                data = ucr_dat_stack, 
                                treat_var = ucr_dat_stack$trip_treat,
                                weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)

did3_ass_coefs[boot,2:11] <- mod_3d_ass$coefficients[2:11]

# 3d ln Robbery
mod_3d_rob <- sparse_impute_did(x = sparse_dat_stack,
                                y = ucr_dat_stack$ln_rob,
                                stage2 = treat_effect ~ trip_treat - 1|0|0|0, 
                                data = ucr_dat_stack, 
                                treat_var = ucr_dat_stack$trip_treat,
                                weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)

did3_lnrob_coefs[boot,1] <- mod_3d_rob$coefficients[1]

mod_3d_rob <- sparse_impute_did(x = sparse_dat_stack,
                                y = ucr_dat_stack$ln_rob,
                                stage2 = treat_effect ~ trip_year_since_treat|0|0|0, 
                                data = ucr_dat_stack, 
                                treat_var = ucr_dat_stack$trip_treat,
                                weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)

did3_lnrob_coefs[boot,2:11] <- mod_3d_rob$coefficients[2:11]


# 3d ln Assault
mod_3d_ass <- sparse_impute_did(x = sparse_dat_stack,
                                y = ucr_dat_stack$ln_ass,
                                stage2 = treat_effect ~ trip_treat - 1|0|0|0, 
                                data = ucr_dat_stack, 
                                treat_var = ucr_dat_stack$trip_treat,
                                weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)

did3_lnass_coefs[boot,1] <- mod_3d_ass$coefficients[1]

mod_3d_ass <- sparse_impute_did(x = sparse_dat_stack,
                                y = ucr_dat_stack$ln_ass,
                                stage2 = treat_effect ~ trip_year_since_treat|0|0|0, 
                                data = ucr_dat_stack, 
                                treat_var = ucr_dat_stack$trip_treat,
                                weights = ucr_dat_stack$boot_weights * ucr_dat_stack$population)

did3_lnass_coefs[boot,2:11] <- mod_3d_ass$coefficients[2:11]

}

did_rob_coefs  %>% colQuantiles(probs = c(.05, .95))
did_ass_coefs %>% colQuantiles(probs = c(.05, .95))
did_lnrob_coefs  %>% colQuantiles(probs = c(.05, .95))
did_lnass_coefs %>% colQuantiles(probs = c(.05, .95))
did_rob_K_coefs  %>% colQuantiles(probs = c(.05, .95))
did_ass_K_coefs %>% colQuantiles(probs = c(.05, .95))
did_lnrob_K_coefs  %>% colQuantiles(probs = c(.05, .95))
did_lnass_K_coefs %>% colQuantiles(probs = c(.05, .95))
did_rob_u_coefs  %>% colQuantiles(probs = c(.05, .95))
did_ass_u_coefs %>% colQuantiles(probs = c(.05, .95))
did_lnrob_u_coefs  %>% colQuantiles(probs = c(.05, .95))
did_lnass_u_coefs %>% colQuantiles(probs = c(.05, .95))
did_rob_o_coefs  %>% colQuantiles(probs = c(.05, .95))
did_ass_o_coefs %>% colQuantiles(probs = c(.05, .95))
did_lnrob_o_coefs  %>% colQuantiles(probs = c(.05, .95))
did_lnass_o_coefs %>% colQuantiles(probs = c(.05, .95))
did3_rob_coefs %>% colQuantiles(probs = c(.05, .95))
did3_ass_coefs %>% colQuantiles(probs = c(.05, .95))
did3_lnrob_coefs %>% colQuantiles(probs = c(.05, .95))
did3_lnass_coefs %>% colQuantiles(probs = c(.05, .95))

did_main_results <- list(did_rob_coefs = did_rob_coefs,
                         did_ass_coefs = did_ass_coefs,
                         did_lnrob_coefs = did_lnrob_coefs,
                         did_lnass_coefs = did_lnass_coefs,
                         did_rob_K_coefs = did_rob_K_coefs,
                         did_ass_K_coefs = did_ass_K_coefs,
                         did_lnrob_K_coefs = did_lnrob_K_coefs,
                         did_lnass_K_coefs = did_lnass_K_coefs,
                         did_rob_u_coefs = did_rob_u_coefs,
                         did_ass_u_coefs = did_ass_u_coefs,
                         did_lnrob_u_coefs = did_lnrob_u_coefs,
                         did_lnass_u_coefs = did_lnass_u_coefs,
                         did_rob_o_coefs = did_rob_o_coefs,
                         did_ass_o_coefs = did_ass_o_coefs,
                         did_lnrob_o_coefs = did_lnrob_o_coefs,
                         did_lnass_o_coefs = did_lnass_o_coefs,
                         did3_rob_coefs = did3_rob_coefs,
                         did3_ass_coefs = did3_ass_coefs,
                         did3_lnrob_coefs = did3_lnrob_coefs,
                         did3_lnass_coefs = did3_lnass_coefs)

save(did_main_results, file = "./Data/did_main_results.RData")








