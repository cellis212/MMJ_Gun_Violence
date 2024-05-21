# Pre-amble ---------------------------------------------------------------

# Clearing Memory
rm(list = ls())

# Loading Packages
library(magrittr)
library(tidyverse)
library(lfe)

# Set Seed
set.seed(42)
nboot <- 10

# Loading Data ------------------------------------------------------------
# Load Data ---------------------------------------------------------------
load(file = "./Data/Uniform Crime Reports/ucr_dat_specs.RData")


# Save it as "dat"
dat <- ucr_dat_stack


# Options -----------------------------------------------------------------
# These are all examples
Model_Type <- c("DID", "Triple Diff")

Additional_Controls <- c("Yes", "No")

Y_Var_Form <- c("Log", "Linear")

Exclude_Rec <- c("Yes", "No")

Exclude_Long_Treat <- c("Yes", "No")

# Population_Weighted <- c("Yes", "No")



# This creates all combinations of your options
# Change this to your options
Models <- expand.grid(Model_Type = Model_Type, 
                      Additional_Controls = Additional_Controls,
                      Y_Var_Form = Y_Var_Form, 
                      Exclude_Rec = Exclude_Rec, 
                      Exclude_Long_Treat = Exclude_Long_Treat)



# If you have too many, you can run it on a random sample
# samp <- sample(1:nrow(Models), num_samps)
# Models <- Models[samp,]

# Place to save estimates 
# Models$treat <- 0
# Models$ub_treat <- 0
# Models$lb_treat <- 0

# If saving multiple estimates, do it like this

Models$treat_1 <- 0
Models$ub_treat_1 <- 0
Models$lb_treat_1 <- 0


Models$treat_2 <- 0
Models$ub_treat_2 <- 0
Models$lb_treat_2 <- 0
# 
# 
# Models$treat_3 <- 0
# Models$ub_treat_3 <- 0
# Models$lb_treat_3 <- 0
# 
# 
# Models$treat_4 <- 0
# Models$ub_treat_4 <- 0
# Models$lb_treat_4 <- 0

# order must go treat, ub, lb

# Number of models
nOpt <- nrow(Models)
nOpt

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

# Run the Models -------------------------------------------------------------------

# Loop over all options
for(i in 1:nrow(Models)){
  
  # Testing
  # i = 1
  
  # Show progress
  print(i/nOpt)
  
  # Placeholder for subsets
  subs <- 1:nrow(dat) 
  
  # subset for DID
  if(Models$Model_Type[i] == "DID"){
    subs <- subs[which(dat$type[subs] == "gun")]
  }
  
  # Drop super long treat
  if(Models$Exclude_Long_Treat[i] == "Yes"){
    subs <- subs[which(dat$year_since_treat[subs] != "6+")]
  }
  
  # Drop rec states
  if(Models$Exclude_Rec[i] == "Yes"){
    subs <- subs[which(dat$treat_rec[subs] == 0)]
  }
  
  # Create sparse dat
  if(Models$Model_Type[i] == "DID" & 
     Models$Additional_Controls[i] == "Yes") {
    sparse_dat <- sparse.model.matrix(as.formula(paste("", "~ ori + year + 
                                                   stand_your_ground +
                                                   ccw_shall_issue +
                                                   ccw_prohibit +
                                                   backgrd_check3 +
                                                   backgrd_check2 +
                                                   backgrd_check1")), data = dat[subs,])
    
  } 
  
  if(Models$Model_Type[i] == "DID" & 
     Models$Additional_Controls[i] == "No") {
    sparse_dat <- sparse.model.matrix(as.formula(paste("", "~ ori + year")), data = dat[subs,])
  }
  
  if(Models$Model_Type[i] == "Triple Diff" & 
     Models$Additional_Controls[i] == "Yes") {
    next()
  }
  
  if(Models$Model_Type[i] == "Triple Diff" & 
     Models$Additional_Controls[i] == "No") {
    sparse_dat <- sparse.model.matrix(as.formula(paste("", "~ type_ori + type_year + year_ori")), data = dat[subs,])
  }
  
  # Get correct yvar
  if(Models$Y_Var_Form[i] == "Log" &
     Models$Model_Type[i] == "DID"){
    yvar_rob <- dat$ln_rob_gun[subs]
    yvar_ass <- dat$ln_ass_gun[subs]
  }
  
  if(Models$Y_Var_Form[i] == "Linear" &
     Models$Model_Type[i] == "DID"){
    yvar_rob <- dat$actual_robbery_with_a_gun[subs]
    yvar_ass <- dat$actual_assault_with_a_gun[subs]
  }
  
  if(Models$Y_Var_Form[i] == "Log" &
     Models$Model_Type[i] == "Triple Diff"){
    yvar_rob <- dat$ln_rob[subs]
    yvar_ass <- dat$ln_ass[subs]
  }
  
  if(Models$Y_Var_Form[i] == "Linear" &
     Models$Model_Type[i] == "Triple Diff"){
    yvar_rob <- dat$robbery[subs]
    yvar_ass <- dat$assault[subs]
  }
  
  
  
  
  # Estimating the Models ------------------------------------------------------
  # Creating the formula using above
  
  coefs_rob <- c()
  coefs_ass <- c()
  dat$boot_weights[subs] <- 1/length(subs)
  
  
  
  # Bootstrap
  for(boot in 1:nboot){
    
    # print(boot/nboot)
    # Get weights
    # This is a clever alternative to resampling
    if(boot > 1){
      bweights_block <- c(rdirichlet(1, rep(1, length(unique(dat$state[subs])))))
      
      for(s in 1:length(unique(dat$state[subs]))){
        
        k <- which(dat$state[subs] == unique(dat$state[subs])[s])
        dat$boot_weights[subs][k] <- bweights_block[s]
        
      }
    }
  
    
    if(Models$Model_Type[i] == "DID") {
      # Robbery
      mod_rob <- sparse_impute_did(x = sparse_dat,
                                   y = yvar_rob,
                                   stage2 = treat_effect ~ treat - 1|0|0|0, 
                                   data = dat[subs,], 
                                   treat_var = dat$treat[subs],
                                   weights = dat$boot_weights[subs])
      
      # Robbery
      mod_ass <- sparse_impute_did(x = sparse_dat,
                                   y = yvar_ass,
                                   stage2 = treat_effect ~ treat - 1|0|0|0, 
                                   data = dat[subs,], 
                                   treat_var = dat$treat[subs],
                                   weights = dat$boot_weights[subs])
    } else {
      # Robbery
      mod_rob <- sparse_impute_did(x = sparse_dat,
                                   y = yvar_rob,
                                   stage2 = treat_effect ~ trip_treat - 1|0|0|0, 
                                   data = dat[subs,], 
                                   treat_var = dat$trip_treat[subs],
                                   weights = dat$boot_weights[subs])
      
      # Robbery
      mod_ass <- sparse_impute_did(x = sparse_dat,
                                   y = yvar_ass,
                                   stage2 = treat_effect ~ trip_treat - 1|0|0|0, 
                                   data = dat[subs,], 
                                   treat_var = dat$trip_treat[subs],
                                   weights = dat$boot_weights[subs])
      
    }
    
    if(Models$Y_Var_Form[i] == "Log") {
      
      coefs_rob[boot] <- mod_rob$coefficients[1]
      coefs_ass[boot] <- mod_ass$coefficients[1]
      
    } else {
      
      coefs_rob[boot] <- mod_rob$coefficients[1]/weighted.mean(yvar_rob, dat$boot_weights[subs], na.rm = T)
      coefs_ass[boot] <- mod_ass$coefficients[1]/weighted.mean(yvar_ass, dat$boot_weights[subs], na.rm = T)
      
    }
  }
  
  # Save with 90% confidence intervals
  Models$treat_1[i] <- coefs_rob[1]  
  Models$lb_treat_1[i] <- quantile(coefs_rob, .05)
  Models$ub_treat_1[i] <- quantile(coefs_rob, .95)
  
  Models$treat_2[i] <- coefs_ass[1]  
  Models$lb_treat_2[i] <- quantile(coefs_ass, .05)
  Models$ub_treat_2[i] <- quantile(coefs_ass, .95)
}


# Saving 
save(Models, file = "./Data/spec_curve_models.RData")

