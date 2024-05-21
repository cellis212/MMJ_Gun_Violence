############################################
#
# This file runs the DiDs
#
#############################################



# Preamble ----------------------------------------------------------------

# Clearing the workspace
rm(list = ls())

# Loading the packages that will be needed. 
require(tidyverse)
require(Hmisc)
require(gtools)
library(MatrixModels)
library(matrixStats)
library(Matrix)
library(lfe)
library(MatchIt)
library(gsynth)



# Install Honest DiD Package
# Install remotes package if not installed
# install.packages("remotes") 
# 
# # Turn off warning-error-conversion, because the tiniest warning stops installation
# Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
# 
# # install from github
# remotes::install_github("asheshrambachan/HonestDiD")

# Load in functions and options file.
source(file = "./Code/Functions_and_Options.R")

# Load Data ---------------------------------------------------------------
load(file = "./Data/ucr_data_stack_cleaned_subset.RData")
load(file = "./Data/did_main_results.RData")
load(file = "./Data/did_first_results.RData")
load(file = "./Data/sui_data_stack_cleaned_2.RData")
load(file = "./Data/sui_data_cleaned_2.RData")
load(file = "./Data/sui_data_cleaned_2_o.RData")

# Split of unstacked data
ucr_stack <- ucr_dat_stack
ucr_dat <- ucr_dat_stack %>% filter(type == "gun")



#  Formatting Functions ---------------------------------------------------

# Coefficient with stars from bootstrap
DIDcoef <- function(model, digits = 3, ci = TRUE, alpha = .05){
  coef <- model[1]
  
  if(ci == FALSE){
  ses <- sd(model)
  
  if(abs(coef/ses) >= 2.58){
    
    return(paste0(round(coef, digits = digits), "^{***}")) 
    
  } else if(abs(coef/ses) >= 1.96 & abs(coef/ses) < 2.58) {
    
    return(paste0(round(coef, digits = digits),"^{**}"))
    
  } else if(abs(coef/ses) >= 1.645 & abs(coef/ses) < 1.96) {
    
    return(paste0(round( coef, digits = digits),"^{*}"))
    
  } else {
    
    return(paste0(round( coef, digits = digits),"^{}"))
  }
  } else {
    ci <- quantile(model, probs = c(alpha, 1-alpha))
    if(sign(ci[1]) == sign(ci[2])){
      return(paste0(round(coef, digits = digits), "^{*}")) 
    } else {
      return(paste0(round( coef, digits = digits),"^{}"))
    }
  }
    
}


# Standard Error pull from bootstrap
DIDse <- function(model, digits = 3){
  
  ses <- sd(model)
  
  if(round(ses,digits = digits) > 0) {
    return(paste0("(", round(ses, digits = digits), ")"))
  } else {
    return(paste0("(<", 1/10^digits, ")"))
  }
}


# Confidence Interval pull from bootstrap
DIDci <- function(model, digits = 3, alpha = .05){
  ci <- quantile(model, probs = c(alpha, 1-alpha))
  return(paste0("[", round(ci[1], digits = digits), " -- ", round(ci[2], digits = digits), "]"))
}

# General formatting
f <- function(x, digits = 3){
  if(x==0){ "0"}
  if(x>100){  format(x, digits = 3, nsmall = 0, big.mark = ",") 
  } else {
    format(x, digits = digits, nsmall = digits, big.mark = ",") 
  }
}


# Confidence Interval pull from bootstrap
DIDci_per <- function(model, digits = 3, alpha = .25, m = 1){
  ci <- 100*quantile(model, probs = c(alpha, 1-alpha))
  return(paste0("[", round(ci[1]/m, digits = digits), "\\% -- ", round(ci[2]/m, digits = digits), "\\%]"))
}


pre_trend_slope <- function(model, digits = 3){
  num_periods <- ncol(model)
  x <- 1:num_periods
  y <- model[1,1:num_periods]
  weights <- 1/colSds(model)
  weights <- weights/sum(weights)
  
  mod <- lm(y~x, weights = weights)
  return(f(mod$coefficients[2]))
}

pre_trend_slope_ci <- function(model, digits = 3){
  num_periods <- ncol(model)
  x <- 1:num_periods
  y <- model[1,1:num_periods]
  weights <- 1/colSds(model)
  weights <- weights/sum(weights)
  
  mod <- lm(y~x, weights = weights) %>% summary()
  return(paste0("[", round(mod$coefficients[2,1] -1.96*mod$coefficients[2,2], digits = digits), " -- ", round(mod$coefficients[2,1] + 1.96*mod$coefficients[2,2], digits = digits),"]"))
}


# get_vcov <- function(b_df)
# {
#   mat <- matrix(0, ncol = ncol(b_df), nrow = ncol(b_df))
#   
#   for(i in 1:ncol(b_df)){
#     for(j in 1:ncol(b_df)){
#       mat[i,j] <- cov(b_df[,i], b_df[,j])
#     }
#   }
#   
#   return(mat)
# }


# # Credible pull from bootstrap (main coefficient goes jn first column)
# DIDcred <- function(model, digits = 2, num_pre = 3, num_post = 6, l = c(0,0,0,1,0,0)){
#   
#   sigma <- get_vcov(model[-1,])
#   betahat <- model[1,]
#   
#   a <- 
#     HonestDiD::createSensitivityResults(betahat = betahat,
#                                         sigma = sigma, 
#                                         numPrePeriods = num_pre,
#                                         numPostPeriods = num_post,
#                                         l_vec = l) %>% 
#     as.data.frame()
#   return(paste(f(a[1,1]), "--", f(a[1,2])))
# } 

# pre_slope <- function(model, digits = 2, num_pre = 4){
#   pre_coefs <- model[1,]




# DIDcred(model = did_main_results$did_rob_coefs[,-c(1)],
#         num_pre = 4,
#         num_post = 5,
#         l = rep(1/5, 5))

# Base Diffs --------------------------------------------------------------
l <- which(ucr_dat$year_since_treat == "-1")
k <- which(sui_dat$year_since_treat == "-1")

# Don't weight by pop
ucr_dat$population <- 1
ucr_stack$population <- 1
sui_dat$total_population <- 1
sui_dat_stack$total_population <- 1
sui_dat_o$total_population <- 1

# Get max pre-trend
max_pt <- c(
  which.max(abs(did_main_results$did_rob_coefs[1, 3:5])),
  which.max(abs(did_main_results$did_ass_coefs[1, 3:5])),
  which.max(abs(did_main_results$did_lnrob_coefs[1, 3:5])),
  which.max(abs(did_main_results$did_lnass_coefs[1, 3:5])))

# Raw Vals
cat("\\textbf{Panel (A): Gun Crime:} &&&&\\\\ \\\\",
    # "Estimated ATT: &", DIDcoef(did_main_results$did_rob_coefs[,1]), 
    # "&", DIDcoef(did_main_results$did_ass_coefs[,1]),
    # "&", DIDcoef(did_main_results$did_lnrob_coefs[,1]),
    # "&", DIDcoef(did_main_results$did_lnass_coefs[,1]),
    # "\\\\",
    # "&", DIDci(did_main_results$did_rob_coefs[,1]), 
    # "&", DIDci(did_main_results$did_ass_coefs[,1]),
    # "&", DIDci(did_main_results$did_lnrob_coefs[,1]),
    # "&", DIDci(did_main_results$did_lnass_coefs[,1]),
    "Estimated ATT: &", DIDcoef(did_main_results$did_rob_coefs[,1]), 
    "&", DIDcoef(did_main_results$did_ass_coefs[,1]),
    "&", DIDcoef(did_main_results$did_lnrob_coefs[,1]),
    "&", DIDcoef(did_main_results$did_lnass_coefs[,1]),
    "\\\\",
    "&", DIDci(did_main_results$did_rob_coefs[,1]), 
    "&", DIDci(did_main_results$did_ass_coefs[,1]),
    "&", DIDci(did_main_results$did_lnrob_coefs[,1]),
    "&", DIDci(did_main_results$did_lnass_coefs[,1]),
    "\\\\ \\\\",    
    "Pre-Legalization Mean ",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_dat$actual_robbery_with_a_gun[l], ucr_dat$population[l])), "}",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_dat$actual_assault_with_a_gun[l], ucr_dat$population[l])), "}",
    "& \\multicolumn{1}{c}{-- }",
    "& \\multicolumn{1}{c}{-- }",
    "\\\\ \\\\",
    "Implied Percent Change ",
    "& \\multicolumn{1}{c}{", f(100*did_main_results$did_rob_coefs[1,1]/weighted.mean(ucr_dat$actual_robbery_with_a_gun[l], ucr_dat$population[l])), "\\%}",
    "& \\multicolumn{1}{c}{", f(100*did_main_results$did_ass_coefs[1,1]/weighted.mean(ucr_dat$actual_assault_with_a_gun[l], ucr_dat$population[l])), "\\%}",
    "& \\multicolumn{1}{c}{ --}",
    "& \\multicolumn{1}{c}{--}",
    # "\\\\ ",
    # "&", DIDci_per(did_main_results$did_rob_coefs[,1], m = weighted.mean(ucr_dat$actual_robbery_with_a_gun[l], ucr_dat$population[l])),
    # "&", DIDci_per(did_main_results$did_ass_coefs[,1], m = weighted.mean(ucr_dat$actual_assault_with_a_gun[l], ucr_dat$population[l])),
    # "&", DIDci_per(did_main_results$did_lnrob_coefs[,1], m = weighted.mean(sui_dat$num_homi[k], sui_dat$total_population[k])),
    # "&", DIDci_per(did_main_results$did_lnass_coefs[,1], m = weighted.mean(sui_dat$num_sui[k], sui_dat$total_population[k])),
    "\\\\ \\\\",
    "Average Pre-trend Violation: &", f(weighted.mean(did_main_results$did_rob_coefs[1,2:5], 1/colSds(did_main_results$did_rob_coefs[,2:5]))), 
    "&", f(weighted.mean(did_main_results$did_ass_coefs[1,2:5], 1/colSds(did_main_results$did_ass_coefs[,2:5]))),
    "&", f(weighted.mean(did_main_results$did_lnrob_coefs[1,2:5], 1/colSds(did_main_results$did_lnrob_coefs[,2:5]))),
    "&", f(weighted.mean(did_main_results$did_lnass_coefs[1,2:5], 1/colSds(did_main_results$did_lnass_coefs[,2:5]))),
    "\\\\ \\\\",    
    # "SD of Pre-trend Violations: &", f(sqrt(wtd.var(did_main_results$did_rob_coefs[1,2:5], 1/colSds(did_main_results$did_rob_coefs[,2:5]), normwt = TRUE))), 
    # "&", f(sqrt(wtd.var(did_main_results$did_ass_coefs[1,2:5], 1/colSds(did_main_results$did_ass_coefs[,2:5]), normwt = TRUE))),
    # "&", f(sqrt(wtd.var(did_main_results$did_lnrob_coefs[1,2:5], 1/colSds(did_main_results$did_lnrob_coefs[,2:5]), normwt = TRUE))),
    # "&", f(sqrt(wtd.var(did_main_results$did_lnass_coefs[1,2:5], 1/colSds(did_main_results$did_lnass_coefs[,2:5]), normwt = TRUE))),
    # "\\\\ \\\\",
    "Slope of Pre-trend Violations: &", pre_trend_slope(did_main_results$did_rob_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_ass_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_lnrob_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_lnass_coefs[,2:5]), 
    "\\\\",
    "&", pre_trend_slope_ci(did_main_results$did_rob_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_ass_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_lnrob_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_lnass_coefs[,2:5]),
    "\\\\ \\\\",
    "Largest Pre-trend Violation: &", DIDcoef(did_main_results$did_rob_coefs[,(max_pt[1]+2)]), 
    "&", DIDcoef(did_main_results$did_ass_coefs[,(max_pt[2]+2)]),
    "&", DIDcoef(did_main_results$did_lnrob_coefs[,(max_pt[3]+2)]),
    "&", DIDcoef(did_main_results$did_lnass_coefs[,(max_pt[4]+2)]),
    "\\\\",
    "&", DIDci(did_main_results$did_rob_coefs[,(max_pt[1]+2)]), 
    "&", DIDci(did_main_results$did_ass_coefs[,(max_pt[2]+2)]),
    "&", DIDci(did_main_results$did_lnrob_coefs[,(max_pt[3]+2)]),
    "&", DIDci(did_main_results$did_lnass_coefs[,(max_pt[4]+2)]),
    "\\\\ \\\\ ",

    
    "\\textbf{Panel (B): Knife Crime:} &&&&\\\\ \\\\",
    "Estimated ATT: &", DIDcoef(did_main_results$did_rob_K_coefs[,1]), 
    "&", DIDcoef(did_main_results$did_ass_K_coefs[,1]),
    "&", DIDcoef(did_main_results$did_lnrob_K_coefs[,1]),
    "&", DIDcoef(did_main_results$did_lnass_K_coefs[,1]),
    "\\\\",
    "&", DIDci(did_main_results$did_rob_K_coefs[,1]), 
    "&", DIDci(did_main_results$did_ass_K_coefs[,1]),
    "&", DIDci(did_main_results$did_lnrob_K_coefs[,1]),
    "&", DIDci(did_main_results$did_lnass_K_coefs[,1]),
    "\\\\ \\\\",  
    "Pre-Legalization Mean ",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_dat$actual_robbery_with_a_knife[l], ucr_dat$population[l])), "}",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_dat$actual_assault_with_a_knife[l], ucr_dat$population[l])), "}",
    "& \\multicolumn{1}{c}{--}",
    "& \\multicolumn{1}{c}{--}",
    "\\\\ \\\\",
    "Implied Percent Change ",
    "& \\multicolumn{1}{c}{", f(100*did_main_results$did_rob_K_coefs[1,1]/weighted.mean(ucr_dat$actual_robbery_with_a_knife[l], ucr_dat$population[l])), "\\%}",
    "& \\multicolumn{1}{c}{", f(100*did_main_results$did_ass_K_coefs[1,1]/weighted.mean(ucr_dat$actual_assault_with_a_knife[l], ucr_dat$population[l])), "\\%}",
    "& \\multicolumn{1}{c}{--}",
    "& \\multicolumn{1}{c}{--}",
    # "\\\\ ",
    # "&", DIDci_per(did_main_results$did_rob_K_coefs[,1], m = weighted.mean(ucr_dat$actual_robbery_with_a_knife[l], ucr_dat$population[l])), 
    # "&", DIDci_per(did_main_results$did_ass_K_coefs[,1], m = weighted.mean(ucr_dat$actual_assault_with_a_knife[l], ucr_dat$population[l])),
    # "&", DIDci_per(did_main_results$did_lnrob_K_coefs[,1], m = weighted.mean(sui_dat_o$num_homi[k], sui_dat$total_population[k])),
    # "&", DIDci_per(did_main_results$did_lnass_K_coefs[,1], m = weighted.mean(sui_dat_o$num_sui[k], sui_dat$total_population[k])),
    "\\\\ \\\\",
    "Average Pre-trend Violation: &", f(weighted.mean(did_main_results$did_rob_K_coefs[1,2:5], 1/colSds(did_main_results$did_rob_K_coefs[,2:5]))), 
    "&", f(weighted.mean(did_main_results$did_ass_K_coefs[1,2:5], 1/colSds(did_main_results$did_ass_K_coefs[,2:5]))),
    "&", f(weighted.mean(did_main_results$did_lnrob_K_coefs[1,2:5], 1/colSds(did_main_results$did_lnrob_K_coefs[,2:5]))),
    "&", f(weighted.mean(did_main_results$did_lnass_K_coefs[1,2:5], 1/colSds(did_main_results$did_lnass_K_coefs[,2:5]))),
    "\\\\ \\\\",    
    # "SD of Pre-trend Violations: &", f(sqrt(wtd.var(did_main_results$did_rob_K_coefs[1,2:5], 1/colSds(did_main_results$did_rob_K_coefs[,2:5]), normwt = TRUE))), 
    # "&", f(sqrt(wtd.var(did_main_results$did_ass_K_coefs[1,2:5], 1/colSds(did_main_results$did_ass_K_coefs[,2:5]), normwt = TRUE))),
    # "&", f(sqrt(wtd.var(did_main_results$did_lnrob_K_coefs[1,2:5], 1/colSds(did_main_results$did_lnrob_K_coefs[,2:5]), normwt = TRUE))),
    # "&", f(sqrt(wtd.var(did_main_results$did_lnass_K_coefs[1,2:5], 1/colSds(did_main_results$did_lnass_K_coefs[,2:5]), normwt = TRUE))),
    # "\\\\ \\\\",
    "Slope of Pre-trend Violations: &", pre_trend_slope(did_main_results$did_rob_K_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_ass_K_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_lnrob_K_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_lnass_K_coefs[,2:5]), 
    "\\\\",
    "&", pre_trend_slope_ci(did_main_results$did_rob_K_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_ass_K_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_lnrob_K_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_lnass_K_coefs[,2:5]),
    "\\\\ \\\\",
    "Largest Pre-trend Violation: &", DIDcoef(did_main_results$did_rob_K_coefs[,(max_pt[1]+2)]), 
    "&", DIDcoef(did_main_results$did_ass_K_coefs[,(max_pt[2]+2)]),
    "&", DIDcoef(did_main_results$did_lnrob_K_coefs[,(max_pt[3]+2)]),
    "&", DIDcoef(did_main_results$did_lnass_K_coefs[,(max_pt[4]+2)]),
    "\\\\",
    "&", DIDci(did_main_results$did_rob_K_coefs[,(max_pt[1]+2)]), 
    "&", DIDci(did_main_results$did_ass_K_coefs[,(max_pt[2]+2)]),
    "&", DIDci(did_main_results$did_lnrob_K_coefs[,(max_pt[3]+2)]),
    "&", DIDci(did_main_results$did_lnass_K_coefs[,(max_pt[4]+2)]),
    "\\\\ \\\\",

    " \\hline ",
    "\\\\[0.5mm]",
    "Gun Law Controls? & Yes & Yes & Yes & Yes",
    "\\\\[0.5mm]",
    "Unit of Observation ",
    "& \\multicolumn{4}{c}{ Jurisdiction - Year}",
    "\\\\[0.5mm]",
    "Number of Observations ",
    "& \\multicolumn{4}{c}{", f(nrow(ucr_dat)), "}",
    "\\\\[0.5mm]", 
    "\\bottomrule",
    sep = "",
    file = paste0(tabPath, "did.tex"))

# Raw Vals
cat("\\textbf{Panel (A): Other Weapon Crime:} &&&&\\\\ \\\\",
    # "Estimated ATT: &", DIDcoef(did_main_results$did_rob_coefs[,1]), 
    # "&", DIDcoef(did_main_results$did_ass_coefs[,1]),
    # "&", DIDcoef(did_main_results$did_lnrob_coefs[,1]),
    # "&", DIDcoef(did_main_results$did_lnass_coefs[,1]),
    # "\\\\",
    # "&", DIDci(did_main_results$did_rob_coefs[,1]), 
    # "&", DIDci(did_main_results$did_ass_coefs[,1]),
    # "&", DIDci(did_main_results$did_lnrob_coefs[,1]),
    # "&", DIDci(did_main_results$did_lnass_coefs[,1]),
    "Estimated ATT: &", DIDcoef(did_main_results$did_rob_o_coefs[,1]), 
    "&", DIDcoef(did_main_results$did_ass_o_coefs[,1]),
    "&", DIDcoef(did_main_results$did_lnrob_o_coefs[,1]),
    "&", DIDcoef(did_main_results$did_lnass_o_coefs[,1]),
    "\\\\",
    "&", DIDci(did_main_results$did_rob_o_coefs[,1]), 
    "&", DIDci(did_main_results$did_ass_o_coefs[,1]),
    "&", DIDci(did_main_results$did_lnrob_o_coefs[,1]),
    "&", DIDci(did_main_results$did_lnass_o_coefs[,1]),
    "\\\\ \\\\",    
    "Pre-Legalization Mean ",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_dat$actual_robbery_other_weapon[l], ucr_dat$population[l])), "}",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_dat$actual_assault_other_weapon[l], ucr_dat$population[l])), "}",
    "& \\multicolumn{1}{c}{-- }",
    "& \\multicolumn{1}{c}{-- }",
    "\\\\ \\\\",
    "Implied Percent Change ",
    "& \\multicolumn{1}{c}{", f(100*did_main_results$did_rob_o_coefs[1,1]/weighted.mean(ucr_dat$actual_robbery_other_weapon[l], ucr_dat$population[l])), "\\%}",
    "& \\multicolumn{1}{c}{", f(100*did_main_results$did_ass_o_coefs[1,1]/weighted.mean(ucr_dat$actual_assault_other_weapon[l], ucr_dat$population[l])), "\\%}",
    "& \\multicolumn{1}{c}{ --}",
    "& \\multicolumn{1}{c}{--}",
    # "\\\\ ",
    # "&", DIDci_per(did_main_results$did_rob_o_coefs[,1], m = weighted.mean(ucr_dat$actual_robbery_with_a_gun[l], ucr_dat$population[l])),
    # "&", DIDci_per(did_main_results$did_ass_o_coefs[,1], m = weighted.mean(ucr_dat$actual_assault_with_a_gun[l], ucr_dat$population[l])),
    # "&", DIDci_per(did_main_results$did_lnrob_o_coefs[,1], m = weighted.mean(sui_dat$num_homi[k], sui_dat$total_population[k])),
    # "&", DIDci_per(did_main_results$did_lnass_o_coefs[,1], m = weighted.mean(sui_dat$num_sui[k], sui_dat$total_population[k])),
    "\\\\ \\\\",
    "Average Pre-trend Violation: &", f(weighted.mean(did_main_results$did_rob_o_coefs[1,2:5], 1/colSds(did_main_results$did_rob_o_coefs[,2:5]))), 
    "&", f(weighted.mean(did_main_results$did_ass_o_coefs[1,2:5], 1/colSds(did_main_results$did_ass_o_coefs[,2:5]))),
    "&", f(weighted.mean(did_main_results$did_lnrob_o_coefs[1,2:5], 1/colSds(did_main_results$did_lnrob_o_coefs[,2:5]))),
    "&", f(weighted.mean(did_main_results$did_lnass_o_coefs[1,2:5], 1/colSds(did_main_results$did_lnass_o_coefs[,2:5]))),
    "\\\\ \\\\",    
    # "SD of Pre-trend Violations: &", f(sqrt(wtd.var(did_main_results$did_rob_o_coefs[1,2:5], 1/colSds(did_main_results$did_rob_o_coefs[,2:5]), normwt = TRUE))), 
    # "&", f(sqrt(wtd.var(did_main_results$did_ass_o_coefs[1,2:5], 1/colSds(did_main_results$did_ass_o_coefs[,2:5]), normwt = TRUE))),
    # "&", f(sqrt(wtd.var(did_main_results$did_lnrob_o_coefs[1,2:5], 1/colSds(did_main_results$did_lnrob_o_coefs[,2:5]), normwt = TRUE))),
    # "&", f(sqrt(wtd.var(did_main_results$did_lnass_o_coefs[1,2:5], 1/colSds(did_main_results$did_lnass_o_coefs[,2:5]), normwt = TRUE))),
    # "\\\\ \\\\",
    "Slope of Pre-trend Violations: &", pre_trend_slope(did_main_results$did_rob_o_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_ass_o_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_lnrob_o_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_lnass_o_coefs[,2:5]), 
    "\\\\",
    "&", pre_trend_slope_ci(did_main_results$did_rob_o_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_ass_o_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_lnrob_o_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_lnass_o_coefs[,2:5]),
    "\\\\ \\\\",
    "Largest Pre-trend Violation: &", DIDcoef(did_main_results$did_rob_o_coefs[,(max_pt[1]+2)]), 
    "&", DIDcoef(did_main_results$did_ass_o_coefs[,(max_pt[2]+2)]),
    "&", DIDcoef(did_main_results$did_lnrob_o_coefs[,(max_pt[3]+2)]),
    "&", DIDcoef(did_main_results$did_lnass_o_coefs[,(max_pt[4]+2)]),
    "\\\\",
    "&", DIDci(did_main_results$did_rob_o_coefs[,(max_pt[1]+2)]), 
    "&", DIDci(did_main_results$did_ass_o_coefs[,(max_pt[2]+2)]),
    "&", DIDci(did_main_results$did_lnrob_o_coefs[,(max_pt[3]+2)]),
    "&", DIDci(did_main_results$did_lnass_o_coefs[,(max_pt[4]+2)]),
    "\\\\ \\\\ ",
    
    
    "\\textbf{Panel (B): Unarmed Crime:} &&&&\\\\ \\\\",
    "Estimated ATT: &", DIDcoef(did_main_results$did_rob_u_coefs[,1]), 
    "&", DIDcoef(did_main_results$did_ass_u_coefs[,1]),
    "&", DIDcoef(did_main_results$did_lnrob_u_coefs[,1]),
    "&", DIDcoef(did_main_results$did_lnass_u_coefs[,1]),
    "\\\\",
    "&", DIDci(did_main_results$did_rob_u_coefs[,1]), 
    "&", DIDci(did_main_results$did_ass_u_coefs[,1]),
    "&", DIDci(did_main_results$did_lnrob_u_coefs[,1]),
    "&", DIDci(did_main_results$did_lnass_u_coefs[,1]),
    "\\\\ \\\\",  
    "Pre-Legalization Mean ",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_dat$actual_robbery_unarmed[l], ucr_dat$population[l])), "}",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_dat$actual_assault_unarmed[l], ucr_dat$population[l])), "}",
    "& \\multicolumn{1}{c}{--}",
    "& \\multicolumn{1}{c}{--}",
    "\\\\ \\\\",
    "Implied Percent Change ",
    "& \\multicolumn{1}{c}{", f(100*did_main_results$did_rob_u_coefs[1,1]/weighted.mean(ucr_dat$actual_robbery_unarmed[l], ucr_dat$population[l])), "\\%}",
    "& \\multicolumn{1}{c}{", f(100*did_main_results$did_ass_u_coefs[1,1]/weighted.mean(ucr_dat$actual_assault_unarmed[l], ucr_dat$population[l])), "\\%}",
    "& \\multicolumn{1}{c}{--}",
    "& \\multicolumn{1}{c}{--}",
    # "\\\\ ",
    # "&", DIDci_per(did_main_results$did_rob_u_coefs[,1], m = weighted.mean(ucr_dat$actual_robbery_with_a_knife[l], ucr_dat$population[l])), 
    # "&", DIDci_per(did_main_results$did_ass_u_coefs[,1], m = weighted.mean(ucr_dat$actual_assault_with_a_knife[l], ucr_dat$population[l])),
    # "&", DIDci_per(did_main_results$did_lnrob_u_coefs[,1], m = weighted.mean(sui_dat_o$num_homi[k], sui_dat$total_population[k])),
    # "&", DIDci_per(did_main_results$did_lnass_u_coefs[,1], m = weighted.mean(sui_dat_o$num_sui[k], sui_dat$total_population[k])),
    "\\\\ \\\\",
    "Average Pre-trend Violation: &", f(weighted.mean(did_main_results$did_rob_u_coefs[1,2:5], 1/colSds(did_main_results$did_rob_u_coefs[,2:5]))), 
    "&", f(weighted.mean(did_main_results$did_ass_u_coefs[1,2:5], 1/colSds(did_main_results$did_ass_u_coefs[,2:5]))),
    "&", f(weighted.mean(did_main_results$did_lnrob_u_coefs[1,2:5], 1/colSds(did_main_results$did_lnrob_u_coefs[,2:5]))),
    "&", f(weighted.mean(did_main_results$did_lnass_u_coefs[1,2:5], 1/colSds(did_main_results$did_lnass_u_coefs[,2:5]))),
    "\\\\ \\\\",    
    # "SD of Pre-trend Violations: &", f(sqrt(wtd.var(did_main_results$did_rob_u_coefs[1,2:5], 1/colSds(did_main_results$did_rob_u_coefs[,2:5]), normwt = TRUE))), 
    # "&", f(sqrt(wtd.var(did_main_results$did_ass_u_coefs[1,2:5], 1/colSds(did_main_results$did_ass_u_coefs[,2:5]), normwt = TRUE))),
    # "&", f(sqrt(wtd.var(did_main_results$did_lnrob_u_coefs[1,2:5], 1/colSds(did_main_results$did_lnrob_u_coefs[,2:5]), normwt = TRUE))),
    # "&", f(sqrt(wtd.var(did_main_results$did_lnass_u_coefs[1,2:5], 1/colSds(did_main_results$did_lnass_u_coefs[,2:5]), normwt = TRUE))),
    # "\\\\ \\\\",
    "Slope of Pre-trend Violations: &", pre_trend_slope(did_main_results$did_rob_u_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_ass_u_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_lnrob_u_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did_lnass_u_coefs[,2:5]), 
    "\\\\",
    "&", pre_trend_slope_ci(did_main_results$did_rob_u_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_ass_u_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_lnrob_u_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did_lnass_u_coefs[,2:5]),
    "\\\\ \\\\",
    "Largest Pre-trend Violation: &", DIDcoef(did_main_results$did_rob_u_coefs[,(max_pt[1]+2)]), 
    "&", DIDcoef(did_main_results$did_ass_u_coefs[,(max_pt[2]+2)]),
    "&", DIDcoef(did_main_results$did_lnrob_u_coefs[,(max_pt[3]+2)]),
    "&", DIDcoef(did_main_results$did_lnass_u_coefs[,(max_pt[4]+2)]),
    "\\\\",
    "&", DIDci(did_main_results$did_rob_u_coefs[,(max_pt[1]+2)]), 
    "&", DIDci(did_main_results$did_ass_u_coefs[,(max_pt[2]+2)]),
    "&", DIDci(did_main_results$did_lnrob_u_coefs[,(max_pt[3]+2)]),
    "&", DIDci(did_main_results$did_lnass_u_coefs[,(max_pt[4]+2)]),
    "\\\\ \\\\",
    
    " \\hline ",
    "\\\\[0.5mm]",
    "Gun Law Controls? & Yes & Yes & Yes & Yes",
    "\\\\[0.5mm]",
    "Unit of Observation ",
    "& \\multicolumn{4}{c}{ Jurisdiction - Year}",
    "\\\\[0.5mm]",
    "Number of Observations ",
    "& \\multicolumn{4}{c}{", f(nrow(ucr_dat)), "}",
    "\\\\[0.5mm]", 
    "\\bottomrule",
    sep = "",
    file = paste0(tabPath, "did_2.tex"))



# # Get max pre-trend
# max_pt <- c(
#   which.max(abs(did_main_results$did_rob_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did_ass_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did_lnrob_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did_lnass_coefs[1, 3:5])))

# l <- which(ucr_dat$trip_year_since_treat == "-1")
# k <- which(sui_dat_stack$trip_year_since_treat == "-1")
# 
# # Get max pre-trend
# max_pt <- c(
#   which.max(abs(did_main_results$did_rob_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did_ass_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did_lnrob_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did_lnass_coefs[1, 3:5])))
# 
# # Raw Vals
# cat("Estimated ATT: &", DIDcoef(did_main_results$did_rob_coefs[,1]), 
#     "&", DIDcoef(did_main_results$did3_ass_coefs[,1]),
#     "&", DIDcoef(did_main_results$did3_lnrob_coefs[,1]),
#     "&", DIDcoef(did_main_results$did3_lnass_coefs[,1]),
#     "\\\\",
#     "&", DIDci(did_main_results$did3_rob_coefs[,1]), 
#     "&", DIDci(did_main_results$did3_ass_coefs[,1]),
#     "&", DIDci(did_main_results$did3_lnrob_coefs[,1]),
#     "&", DIDci(did_main_results$did3_lnass_coefs[,1]),
#     "\\\\ \\\\",    
#     "Average Pre-trend Violations: &", f(weighted.mean(did_main_results$did3_rob_coefs[1,2:5], 1/colSds(did_main_results$did3_rob_coefs[,2:5]))), 
#     "&", f(weighted.mean(did_main_results$did3_ass_coefs[1,2:5], 1/colSds(did_main_results$did3_ass_coefs[,2:5]))),
#     "&", f(weighted.mean(did_main_results$did3_lnrob_coefs[1,2:5], 1/colSds(did_main_results$did3_lnrob_coefs[,2:5]))),
#     "&", f(weighted.mean(did_main_results$did3_lnass_coefs[1,2:5], 1/colSds(did_main_results$did3_lnass_coefs[,2:5]))),
#     "\\\\ \\\\",    
#     "SD of Pre-trend Violations: &", f(sqrt(wtd.var(did_main_results$did_rob_coefs[1,2:5], 1/colSds(did_main_results$did3_rob_coefs[,2:5]), normwt = TRUE))), 
#     "&", f(sqrt(wtd.var(did_main_results$did3_ass_coefs[1,2:5], 1/colSds(did_main_results$did3_ass_coefs[,2:5]), normwt = TRUE))),
#     "&", f(sqrt(wtd.var(did_main_results$did3_lnrob_coefs[1,2:5], 1/colSds(did_main_results$did3_lnrob_coefs[,2:5]), normwt = TRUE))),
#     "&", f(sqrt(wtd.var(did_main_results$did3_lnass_coefs[1,2:5], 1/colSds(did_main_results$did3_lnass_coefs[,2:5]), normwt = TRUE))),
#     "\\\\ \\\\",
#     "Largest Pre-trend Violation: &", DIDcoef(did_main_results$did3_rob_coefs[,(max_pt[1]+2)]), 
#     "&", DIDcoef(did_main_results$did3_ass_coefs[,(max_pt[2]+2)]),
#     "&", DIDcoef(did_main_results$did3_lnrob_coefs[,(max_pt[3]+2)]),
#     "&", DIDcoef(did_main_results$did3_lnass_coefs[,(max_pt[4]+2)]),
#     "\\\\",
#     "&", DIDci(did_main_results$did3_rob_coefs[,(max_pt[1]+2)]), 
#     "&", DIDci(did_main_results$did3_ass_coefs[,(max_pt[2]+2)]),
#     "&", DIDci(did_main_results$did3_lnrob_coefs[,(max_pt[3]+2)]),
#     "&", DIDci(did_main_results$did3_lnass_coefs[,(max_pt[4]+2)]),
#     # "\\\\ \\\\",
#     # "Non-Parallel Trends Interval: &",  DIDcred(did_main_results$did3_rob_coefs[,-c(1,2)]),
#     # "&", DIDcred(did_main_results$did3_ass_coefs[,-c(1,2)]),
#     # "&", DIDcred(did_main_results$did3_lnrob_coefs[,-c(1,2)]),
#     # "&", DIDcred(did_main_results$did3_lnass_coefs[,-c(1,2)]),
#     "\\\\ ",
#     " \\hline ",
#     "\\\\[0.5mm]",
#     "Pre-Legalization Mean ",
#     "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_dat$robbery[l], ucr_dat$population[l])), "}",
#     "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_dat$assaults[l], ucr_dat$population[l])), "}",
#     "& \\multicolumn{1}{c}{", f(weighted.mean(sui_dat_stack$num_homi[k], sui_dat_stack$total_population[k])), "}",
#     "& \\multicolumn{1}{c}{", f(weighted.mean(sui_dat_stack$num_sui[k], sui_dat_stack$total_population[k])), "}",
#     "\\\\[0.5mm]",
#     # "Implied Percent Change ",
#     # "& \\multicolumn{1}{c}{", f(DIDcoef(did_main_results$did3_rob_coefs[,1])/weighted.mean(ucr_dat$robbery[l], ucr_dat$population[l])), "}",
#     # "& \\multicolumn{1}{c}{", f(DIDcoef(did_main_results$did3_ass_coefs[,1])/weighted.mean(ucr_dat$assaults[l], ucr_dat$population[l])), "}",
#     # "& \\multicolumn{1}{c}{", f(DIDcoef(did_main_results$did3_lnrob_coefs[,1])/weighted.mean(sui_dat_stack$num_homi[k], sui_dat_stack$total_population[k])), "}",
#     # "& \\multicolumn{1}{c}{", f(DIDcoef(did_main_results$did3_lnass_coefs[,1])/weighted.mean(sui_dat_stack$num_sui[k], sui_dat_stack$total_population[k])), "}",
#     # "\\\\[0.5mm]",
#     "Unit of Observation ",
#     "& \\multicolumn{2}{c}{ Jurisdiction}",
#     "& \\multicolumn{2}{c}{ State}",
#     "\\\\[0.5mm]",
#     "Number of Observations ",
#     "& \\multicolumn{2}{c}{", f(nrow(ucr_dat)), "}",
#     "& \\multicolumn{2}{c}{", f(nrow(sui_dat_stack)), "}",
#     "\\\\[0.5mm]",
#     "\\bottomrule",
#     sep = "",
#     file = paste0(tabPath, "did_3d_logs.tex"))


# Triple Ds ---------------------------------------------------------------
# ucr_stack <- ucr_stack %>% 
#   group_by(type_ori) %>% 
#   summarize(ever_treat = max(trip_treat)) %>% 
#   right_join(ucr_stack)

l <- which(ucr_stack$trip_year_since_treat == "-1g")
k <- which(sui_dat_stack$trip_year_since_treat == "-1g")

# Get max pre-trend
max_pt <- c(
  which.max(abs(did_main_results$did3_rob_coefs[1, 3:5])),
  which.max(abs(did_main_results$did3_ass_coefs[1, 3:5])),
  which.max(abs(did_main_results$did3_lnrob_coefs[1, 3:5])),
  which.max(abs(did_main_results$did3_lnass_coefs[1, 3:5])))

# Raw Vals
cat("Estimated ATT: &", DIDcoef(did_main_results$did3_rob_coefs[,1]), 
    "&", DIDcoef(did_main_results$did3_ass_coefs[,1]),
    "&", DIDcoef(did_main_results$did3_lnrob_coefs[,1]),
    "&", DIDcoef(did_main_results$did3_lnass_coefs[,1]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,1]), 
    "&", DIDci(did_main_results$did3_ass_coefs[,1]),
    "&", DIDci(did_main_results$did3_lnrob_coefs[,1]),
    "&", DIDci(did_main_results$did3_lnass_coefs[,1]),
    "\\\\ \\\\",    
    "Pre-Legalization Mean ",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_stack$actual_robbery_with_a_gun[l], ucr_stack$population[l])), "}",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_stack$actual_assault_with_a_gun[l], ucr_stack$population[l])), "}",
    "& \\multicolumn{1}{c}{--}",
    "& \\multicolumn{1}{c}{--}",
    "\\\\ \\\\",
    "Implied Percent Change ",
    "& \\multicolumn{1}{c}{", f(100*did_main_results$did3_rob_coefs[1,1]/weighted.mean(ucr_stack$actual_robbery_with_a_gun[l], ucr_stack$population[l])), "\\%}",
    "& \\multicolumn{1}{c}{", f(100*did_main_results$did3_ass_coefs[1,1]/weighted.mean(ucr_stack$actual_assault_with_a_gun[l], ucr_stack$population[l])), "\\%}",
    "& \\multicolumn{1}{c}{--}",
    "& \\multicolumn{1}{c}{--}",
    # "\\\\",
    # "&", DIDci_per(did_main_results$did3_rob_coefs[,1], m = weighted.mean(ucr_stack$robbery[l], ucr_stack$population[l])),
    # "&", DIDci_per(did_main_results$did3_ass_coefs[,1], m = weighted.mean(ucr_stack$assaults[l], ucr_stack$population[l])),
    # "&", DIDci_per(did_main_results$did3_lnrob_coefs[,1], m = weighted.mean(sui_dat_stack$num_homi[k], sui_dat_stack$total_population[k])),
    # "&", DIDci_per(did_main_results$did3_lnass_coefs[,1], m = weighted.mean(sui_dat_stack$num_sui[k], sui_dat_stack$total_population[k])),
    "\\\\ \\\\",
    "Average Pre-trend Violation: &", f(weighted.mean(did_main_results$did3_rob_coefs[1,2:5], 1/colSds(did_main_results$did3_rob_coefs[,2:5]))), 
    "&", f(weighted.mean(did_main_results$did3_ass_coefs[1,2:5], 1/colSds(did_main_results$did3_ass_coefs[,2:5]))),
    "&", f(weighted.mean(did_main_results$did3_lnrob_coefs[1,2:5], 1/colSds(did_main_results$did3_lnrob_coefs[,2:5]))),
    "&", f(weighted.mean(did_main_results$did3_lnass_coefs[1,2:5], 1/colSds(did_main_results$did3_lnass_coefs[,2:5]))),
    "\\\\ \\\\",    
    # "SD of Pre-trend Violations: &", f(sqrt(wtd.var(did_main_results$did3_rob_coefs[1,2:5], 1/colSds(did_main_results$did3_rob_coefs[,2:5]), normwt = TRUE))), 
    # "&", f(sqrt(wtd.var(did_main_results$did3_ass_coefs[1,2:5], 1/colSds(did_main_results$did3_ass_coefs[,2:5]), normwt = TRUE))),
    # "&", f(sqrt(wtd.var(did_main_results$did3_lnrob_coefs[1,2:5], 1/colSds(did_main_results$did3_lnrob_coefs[,2:5]), normwt = TRUE))),
    # "&", f(sqrt(wtd.var(did_main_results$did3_lnass_coefs[1,2:5], 1/colSds(did_main_results$did3_lnass_coefs[,2:5]), normwt = TRUE))),
    # "\\\\ \\\\",
    "Slope of Pre-trend Violations: &", pre_trend_slope(did_main_results$did3_rob_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did3_ass_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did3_lnrob_coefs[,2:5]), 
    "&", pre_trend_slope(did_main_results$did3_lnass_coefs[,2:5]), 
    "\\\\",
    "&", pre_trend_slope_ci(did_main_results$did3_rob_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did3_ass_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did3_lnrob_coefs[,2:5]), 
    "&", pre_trend_slope_ci(did_main_results$did3_lnass_coefs[,2:5]),
    "\\\\ \\\\",
    "Largest Pre-trend Violation: &", DIDcoef(did_main_results$did3_rob_coefs[,(max_pt[1]+2)]), 
    "&", DIDcoef(did_main_results$did3_ass_coefs[,(max_pt[2]+2)]),
    "&", DIDcoef(did_main_results$did3_lnrob_coefs[,(max_pt[3]+2)]),
    "&", DIDcoef(did_main_results$did3_lnass_coefs[,(max_pt[4]+2)]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,(max_pt[1]+2)]), 
    "&", DIDci(did_main_results$did3_ass_coefs[,(max_pt[2]+2)]),
    "&", DIDci(did_main_results$did3_lnrob_coefs[,(max_pt[3]+2)]),
    "&", DIDci(did_main_results$did3_lnass_coefs[,(max_pt[4]+2)]),
    # "\\\\ \\\\",
    # "Non-Parallel Trends Interval: &",  DIDcred(did_main_results$did3_rob_coefs[,-c(1,2)]),
    # "&", DIDcred(did_main_results$did3_ass_coefs[,-c(1,2)]),
    # "&", DIDcred(did_main_results$did3_lnrob_coefs[,-c(1,2)]),
    # "&", DIDcred(did_main_results$did3_lnass_coefs[,-c(1,2)]),
    "\\\\",
    " \\hline ",
    "\\\\",
    "Gun Law Controls? & No & No & No & No",
    "\\\\[0.5mm]",  
    "Non-Gun Crime Type? & All & All & All & All",
    "\\\\[0.5mm]",    
    "Unit of Observation ",
    "& \\multicolumn{4}{c}{ Jurisdiction - Year - Type of Offense}",
    "\\\\[0.5mm]",
    "Number of Observations ",
    "& \\multicolumn{4}{c}{", f(nrow(ucr_stack)), "}",
    "\\\\[0.5mm]",
    "\\bottomrule",
    sep = "",
    file = paste0(tabPath, "did_3d.tex"))

# Logs

# # Get max pre-trend
# max_pt <- c(
#   which.max(abs(did_main_results$did3_rob_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did3_ass_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did3_lnrob_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did3_lnass_coefs[1, 3:5])))
# 
# l <- which(ucr_stack$trip_year_since_treat == "-1g")
# k <- which(sui_dat_stack$trip_year_since_treat == "-1g")
# 
# # Get max pre-trend
# max_pt <- c(
#   which.max(abs(did_main_results$did3_rob_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did3_ass_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did3_lnrob_coefs[1, 3:5])),
#   which.max(abs(did_main_results$did3_lnass_coefs[1, 3:5])))
# 
# # Raw Vals
# cat("Estimated ATT: &", DIDcoef(did_main_results$did3_rob_coefs[,1]), 
#     "&", DIDcoef(did_main_results$did3_ass_coefs[,1]),
#     "&", DIDcoef(did_main_results$did3_lnrob_coefs[,1]),
#     "&", DIDcoef(did_main_results$did3_lnass_coefs[,1]),
#     "\\\\",
#     "&", DIDci(did_main_results$did3_rob_coefs[,1]), 
#     "&", DIDci(did_main_results$did3_ass_coefs[,1]),
#     "&", DIDci(did_main_results$did3_lnrob_coefs[,1]),
#     "&", DIDci(did_main_results$did3_lnass_coefs[,1]),
#     "\\\\ \\\\",    
#     "Average Pre-trend Violations: &", f(weighted.mean(did_main_results$did3_rob_coefs[1,2:5], 1/colSds(did_main_results$did3_rob_coefs[,2:5]))), 
#     "&", f(weighted.mean(did_main_results$did3_ass_coefs[1,2:5], 1/colSds(did_main_results$did3_ass_coefs[,2:5]))),
#     "&", f(weighted.mean(did_main_results$did3_lnrob_coefs[1,2:5], 1/colSds(did_main_results$did3_lnrob_coefs[,2:5]))),
#     "&", f(weighted.mean(did_main_results$did3_lnass_coefs[1,2:5], 1/colSds(did_main_results$did3_lnass_coefs[,2:5]))),
#     "\\\\ \\\\",    
#     "SD of Pre-trend Violations: &", f(sqrt(wtd.var(did_main_results$did_rob_coefs[1,2:5], 1/colSds(did_main_results$did3_rob_coefs[,2:5]), normwt = TRUE))), 
#     "&", f(sqrt(wtd.var(did_main_results$did3_ass_coefs[1,2:5], 1/colSds(did_main_results$did3_ass_coefs[,2:5]), normwt = TRUE))),
#     "&", f(sqrt(wtd.var(did_main_results$did3_lnrob_coefs[1,2:5], 1/colSds(did_main_results$did3_lnrob_coefs[,2:5]), normwt = TRUE))),
#     "&", f(sqrt(wtd.var(did_main_results$did3_lnass_coefs[1,2:5], 1/colSds(did_main_results$did3_lnass_coefs[,2:5]), normwt = TRUE))),
#     "\\\\ \\\\",
#     "Largest Pre-trend Violation: &", DIDcoef(did_main_results$did3_rob_coefs[,(max_pt[1]+2)]), 
#     "&", DIDcoef(did_main_results$did3_ass_coefs[,(max_pt[2]+2)]),
#     "&", DIDcoef(did_main_results$did3_lnrob_coefs[,(max_pt[3]+2)]),
#     "&", DIDcoef(did_main_results$did3_lnass_coefs[,(max_pt[4]+2)]),
#     "\\\\",
#     "&", DIDci(did_main_results$did3_rob_coefs[,(max_pt[1]+2)]), 
#     "&", DIDci(did_main_results$did3_ass_coefs[,(max_pt[2]+2)]),
#     "&", DIDci(did_main_results$did3_lnrob_coefs[,(max_pt[3]+2)]),
#     "&", DIDci(did_main_results$did3_lnass_coefs[,(max_pt[4]+2)]),
#     # "\\\\ \\\\",
#     # "Non-Parallel Trends Interval: &",  DIDcred(did_main_results$did3_rob_coefs[,-c(1,2)]),
#     # "&", DIDcred(did_main_results$did3_ass_coefs[,-c(1,2)]),
#     # "&", DIDcred(did_main_results$did3_lnrob_coefs[,-c(1,2)]),
#     # "&", DIDcred(did_main_results$did3_lnass_coefs[,-c(1,2)]),
#     "\\\\ ",
#     " \\hline ",
#     "\\\\[0.5mm]",
#     "Pre-Legalization Mean ",
#     "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_stack$robbery[l], ucr_stack$population[l])), "}",
#     "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_stack$assaults[l], ucr_stack$population[l])), "}",
#     "& \\multicolumn{1}{c}{", f(weighted.mean(sui_dat_stack$num_homi[k], sui_dat_stack$total_population[k])), "}",
#     "& \\multicolumn{1}{c}{", f(weighted.mean(sui_dat_stack$num_sui[k], sui_dat_stack$total_population[k])), "}",
#     "\\\\[0.5mm]",
#     # "Implied Percent Change ",
#     # "& \\multicolumn{1}{c}{", f(DIDcoef(did_main_results$did3_rob_coefs[,1])/weighted.mean(ucr_stack$robbery[l], ucr_stack$population[l])), "}",
#     # "& \\multicolumn{1}{c}{", f(DIDcoef(did_main_results$did3_ass_coefs[,1])/weighted.mean(ucr_stack$assaults[l], ucr_stack$population[l])), "}",
#     # "& \\multicolumn{1}{c}{", f(DIDcoef(did_main_results$did3_lnrob_coefs[,1])/weighted.mean(sui_dat_stack$num_homi[k], sui_dat_stack$total_population[k])), "}",
#     # "& \\multicolumn{1}{c}{", f(DIDcoef(did_main_results$did3_lnass_coefs[,1])/weighted.mean(sui_dat_stack$num_sui[k], sui_dat_stack$total_population[k])), "}",
#     # "\\\\[0.5mm]",
#     "Unit of Observation ",
#     "& \\multicolumn{2}{c}{ Jurisdiction}",
#     "& \\multicolumn{2}{c}{ State}",
#     "\\\\[0.5mm]",
#     "Number of Observations ",
#     "& \\multicolumn{2}{c}{", f(nrow(ucr_stack)), "}",
#     "& \\multicolumn{2}{c}{", f(nrow(sui_dat_stack)), "}",
#     "\\\\[0.5mm]",
#     "\\bottomrule",
#     sep = "",
#     file = paste0(tabPath, "did_3d_logs.tex"))


# Base Diffs Event --------------------------------------------------------------
# # Ever Treat variable
# ucr_dat <- ucr_dat %>% 
#   group_by(state) %>% 
#   summarize(ever_treat = max(treat)) %>% 
#   right_join(ucr_dat)
# 
# l <- which(ucr_dat$treat == 0 &
#              ucr_dat$ever_treat == 1)
# 
# cat("Treat -5+&", DIDcoef(did_main_results$did_rob_coefs[,2]), 
#     "&", DIDcoef(did_main_results$did_rob_coefs[,2]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,2]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,2]),
#     "\\\\",
#     "&", DIDci(did_main_results$did_rob_coefs[,2]), 
#     "&", DIDci(did_main_results$did_rob_coefs[,2]),
#     "&", DIDci(did_main_results$did_ass_coefs[,2]),
#     "&", DIDci(did_main_results$did_ass_coefs[,2]),
#     "\\\\ \\\\", 
#     "Treat -4 &", DIDcoef(did_main_results$did_rob_coefs[,3]), 
#     "&", DIDcoef(did_main_results$did_rob_coefs[,3]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,3]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,3]),
#     "\\\\",
#     "&", DIDci(did_main_results$did_rob_coefs[,3]), 
#     "&", DIDci(did_main_results$did_rob_coefs[,3]),
#     "&", DIDci(did_main_results$did_ass_coefs[,3]),
#     "&", DIDci(did_main_results$did_ass_coefs[,3]),
#     "\\\\ \\\\",
#     "Treat -3 &", DIDcoef(did_main_results$did_rob_coefs[,4]), 
#     "&", DIDcoef(did_main_results$did_rob_coefs[,4]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,4]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,4]),
#     "\\\\",
#     "&", DIDci(did_main_results$did_rob_coefs[,4]), 
#     "&", DIDci(did_main_results$did_rob_coefs[,4]),
#     "&", DIDci(did_main_results$did_ass_coefs[,4]),
#     "&", DIDci(did_main_results$did_ass_coefs[,4]),
#     "\\\\ \\\\",
#     "Treat -2 &", DIDcoef(did_main_results$did_rob_coefs[,5]), 
#     "&", DIDcoef(did_main_results$did_rob_coefs[,5]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,5]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,5]),
#     "\\\\",
#     "&", DIDci(did_main_results$did_rob_coefs[,5]), 
#     "&", DIDci(did_main_results$did_rob_coefs[,5]),
#     "&", DIDci(did_main_results$did_ass_coefs[,5]),
#     "&", DIDci(did_main_results$did_ass_coefs[,5]),
#     "\\\\ \\\\",
#     "Treat 0 &", DIDcoef(did_main_results$did_rob_coefs[,6]), 
#     "&", DIDcoef(did_main_results$did_rob_coefs[,6]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,6]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,6]),
#     "\\\\",
#     "&", DIDci(did_main_results$did_rob_coefs[,6]), 
#     "&", DIDci(did_main_results$did_rob_coefs[,6]),
#     "&", DIDci(did_main_results$did_ass_coefs[,6]),
#     "&", DIDci(did_main_results$did_ass_coefs[,6]),
#     "\\\\ \\\\",
#     "Treat +1 &", DIDcoef(did_main_results$did_rob_coefs[,7]), 
#     "&", DIDcoef(did_main_results$did_rob_coefs[,7]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,7]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,7]),
#     "\\\\",
#     "&", DIDci(did_main_results$did_rob_coefs[,7]), 
#     "&", DIDci(did_main_results$did_rob_coefs[,7]),
#     "&", DIDci(did_main_results$did_ass_coefs[,7]),
#     "&", DIDci(did_main_results$did_ass_coefs[,7]),
#     "\\\\ \\\\",
#     "Treat +2 &", DIDcoef(did_main_results$did_rob_coefs[,8]), 
#     "&", DIDcoef(did_main_results$did_rob_coefs[,8]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,8]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,8]),
#     "\\\\",
#     "&", DIDci(did_main_results$did_rob_coefs[,8]), 
#     "&", DIDci(did_main_results$did_rob_coefs[,8]),
#     "&", DIDci(did_main_results$did_ass_coefs[,8]),
#     "&", DIDci(did_main_results$did_ass_coefs[,8]),
#     "\\\\ \\\\",
#     "Treat +3 &", DIDcoef(did_main_results$did_rob_coefs[,9]), 
#     "&", DIDcoef(did_main_results$did_rob_coefs[,9]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,9]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,9]),
#     "\\\\",
#     "&", DIDci(did_main_results$did_rob_coefs[,9]), 
#     "&", DIDci(did_main_results$did_rob_coefs[,9]),
#     "&", DIDci(did_main_results$did_ass_coefs[,9]),
#     "&", DIDci(did_main_results$did_ass_coefs[,9]),
#     "\\\\ \\\\",
#     "Treat +4 &", DIDcoef(did_main_results$did_rob_coefs[,10]), 
#     "&", DIDcoef(did_main_results$did_rob_coefs[,10]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,10]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,10]),
#     "\\\\",
#     "&", DIDci(did_main_results$did_rob_coefs[,10]), 
#     "&", DIDci(did_main_results$did_rob_coefs[,10]),
#     "&", DIDci(did_main_results$did_ass_coefs[,10]),
#     "&", DIDci(did_main_results$did_ass_coefs[,10]),
#     "\\\\ \\\\",
#     "Treat +5+ &", DIDcoef(did_main_results$did_rob_coefs[,11]),
#     "&", DIDcoef(did_main_results$did_rob_coefs[,11]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,11]),
#     "&", DIDcoef(did_main_results$did_ass_coefs[,11]),
#     "\\\\",
#     "&", DIDci(did_main_results$did_rob_coefs[,11]),
#     "&", DIDci(did_main_results$did_rob_coefs[,11]),
#     "&", DIDci(did_main_results$did_ass_coefs[,11]),
#     "&", DIDci(did_main_results$did_ass_coefs[,11]),
#     "\\\\",
#     " \\hline ",
#     "\\\\[0.5mm]",
#     "Pre-Legalization Mean & \\multicolumn{1}{c}{", f(mean(ucr_dat$actual_robbery_with_a_gun[l])), "}",
#     "&& \\multicolumn{1}{c}{", f(mean(ucr_dat$actual_assault_with_a_gun[l])), "} &\\\\[0.5mm]
#     \\bottomrule",
#     sep = "",
#     file = paste0(tabPath, "did_base_event.tex"))
# 
# 

# Triple Ds ---------------------------------------------------------------
ucr_stack <- ucr_stack %>% 
  group_by(type_ori) %>% 
  summarise(ever_treat = max(trip_treat)) %>% 
  right_join(ucr_stack)

l <- which(ucr_stack$trip_treat == 0 &
             ucr_stack$ever_treat == 1 )

# Raw Values
cat("Treat -5+&", DIDcoef(did_main_results$did3_rob_coefs[,2]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,2]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,2]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,2]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,2]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,2]),
    "&", DIDci(did_main_results$did3_ass_coefs[,2]),
    "&", DIDci(did_main_results$did3_ass_coefs[,2]),
    "\\\\ \\\\", 
    "Treat -4 &", DIDcoef(did_main_results$did3_rob_coefs[,3]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,3]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,3]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,3]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,3]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,3]),
    "&", DIDci(did_main_results$did3_ass_coefs[,3]),
    "&", DIDci(did_main_results$did3_ass_coefs[,3]),
    "\\\\ \\\\",
    "Treat -3 &", DIDcoef(did_main_results$did3_rob_coefs[,4]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,4]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,4]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,4]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,4]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,4]),
    "&", DIDci(did_main_results$did3_ass_coefs[,4]),
    "&", DIDci(did_main_results$did3_ass_coefs[,4]),
    "\\\\ \\\\",
    "Treat -2 &", DIDcoef(did_main_results$did3_rob_coefs[,5]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,5]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,5]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,5]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,5]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,5]),
    "&", DIDci(did_main_results$did3_ass_coefs[,5]),
    "&", DIDci(did_main_results$did3_ass_coefs[,5]),
    "\\\\ \\\\",
    "Treat 0 &", DIDcoef(did_main_results$did3_rob_coefs[,6]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,6]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,6]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,6]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,6]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,6]),
    "&", DIDci(did_main_results$did3_ass_coefs[,6]),
    "&", DIDci(did_main_results$did3_ass_coefs[,6]),
    "\\\\ \\\\",
    "Treat +1 &", DIDcoef(did_main_results$did3_rob_coefs[,7]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,7]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,7]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,7]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,7]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,7]),
    "&", DIDci(did_main_results$did3_ass_coefs[,7]),
    "&", DIDci(did_main_results$did3_ass_coefs[,7]),
    "\\\\ \\\\",
    "Treat +2 &", DIDcoef(did_main_results$did3_rob_coefs[,8]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,8]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,8]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,8]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,8]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,8]),
    "&", DIDci(did_main_results$did3_ass_coefs[,8]),
    "&", DIDci(did_main_results$did3_ass_coefs[,8]),
    "\\\\ \\\\",
    "Treat +3 &", DIDcoef(did_main_results$did3_rob_coefs[,9]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,9]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,9]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,9]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,9]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,9]),
    "&", DIDci(did_main_results$did3_ass_coefs[,9]),
    "&", DIDci(did_main_results$did3_ass_coefs[,9]),
    "\\\\ \\\\",
    "Treat + 4 &", DIDcoef(did_main_results$did3_rob_coefs[,10]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,10]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,10]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,10]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,10]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,10]),
    "&", DIDci(did_main_results$did3_ass_coefs[,10]),
    "&", DIDci(did_main_results$did3_ass_coefs[,10]),
    "\\\\ \\\\",
    "Treat +5+ &", DIDcoef(did_main_results$did3_rob_coefs[,11]),
    "&", DIDcoef(did_main_results$did3_rob_coefs[,11]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,11]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,11]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,11]),
    "&", DIDci(did_main_results$did3_rob_coefs[,11]),
    "&", DIDci(did_main_results$did3_ass_coefs[,11]),
    "&", DIDci(did_main_results$did3_ass_coefs[,11]),
    "\\\\",
    " \\hline ",
    "\\\\[0.5mm]",
    "Pre-Legalization Mean & \\multicolumn{1}{c}{", f(mean(ucr_stack$robbery[l])), "}",
    "&& \\multicolumn{1}{c}{", f(mean(ucr_stack$assault[l])), "} &\\\\[0.5mm]
    \\bottomrule",
    sep = "",
    file = paste0(tabPath, "did_3d_event.tex"))

# Logs
cat("Treat -5+&", DIDcoef(did_main_results$did3_rob_coefs[,2]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,2]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,2]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,2]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,2]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,2]),
    "&", DIDci(did_main_results$did3_ass_coefs[,2]),
    "&", DIDci(did_main_results$did3_ass_coefs[,2]),
    "\\\\ \\\\", 
    "Treat -4 &", DIDcoef(did_main_results$did3_rob_coefs[,3]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,3]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,3]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,3]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,3]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,3]),
    "&", DIDci(did_main_results$did3_ass_coefs[,3]),
    "&", DIDci(did_main_results$did3_ass_coefs[,3]),
    "\\\\ \\\\",
    "Treat -3 &", DIDcoef(did_main_results$did3_rob_coefs[,4]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,4]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,4]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,4]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,4]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,4]),
    "&", DIDci(did_main_results$did3_ass_coefs[,4]),
    "&", DIDci(did_main_results$did3_ass_coefs[,4]),
    "\\\\ \\\\",
    "Treat -2 &", DIDcoef(did_main_results$did3_rob_coefs[,5]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,5]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,5]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,5]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,5]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,5]),
    "&", DIDci(did_main_results$did3_ass_coefs[,5]),
    "&", DIDci(did_main_results$did3_ass_coefs[,5]),
    "\\\\ \\\\",
    "Treat 0 &", DIDcoef(did_main_results$did3_rob_coefs[,6]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,6]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,6]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,6]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,6]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,6]),
    "&", DIDci(did_main_results$did3_ass_coefs[,6]),
    "&", DIDci(did_main_results$did3_ass_coefs[,6]),
    "\\\\ \\\\",
    "Treat +1 &", DIDcoef(did_main_results$did3_rob_coefs[,7]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,7]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,7]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,7]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,7]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,7]),
    "&", DIDci(did_main_results$did3_ass_coefs[,7]),
    "&", DIDci(did_main_results$did3_ass_coefs[,7]),
    "\\\\ \\\\",
    "Treat +2 &", DIDcoef(did_main_results$did3_rob_coefs[,8]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,8]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,8]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,8]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,8]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,8]),
    "&", DIDci(did_main_results$did3_ass_coefs[,8]),
    "&", DIDci(did_main_results$did3_ass_coefs[,8]),
    "\\\\ \\\\",
    "Treat +3 &", DIDcoef(did_main_results$did3_rob_coefs[,9]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,9]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,9]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,9]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,9]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,9]),
    "&", DIDci(did_main_results$did3_ass_coefs[,9]),
    "&", DIDci(did_main_results$did3_ass_coefs[,9]),
    "\\\\ \\\\",
    "Treat + 4 &", DIDcoef(did_main_results$did3_rob_coefs[,10]), 
    "&", DIDcoef(did_main_results$did3_rob_coefs[,10]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,10]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,10]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,10]), 
    "&", DIDci(did_main_results$did3_rob_coefs[,10]),
    "&", DIDci(did_main_results$did3_ass_coefs[,10]),
    "&", DIDci(did_main_results$did3_ass_coefs[,10]),
    "\\\\ \\\\",
    "Treat +5+ &", DIDcoef(did_main_results$did3_rob_coefs[,11]),
    "&", DIDcoef(did_main_results$did3_rob_coefs[,11]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,11]),
    "&", DIDcoef(did_main_results$did3_ass_coefs[,11]),
    "\\\\",
    "&", DIDci(did_main_results$did3_rob_coefs[,11]),
    "&", DIDci(did_main_results$did3_rob_coefs[,11]),
    "&", DIDci(did_main_results$did3_ass_coefs[,11]),
    "&", DIDci(did_main_results$did3_ass_coefs[,11]),
    "\\\\",
    " \\hline ",
    "\\\\[0.5mm]",
    "Pre-Legalization Mean & \\multicolumn{1}{c}{", f(mean(ucr_stack$robbery[l])), "}",
    "&& \\multicolumn{1}{c}{", f(mean(ucr_stack$assault[l])), "} &\\\\[0.5mm]
    \\bottomrule",
    sep = "",
    file = paste0(tabPath, "did_3d_event_ln.tex"))




# Dispensary --------------------------------------------------------------
load(file = "./Data/did_disp_results.RData")


l <- which(ucr_stack$trip_year_since_treat == "-1g" & ucr_stack$disp25 == 1)
# k <- which(sui_dat_stack$trip_year_since_treat == "-1g")

# Get max pre-trend
max_pt <- c(
  which.max(abs(did_disp_results$did3_rob_coefs[1, 3:5])),
  which.max(abs(did_disp_results$did3_ass_coefs[1, 3:5])),
  which.max(abs(did_disp_results$did3_lnrob_coefs[1, 3:5])),
  which.max(abs(did_disp_results$did3_lnass_coefs[1, 3:5])))

# Raw Vals
cat("Legalization: &", DIDcoef(did_disp_results$did3_rob_coefs[,1]), 
    "&", DIDcoef(did_disp_results$did3_ass_coefs[,1]),
    "&", DIDcoef(did_disp_results$did3_lnrob_coefs[,1]),
    "&", DIDcoef(did_disp_results$did3_lnass_coefs[,1]),
    "\\\\",
    "&", DIDci(did_disp_results$did3_rob_coefs[,1]), 
    "&", DIDci(did_disp_results$did3_ass_coefs[,1]),
    "&", DIDci(did_disp_results$did3_lnrob_coefs[,1]),
    "&", DIDci(did_disp_results$did3_lnass_coefs[,1]),
    "\\\\ \\\\", 
    "Legalization*Dispensary: &", DIDcoef(did_disp_results$did3_rob_coefs[,2]), 
    "&", DIDcoef(did_disp_results$did3_ass_coefs[,2]),
    "&", DIDcoef(did_disp_results$did3_lnrob_coefs[,2]),
    "&", DIDcoef(did_disp_results$did3_lnass_coefs[,2]),
    "\\\\",
    "&", DIDci(did_disp_results$did3_rob_coefs[,2]), 
    "&", DIDci(did_disp_results$did3_ass_coefs[,2]),
    "&", DIDci(did_disp_results$did3_lnrob_coefs[,2]),
    "&", DIDci(did_disp_results$did3_lnass_coefs[,2]),
    "\\\\ \\\\",   
    "Pre-Legalization Mean (Dispensaries)",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_stack$actual_robbery_with_a_gun[l], ucr_stack$population[l])), "}",
    "& \\multicolumn{1}{c}{", f(weighted.mean(ucr_stack$actual_assault_with_a_gun[l], ucr_stack$population[l])), "}",
    "& \\multicolumn{1}{c}{--}",
    "& \\multicolumn{1}{c}{--}",
    "\\\\ \\\\",
    # "Implied Percent Change ",
    # "& \\multicolumn{1}{c}{", f(100*did_disp_results$did3_rob_coefs[1,1]/weighted.mean(ucr_stack$actual_robbery_with_a_gun[l], ucr_stack$population[l])), "\\%}",
    # "& \\multicolumn{1}{c}{", f(100*did_disp_results$did3_ass_coefs[1,1]/weighted.mean(ucr_stack$actual_assault_with_a_gun[l], ucr_stack$population[l])), "\\%}",
    # "& \\multicolumn{1}{c}{--}",
    # "& \\multicolumn{1}{c}{--}",
    # "\\\\",
    # "&", DIDci_per(did_disp_results$did3_rob_coefs[,1], m = weighted.mean(ucr_stack$robbery[l], ucr_stack$population[l])),
    # "&", DIDci_per(did_disp_results$did3_ass_coefs[,1], m = weighted.mean(ucr_stack$assaults[l], ucr_stack$population[l])),
    # "&", DIDci_per(did_disp_results$did3_lnrob_coefs[,1], m = weighted.mean(sui_dat_stack$num_homi[k], sui_dat_stack$total_population[k])),
    # "&", DIDci_per(did_disp_results$did3_lnass_coefs[,1], m = weighted.mean(sui_dat_stack$num_sui[k], sui_dat_stack$total_population[k])),
    # "\\\\ \\\\",
    # "Average Pre-trend Violation: &", f(weighted.mean(did_disp_results$did3_rob_coefs[1,2:5], 1/colSds(did_disp_results$did3_rob_coefs[,2:5]))), 
    # "&", f(weighted.mean(did_disp_results$did3_ass_coefs[1,2:5], 1/colSds(did_disp_results$did3_ass_coefs[,2:5]))),
    # "&", f(weighted.mean(did_disp_results$did3_lnrob_coefs[1,2:5], 1/colSds(did_disp_results$did3_lnrob_coefs[,2:5]))),
    # "&", f(weighted.mean(did_disp_results$did3_lnass_coefs[1,2:5], 1/colSds(did_disp_results$did3_lnass_coefs[,2:5]))),
    # "\\\\ \\\\",    
    # # "SD of Pre-trend Violations: &", f(sqrt(wtd.var(did_disp_results$did3_rob_coefs[1,2:5], 1/colSds(did_disp_results$did3_rob_coefs[,2:5]), normwt = TRUE))), 
    # # "&", f(sqrt(wtd.var(did_disp_results$did3_ass_coefs[1,2:5], 1/colSds(did_disp_results$did3_ass_coefs[,2:5]), normwt = TRUE))),
    # # "&", f(sqrt(wtd.var(did_disp_results$did3_lnrob_coefs[1,2:5], 1/colSds(did_disp_results$did3_lnrob_coefs[,2:5]), normwt = TRUE))),
    # # "&", f(sqrt(wtd.var(did_disp_results$did3_lnass_coefs[1,2:5], 1/colSds(did_disp_results$did3_lnass_coefs[,2:5]), normwt = TRUE))),
    # # "\\\\ \\\\",
    # "Slope of Pre-trend Violations: &", pre_trend_slope(did_disp_results$did3_rob_coefs[,2:5]), 
    # "&", pre_trend_slope(did_disp_results$did3_ass_coefs[,2:5]), 
    # "&", pre_trend_slope(did_disp_results$did3_lnrob_coefs[,2:5]), 
    # "&", pre_trend_slope(did_disp_results$did3_lnass_coefs[,2:5]), 
    # "\\\\",
    # "&", pre_trend_slope_ci(did_disp_results$did3_rob_coefs[,2:5]), 
    # "&", pre_trend_slope_ci(did_disp_results$did3_ass_coefs[,2:5]), 
    # "&", pre_trend_slope_ci(did_disp_results$did3_lnrob_coefs[,2:5]), 
    # "&", pre_trend_slope_ci(did_disp_results$did3_lnass_coefs[,2:5]),
    # "\\\\ \\\\",
    # "Largest Pre-trend Violation: &", DIDcoef(did_disp_results$did3_rob_coefs[,(max_pt[1]+2)]), 
    # "&", DIDcoef(did_disp_results$did3_ass_coefs[,(max_pt[2]+2)]),
    # "&", DIDcoef(did_disp_results$did3_lnrob_coefs[,(max_pt[3]+2)]),
    # "&", DIDcoef(did_disp_results$did3_lnass_coefs[,(max_pt[4]+2)]),
    # "\\\\",
    # "&", DIDci(did_disp_results$did3_rob_coefs[,(max_pt[1]+2)]), 
    # "&", DIDci(did_disp_results$did3_ass_coefs[,(max_pt[2]+2)]),
    # "&", DIDci(did_disp_results$did3_lnrob_coefs[,(max_pt[3]+2)]),
    # "&", DIDci(did_disp_results$did3_lnass_coefs[,(max_pt[4]+2)]),
    # "\\\\ \\\\",
    # "Non-Parallel Trends Interval: &",  DIDcred(did_disp_results$did3_rob_coefs[,-c(1,2)]),
    # "&", DIDcred(did_disp_results$did3_ass_coefs[,-c(1,2)]),
    # "&", DIDcred(did_disp_results$did3_lnrob_coefs[,-c(1,2)]),
    # "&", DIDcred(did_disp_results$did3_lnass_coefs[,-c(1,2)]),
    "\\\\",
    " \\hline ",
    "\\\\",
    "Gun Law Controls? & No & No & No & No",
    "\\\\[0.5mm]",  
    "Non-Gun Crime Type? & Knife & Knife & Knife & Knife",
    "\\\\[0.5mm]",    
    "Unit of Observation ",
    "& \\multicolumn{4}{c}{ Jurisdiction - Year - Type of Offense}",
    "\\\\[0.5mm]",
    "Number of Observations ",
    "& \\multicolumn{4}{c}{", f(nrow(ucr_stack)), "}",
    "\\\\[0.5mm]",
    "\\bottomrule",
    sep = "",
    file = paste0(tabPath, "did_disp.tex"))