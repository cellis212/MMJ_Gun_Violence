############################################
#
# This file loads and cleans the trace data
# And merges it with the legalization variables
#
#############################################


# Preamble ----------------------------------------------------------------

# Clearing the workspace
rm(list = ls())

# Loading the packages that will be needed. 
require(tidyverse)
require(lubridate)
library(reticulate)
require(magrittr)
require(DescTools)
library(readxl)

# Load in functions and options file.
source(file = "./Code/Functions_and_Options.R")

# Loading in Data ---------------------------------------------------------
# Timeline
timeline <- read.csv(file = "./Data/Timeline of MJ Laws.csv", fileEncoding="UTF-8-BOM")

# Gun Trace Data
trace_dat <- read_xlsx("./Data/Aggregated Trace Data.xlsx")

# Loading Laws ------------------------------------------------------------
law_dat <- read.csv("./Data/Gun Law Enactments and Repeals.csv")

# Getting states in the same form
timeline$StateName <- tolower(timeline$stname)
trace_dat$StateName <- tolower(trace_dat$StateName)


# Merging -----------------------------------------------------------------

# Join the trace data with the timeline data
trace_dat <- 
  trace_dat %>% 
  inner_join(timeline, by = c("StateName" = "StateName"))


# Join the trace data with the law data
# Function to turn NA into 9999
na_to_9999 <- function(x) {
  x[is.na(x)] <- 9999
  return(x)
}


# Loop over law_dat columns and turn NAs into 9999
for (i in 1:ncol(law_dat)) {
  law_dat[,i] <- na_to_9999(law_dat[,i])
}

# Merge with law data
trace_dat <- 
  trace_dat %>% 
  inner_join(law_dat, by = c("State" = "state"))

# Create binary variables for laws in each year
trace_dat <- 
  trace_dat %>% 
  mutate(registration_required = (year >= registration_required_enact & year < registration_required_repeal),
         stand_your_ground = (year >= stand_your_ground_enact & year < stand_your_ground_repeal),
         ccw_shall_issue = (year >= ccw_shall_issue_enact & year < ccw_shall_issue_repeal),
         ccw_prohibit = (year >= ccw_prohibit_enact & year < ccw_prohibit_repeal),
         backgrd_check3 = (year >= backgrd_check3_enact & year < backgrd_check3_repeal),
         backgrd_check2 = (year >= backgrd_check2_enact & year < backgrd_check2_repeal),
         backgrd_check1 = (year >= backgrd_check1_enact & year < backgrd_check1_repeal))


# New Vars ----------------------------------------------------------------
trace_dat$treat <- case_when(trace_dat$year > year(mdy(trace_dat$First.MMJ.Dispensary.Opening.Date)) ~ 1,
                           trace_dat$year == year(mdy(trace_dat$First.MMJ.Dispensary.Opening.Date)) ~ (decimal_date(mdy(trace_dat$First.MMJ.Dispensary.Opening.Date)) - year(mdy(trace_dat$First.MMJ.Dispensary.Opening.Date))),
                           trace_dat$year < year(mdy(trace_dat$First.MMJ.Dispensary.Opening.Date)) ~ 0,
                           trace_dat$First.MMJ.Dispensary.Opening.Date == "" ~ 0)

trace_dat$treat_rec <- case_when(trace_dat$year > year(mdy(trace_dat$First.Recreational.Sale.Date)) ~ 1,
                               trace_dat$year == year(mdy(trace_dat$First.Recreational.Sale.Date)) ~ (decimal_date(mdy(trace_dat$First.Recreational.Sale.Date)) - year(mdy(trace_dat$First.Recreational.Sale.Date))),
                               trace_dat$year < year(mdy(trace_dat$First.Recreational.Sale.Date)) ~ 0,
                               trace_dat$First.Recreational.Sale.Date == "" ~ 0)


# Drop missing treat, DC, and CA
trace_dat <- trace_dat %>% 
  filter(!(is.na(treat)),
         State != "DC") 

# Winsorize
# trace_dat$actual_assault_with_a_gun_win <- trace_dat$actual_assault_with_a_gun %>% Winsorize()
# trace_dat$actual_robbery_with_a_gun_win <- trace_dat$actual_robbery_with_a_gun %>% Winsorize()
# trace_dat$actual_assault_with_a_knife_win <- trace_dat$actual_assault_with_a_knife %<>% Winsorize()
# trace_dat$actual_robbery_with_a_knife_win <- trace_dat$actual_robbery_with_a_knife %<>% Winsorize()
# trace_dat$actual_assault_unarmed_win <- trace_dat$actual_assault_actual_assault_unarmed %<>% Winsorize()
# trace_dat$actual_robbery_unarmed_win <- trace_dat$actual_robbery_unarmed %<>% Winsorize()
# trace_dat$actual_assault_other_weapon_win <- trace_dat$actual_assault_other_weapon %<>% Winsorize()
# trace_dat$actual_robbery_other_weapon_win <- trace_dat$actual_robbery_other_weapon %<>% Winsorize()

# Create Log Vars
trace_dat <- trace_dat %>%
  mutate(year1 = Under3 + month3to7 + month7toyear1) %>%
  mutate(year1 = Winsorize(year1),
         year3 = Winsorize(year3plus)) %>% 
  mutate(
    year1_ln = log(1 + year1),
    Under3_ln = log(1 + Under3),
    month3to7_ln = log(1 + month3to7),
    month7toyear1_ln = log(1 + month7toyear1),
    year1to2_ln = log(1 + year1to2),
    year2to3_ln = log(1 + year2to3),
    year3_ln = log(1 + year3plus)
  )

# Get first treat year
trace_dat$first_treat <- year(mdy(trace_dat$First.MMJ.Dispensary.Opening.Date))

trace_dat$year_since_treat <- 
  case_when(trace_dat$year - trace_dat$first_treat <= -5 ~ "-5+",
            trace_dat$year - trace_dat$first_treat >= 3 ~ "3+",
            is.na(trace_dat$first_treat) ~ "-5+",
            TRUE ~ as.character(trace_dat$year - trace_dat$first_treat)) %>% 
  factor(levels = c("-1", "NT", "-5+", "-4", "-3", "-2", "0", "1", "2", "3+"))


# Create a copy of the original dataset, trace_dat, and assign it to a new variable, temp1
temp1 <- 
  trace_dat %>% 
  mutate(type = "year1",
         guns = year1_ln)

# Create another copy of the original dataset, trace_dat, and assign it to a new variable, temp2
temp2 <-
trace_dat %>% 
  mutate(type = "year3",
         guns = year3_ln)

# Combine the rows of temp1 - temp4 into a single dataframe, trace_dat_stack
trace_dat_stack <- rbind(temp1, temp2)


# Create a new variable 'trip_treat' in the trace_dat_stack dataframe.
# If the 'treat' column is 1 and the 'type' column is 'gun', set 'trip_treat' to 1, otherwise set it to 0.
trace_dat_stack$trip_treat <- if_else(trace_dat_stack$treat == 1 &
                                      trace_dat_stack$type == "year1", 1, 0)

# Create a new variable 'trip_year_since_treat' in the trace_dat_stack dataframe.
# If the 'type' column is 'gun', concatenate the 'year_since_treat' column with 'g', otherwise concatenate it with 'o'.
# The resulting variable is then converted into a factor with specific levels.
trace_dat_stack$trip_year_since_treat <- if_else(trace_dat_stack$type == "year1", 
                                               paste0(trace_dat_stack$year_since_treat, "new"), 
                                               paste0(trace_dat_stack$year_since_treat, "old")) %>%
  factor(levels = c("-1new", "-5+new", "-4new", "-3new", "-2new", "0new", "1new", "2new", "3+new",
                    "-5+old", "-4old", "-3old", "-2old", "0old", "1old", "2old", "3+old",  "-1old"))


# Get interaction FEs
trace_dat_stack$year_state <- paste(trace_dat_stack$year, trace_dat_stack$State)
trace_dat_stack$type_state <- paste(trace_dat_stack$type, trace_dat_stack$State)
trace_dat_stack$type_year <- paste(trace_dat_stack$type, trace_dat_stack$year)



# Save data for spec curve
save(trace_dat_stack, file = "./Data/Uniform Crime Reports/trace_dat_specs.RData")


# Saving
save(trace_dat_stack, file = "./Data/trace_data_stack_cleaned.RData")
write.csv(trace_dat_stack, file = "./Data/trace_data_stack_cleaned.csv")

# Subsets ---------------------------------------------------

# Drop long treat
trace_dat_stack %<>% filter(year_since_treat != "3+")

# Drop rec states
trace_dat_stack <-
  trace_dat_stack %>%
  filter(treat_rec == 0)

# Drop Hawaii and Alaska
trace_dat_stack <- trace_dat_stack %>% filter(!(StateName %in% c("hawaii", "alaska")))

# Balanced panel
# trace_dat_stack <- 
#   trace_dat_stack %>% 
#   group_by(State) %>%
#   filter(n() == 26)
  
  
# Saving
save(trace_dat_stack, file = "./Data/trace_data_stack_cleaned_subset.RData")
write.csv(trace_dat_stack, file = "./Data/trace_data_stack_cleaned_subset.csv")
