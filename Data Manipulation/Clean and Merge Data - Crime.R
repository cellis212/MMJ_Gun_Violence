############################################
#
# This file loads and cleans the UCR data
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

# Load in functions and options file.
source(file = "./Code/Functions_and_Options.R")

# Loading in Data ---------------------------------------------------------
# Timeline
timeline <- read.csv(file = "./Data/Timeline of MJ Laws.csv", fileEncoding="UTF-8-BOM")

# UCR Data
load(file = "./Data/Uniform Crime Reports/offenses_known_yearly_1960_2020.RData")


# Loading Laws ------------------------------------------------------------
law_dat <- read.csv("./Data/Gun Law Enactments and Repeals.csv")


# Merging -----------------------------------------------------------------

# Join the UCR data with the timeline data
ucr_dat <- 
  offenses_known_yearly_1960_2020 %>% 
  inner_join(timeline, by = c("state_abb" = "State"))


ucr_dat %>% 
  group_by(year) %>% 
  summarise(rob_per = sum(actual_robbery_with_a_gun)/sum(actual_robbery_with_a_gun + actual_robbery_with_a_knife + actual_robbery_other_weapon + actual_robbery_unarmed),
            ass_per = sum(actual_assault_with_a_gun)/sum(actual_assault_with_a_gun + actual_assault_with_a_knife + actual_assault_other_weapon + actual_robbery_unarmed),
            totn = sum(actual_robbery_with_a_gun + actual_assault_with_a_gun)) %>% tail()


# Join the UCR data with the law data
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
ucr_dat <- 
  ucr_dat %>% 
  inner_join(law_dat, by = c("state_abb" = "state"))

# Create binary variables for laws in each year
ucr_dat <- 
  ucr_dat %>% 
  mutate(registration_required = (year >= registration_required_enact & year < registration_required_repeal),
         stand_your_ground = (year >= stand_your_ground_enact & year < stand_your_ground_repeal),
         ccw_shall_issue = (year >= ccw_shall_issue_enact & year < ccw_shall_issue_repeal),
         ccw_prohibit = (year >= ccw_prohibit_enact & year < ccw_prohibit_repeal),
         backgrd_check3 = (year >= backgrd_check3_enact & year < backgrd_check3_repeal),
         backgrd_check2 = (year >= backgrd_check2_enact & year < backgrd_check2_repeal),
         backgrd_check1 = (year >= backgrd_check1_enact & year < backgrd_check1_repeal))



# Drop oris with no variation
# no_var <- ucr_dat %>%
#   group_by(ori) %>%
#   summarise(dist = max(actual_robbery_with_a_gun) - min(actual_robbery_with_a_gun)) %>%
#   filter(dist == 0)
# 
# ucr_dat <- ucr_dat %>% filter(!(ori %in% no_var$ori))



# New Vars ----------------------------------------------------------------
ucr_dat$treat <- case_when(ucr_dat$year > year(mdy(ucr_dat$First.MMJ.Dispensary.Opening.Date)) ~ 1,
                           ucr_dat$year == year(mdy(ucr_dat$First.MMJ.Dispensary.Opening.Date)) ~ (decimal_date(mdy(ucr_dat$First.MMJ.Dispensary.Opening.Date)) - year(mdy(ucr_dat$First.MMJ.Dispensary.Opening.Date))),
                           ucr_dat$year < year(mdy(ucr_dat$First.MMJ.Dispensary.Opening.Date)) ~ 0,
                           ucr_dat$First.MMJ.Dispensary.Opening.Date == "" ~ 0)

ucr_dat$treat_rec <- case_when(ucr_dat$year > year(mdy(ucr_dat$First.Recreational.Sale.Date)) ~ 1,
                               ucr_dat$year == year(mdy(ucr_dat$First.Recreational.Sale.Date)) ~ (decimal_date(mdy(ucr_dat$First.Recreational.Sale.Date)) - year(mdy(ucr_dat$First.Recreational.Sale.Date))),
                               ucr_dat$year < year(mdy(ucr_dat$First.Recreational.Sale.Date)) ~ 0,
                               ucr_dat$First.Recreational.Sale.Date == "" ~ 0)


# Drop missing treat, DC, and CA
ucr_dat <- ucr_dat %>% 
  filter(!(is.na(treat)),
         state_abb != "DC") 

# Keep non-negative and year > 1998 and drop data with missing months and zero pop
ucr_dat <- ucr_dat %>% filter(actual_assault_with_a_gun >= 0 &
                                actual_robbery_with_a_gun >= 0 &
                                actual_assault_with_a_knife >= 0 & 
                                actual_robbery_with_a_knife >= 0 &
                                actual_assault_other_weapon >= 0 & 
                                actual_robbery_other_weapon >= 0 &
                                actual_assault_unarmed >= 0 & 
                                actual_robbery_unarmed >= 0 &
                                year > 1998 &
                                year < 2020 &
                                number_of_months_missing == 0 &
                                population > 0)

# Winsorize
ucr_dat$actual_assault_with_a_gun_win <- ucr_dat$actual_assault_with_a_gun %>% Winsorize()
ucr_dat$actual_robbery_with_a_gun_win <- ucr_dat$actual_robbery_with_a_gun %>% Winsorize()
ucr_dat$actual_assault_with_a_knife_win <- ucr_dat$actual_assault_with_a_knife %<>% Winsorize()
ucr_dat$actual_robbery_with_a_knife_win <- ucr_dat$actual_robbery_with_a_knife %<>% Winsorize()
ucr_dat$actual_assault_unarmed_win <- ucr_dat$actual_assault_actual_assault_unarmed %<>% Winsorize()
ucr_dat$actual_robbery_unarmed_win <- ucr_dat$actual_robbery_unarmed %<>% Winsorize()
ucr_dat$actual_assault_other_weapon_win <- ucr_dat$actual_assault_other_weapon %<>% Winsorize()
ucr_dat$actual_robbery_other_weapon_win <- ucr_dat$actual_robbery_other_weapon %<>% Winsorize()

# Create Log Vars
ucr_dat <- ucr_dat %>%
  mutate(
    ln_ass_gun = log(1 + actual_assault_with_a_gun),
    ln_rob_gun = log(1 + actual_robbery_with_a_gun),
    ln_ass_knife = log(1 + actual_assault_with_a_knife),
    ln_rob_knife = log(1 + actual_robbery_with_a_knife),
    ln_ass_unarmed = log(1 + actual_assault_unarmed),
    ln_rob_unarmed = log(1 + actual_robbery_unarmed),
    ln_ass_other_weapon = log(1 + actual_assault_other_weapon),
    ln_rob_other_weapon = log(1 + actual_robbery_other_weapon)
  )

# Get first treat year
ucr_dat$first_treat <- year(mdy(ucr_dat$First.MMJ.Dispensary.Opening.Date))

ucr_dat$year_since_treat <- 
  case_when(ucr_dat$year - ucr_dat$first_treat <= -5 ~ "-5+",
            ucr_dat$year - ucr_dat$first_treat >= 6 ~ "6+",
            is.na(ucr_dat$first_treat) ~ "-5+",
            TRUE ~ as.character(ucr_dat$year - ucr_dat$first_treat)) %>% 
  factor(levels = c("-1", "-5+", "-4", "-3", "-2", "0", "1", "2", "3", "4","5","6+"))


# Create a copy of the original dataset, ucr_dat, and assign it to a new variable, temp1
temp1 <- 
  ucr_dat %>% 
  mutate(type = "gun",
         assault = actual_assault_with_a_gun_win,
         robbery = actual_robbery_with_a_gun_win,
         ln_ass = log(1 + actual_assault_with_a_gun),
         ln_rob = log(1 + actual_robbery_with_a_gun))

# Create another copy of the original dataset, ucr_dat, and assign it to a new variable, temp2
temp2 <-
ucr_dat %>% 
  mutate(type = "knife",
         assault = actual_assault_with_a_knife_win,
         robbery = actual_robbery_with_a_knife_win,
         ln_ass = log(1 + actual_assault_with_a_knife),
         ln_rob = log(1 + actual_robbery_with_a_knife))

# Create another copy of the original dataset, ucr_dat, and assign it to a new variable, temp2
temp3 <-
  ucr_dat %>% 
  mutate(type = "unarmed",
         assault = actual_assault_unarmed_win,
         robbery = actual_robbery_unarmed_win,
         ln_ass = log(1 + actual_assault_unarmed),
         ln_rob = log(1 + actual_robbery_unarmed))

# Create another copy of the original dataset, ucr_dat, and assign it to a new variable, temp2
temp4 <-
  ucr_dat %>% 
  mutate(type = "other",
         assault = actual_assault_other_weapon_win,
         robbery = actual_robbery_other_weapon_win,
         ln_ass = log(1 + actual_assault_other_weapon),
         ln_rob = log(1 + actual_robbery_other_weapon))

# Combine the rows of temp1 - temp4 into a single dataframe, ucr_dat_stack
ucr_dat_stack <- rbind(temp1, temp2, temp3, temp4)


# Create a new variable 'trip_treat' in the ucr_dat_stack dataframe.
# If the 'treat' column is 1 and the 'type' column is 'gun', set 'trip_treat' to 1, otherwise set it to 0.
ucr_dat_stack$trip_treat <- if_else(ucr_dat_stack$treat == 1 &
                                      ucr_dat_stack$type == "gun", 1, 0)

# Create a new variable 'trip_year_since_treat' in the ucr_dat_stack dataframe.
# If the 'type' column is 'gun', concatenate the 'year_since_treat' column with 'g', otherwise concatenate it with 'o'.
# The resulting variable is then converted into a factor with specific levels.
ucr_dat_stack$trip_year_since_treat <- if_else(ucr_dat_stack$type == "gun", 
                                               paste0(ucr_dat_stack$year_since_treat, "g"), 
                                               paste0(ucr_dat_stack$year_since_treat, "o")) %>%
  factor(levels = c("-1g", "-5+g", "-4g", "-3g", "-2g", "0g", "1g", "2g", "3g", "4g","5g", "6+g",
                    "-5+o", "-4o", "-3o", "-2o", "0o", "1o", "2o", "3o", "4o","5o", "6+o", "-1o"))
# factor(levels = c("-1g",  "-4+g", "-3g", "-2g", "0g", "1g", "2g", "3g", "4+g",
#                   "-4+o", "-3o", "-2o", "0o", "1o", "2o", "3o", "4+o", "-1o"))


# Get interaction FEs
ucr_dat_stack$year_ori <- paste(ucr_dat_stack$year, ucr_dat_stack$ori)
ucr_dat_stack$type_ori <- paste(ucr_dat_stack$type, ucr_dat_stack$ori)
ucr_dat_stack$type_year <- paste(ucr_dat_stack$type, ucr_dat_stack$year)



# Save data for spec curve
save(ucr_dat_stack, file = "./Data/Uniform Crime Reports/ucr_dat_specs.RData")

# Dispensary Data ---------------------------------------------------------
require(haven)
require(readxl)
require(lubridate)

# Matching dispensaries to zips within 25 miles ---------------------------
# Loading ZipDistance data
Distance <- read.csv(file = "./Data/ZipDistances_2010_FieldKey/gaz2016zcta5distance25miles.csv")
# Distance <- Distance %>% filter(mi_to_zcta5 < 10)

# Initializing
disZip <- data.frame("zip" = unique(Distance$zip1), "year" = 2010, "disp25" = 0, "everDisp" = 0)
temp <- disZip
for(year in 2011:2019){
  temp$year <- year
  disZip <- rbind(disZip, temp)
}

# Loading dispensary data
disData <- read_excel("./Data/NewStates_Dispensary_2019update.xlsx")
disData$zip <- as.numeric(disData$zip)

ucr_dat_stack %>% names()

# Adding in new dis data
newDis <- read_stata(file = "./Data/WEEDMAPS OKPA_ZIPS.dta")
l <- which(newDis$opening_date == "Coming Soon")
newDis$opening_date[l] <- "01/01/9999"
newDis %<>% filter(opening_date != "")
newDis$year <- substr(newDis$opening_date, nchar(newDis$opening_date) - 1, nchar(newDis$opening_date)) %>% as.numeric()
newDis$year <- newDis$year + 2000

# Ever have a dispensary
l <- which(Distance$zip1 %in% disData$zip |
             Distance$zip1 %in% newDis$zip)

j <- which(disZip$zip %in% Distance$zip1[l] |
             disZip$zip %in% Distance$zip2[l])

disZip$everDisp[j] <- 1

# Open dispensary within 25 miles
for(i in 1:nrow(disData)){
  l <- which(Distance$zip1 == disData$zip[i])
  j <- which((disZip$zip %in% Distance$zip1[l] |
                disZip$zip %in% Distance$zip2[l]) &
               disZip$year >= disData$Yearopened[i])


  disZip$disp25[j] <- 1
  print(i/nrow(disData))
}

for(i in 1:nrow(newDis)){
  l <- which(Distance$zip1 == newDis$zip[i])
  j <- which((disZip$zip %in% Distance$zip1[l] |
                disZip$zip %in% Distance$zip2[l]) &
               disZip$year >= newDis$year[i])

  disZip$disp25[j] <- 1


  print(i/nrow(newDis))
}



ucr_dat_stack <- ucr_dat_stack %>% left_join(disZip, by = c("zip_code" = "zip", "year"="year")) %>% mutate(disp_eligible = if_else(first_treat < 2014, 0, 1, 1))


############## Getting the registered patients number ################
regist <- read.csv(file = "./Data/Medical MJ Registered Patients/Medical Marijuana Registered Patients.csv",
                   stringsAsFactors = FALSE)




# Making numeric
for(i in 3:9){
  l <- which(regist[, i] == "")
  regist[l, i] <- "0"
  regist[, i] <- gsub(",", "", regist[, i]) # dropping commas
  regist[, i] <- gsub(" ", "", regist[, i]) # dropping spaces
  regist[, i] <- gsub("missing", "-1", regist[, i]) # changing missing to -1
  regist[, i] <- as.numeric(regist[, i])
}

# Imputing missing data
regist$Registered.2014[39] <- regist$Registered.2015[39] + (regist$Registered.2015[39] - regist$Registered.2016[39]) # Oregon special case 
regist$x2017[19] <- regist$x2018[19]/2 # Louisiana special case

l <- which(regist$Registered.2013 == -1)
regist$X2013[l] <- regist$x2014[l]/2

l <- which(regist$Registered.2014 == -1)
regist$Registered.2014 <- (regist$Registered.2013[l] + regist$Registered.2015[l])/2

l <- which(regist$Registered.2015 == -1)
regist$Registered.2015[l] <- (regist$Registered.2014[l] + regist$Registered.2016[l])/2

l <- which(regist$Registered.2016 == -1)
regist$Registered.2016[l] <- (regist$Registered.2015[l] + regist$Registered.2017[l])/2

l <- which(regist$Registered.2017 == -1)
regist$Registered.2017[l] <- (regist$Registered.2016[l] + regist$Registered.2018[l])/2

l <- which(regist$Registered.2018 == -1)
regist$Registered.2018[l] <- regist$Registered.2017[l] + (regist$Registered.2017[l] - regist$Registered.2016[l])

ucr_dat_stack$registered <- 0
for(i in unique(ucr_dat_stack$State)){
  for(year in 2014:2019){
    m <- which(regist$State == i)
    l <- which(ucr_dat_stack$state_abb == i &
                 ucr_dat_stack$year == year)
    
    if(length(m) == 0) {
      print(i) 
      print(year)
    }
    if(length(l) == 0) {
      print(i) 
      print(year)
    }
    
    if(year == 2014) ucr_dat_stack$registered[l] <- regist$Registered.2014[m]
    if(year == 2015) ucr_dat_stack$registered[l] <- regist$Registered.2015[m]
    if(year == 2016) ucr_dat_stack$registered[l] <- regist$Registered.2016[m]
    if(year == 2017) ucr_dat_stack$registered[l] <- regist$Registered.2017[m]
    if(year == 2018) ucr_dat_stack$registered[l] <- regist$Registered.2018[m]
    if(year == 2019) ucr_dat_stack$registered[l] <- regist$Registered.2019[m]
    
  }
}

# Logged number registered
l <- which(ucr_dat_stack$registered > 0)
ucr_dat_stack$lregis <- 0
ucr_dat_stack$lregis[l] <- log(ucr_dat_stack$registered[l])




# Saving
save(ucr_dat_stack, file = "./Data/ucr_data_stack_cleaned.RData")
write.csv(ucr_dat_stack, file = "./Data/ucr_data_stack_cleaned.csv")

# Subsets ---------------------------------------------------

load("./Data/ucr_data_stack_cleaned.RData")


# Make NA zips have no Dispensaries
ucr_dat_stack$everDisp <- if_else(is.na(ucr_dat_stack$everDisp), 0, ucr_dat_stack$everDisp)
ucr_dat_stack$disp25 <- if_else(is.na(ucr_dat_stack$disp25), 0, ucr_dat_stack$disp25)

# Remove missing data controls
ucr_dat_stack <- 
  ucr_dat_stack %>% 
  select(treat, ori, year, 
         registration_required, stand_your_ground, ccw_shall_issue,
         ccw_prohibit, backgrd_check3, backgrd_check2, backgrd_check1, 
         population, year_since_treat, 
         actual_robbery_with_a_gun, actual_assault_with_a_gun,
         actual_robbery_with_a_knife, actual_assault_with_a_knife,
         ln_rob_gun, ln_ass_gun, ln_rob_knife, ln_ass_knife, 
         ln_rob, ln_ass, robbery, assault, state,
         actual_robbery_unarmed, actual_assault_unarmed,
         ln_rob_unarmed, ln_ass_unarmed, actual_robbery_other_weapon, 
         actual_assault_other_weapon,
         ln_rob_other_weapon, ln_ass_other_weapon, lregis, 
         disp_eligible, disp25, everDisp, treat_rec, type, trip_treat, 
         trip_year_since_treat, year_since_treat, type_ori, type_year, year_ori,
         agency_type) %>% 
  filter(complete.cases(.))

# Save for spec curve
save(ucr_dat_stack, file = "./Data/Uniform Crime Reports/ucr_dat_specs.RData")

# Drop super long treat
ucr_dat_stack %<>% filter(year_since_treat != "6+")

# Drop rec states
ucr_dat_stack <-
  ucr_dat_stack %>%
  filter(treat_rec == 0)

# Drop Hawaii and Alaska
ucr_dat_stack <- ucr_dat_stack %>% filter(!(state %in% c("hawaii", "alaska")))

# Only local police and drop other and unarmed (eventually go back and remove from earlier code).
ucr_dat_stack <- ucr_dat_stack %>% filter(type != "other" &
                                            type != "unarmed")

# Saving
save(ucr_dat_stack, file = "./Data/ucr_data_stack_cleaned_subset.RData")
write.csv(ucr_dat_stack, file = "./Data/ucr_data_stack_cleaned_subset.csv")

# # Load suicide data -------------------------------------------------
sui_dat <- read.csv(file = "./Data/Firearm_suicide_homicide_dataset.csv", fileEncoding="UTF-8-BOM")
sui_dat %>% group_by(year) %>% summarise(num = sum(firearm_homicides + firearm_suicides)) %>% tail()

# 
# 
# 
# # Merge with UCR data for controls and treat vars
# sui_dat <- 
#   ucr_dat %>%
#   select(registration_required, stand_your_ground, ccw_shall_issue, ccw_prohibit, 
#          backgrd_check3, backgrd_check2, backgrd_check1,
#          treat, treat_rec, year_since_treat, stname, year, first_treat) %>% 
#   group_by(stname, year) %>% 
#   summarise(across(everything(), mean, na.rm = TRUE)) %>% 
#   right_join(sui_dat, by = c("stname" = "state", "year" = "year")) %>% 
#   filter(year > 1998 & 
#            stname != "District of Columbia" &
#            stname != "California")
# 
# # New Vars
# sui_dat <- 
#   sui_dat %>%
#   mutate(num_sui = firearm_suicides,
#          num_homi = firearm_homicides,
#          num_sui_o = total_suicides - firearm_suicides,
#          num_homi_o = nonfirearm_homicides,
#          ln_num_sui = log(firearm_suicides + 1),
#          ln_num_homi = log(firearm_homicides + 1),
#          ln_num_sui_o = log(total_suicides - firearm_suicides + 1),
#          ln_num_homi_o = log(nonfirearm_homicides + 1),
#          state = stname)
# 
# # Fix year since treat
# sui_dat$year_since_treat <- 
#   case_when(sui_dat$year - sui_dat$first_treat <= -5 ~ "-5+",
#             sui_dat$year - sui_dat$first_treat >= 6 ~ "6+",
#             is.na(sui_dat$first_treat) ~ "-5+",
#             TRUE ~ as.character(sui_dat$year - sui_dat$first_treat)) %>% 
#   factor(levels = c("-1", "-5+", "-4", "-3", "-2", "0", "1", "2", "3", "4","5","6+"))
#   
# # Create a copy of the original dataset, sui_dat, and assign it to a new variable, temp1
# temp1 <- 
#   sui_dat %>% 
#   mutate(type = "gun",
#          num_homi = firearm_homicides,
#          num_sui = firearm_suicides,
#          ln_num_homi = log(1 + firearm_homicides),
#          ln_num_sui = log(1 + firearm_suicides))
# 
# # Create another copy of the original dataset, sui_dat, and assign it to a new variable, temp2
# temp2 <-
#   sui_dat %>% 
#   mutate(type = "not gun",
#          num_homi = nonfirearm_homicides,
#          num_sui =  total_suicides - firearm_suicides,
#          ln_num_homi = log(1 + nonfirearm_homicides),
#          ln_num_sui = log(1 + total_suicides - firearm_suicides))
# 
# # Combine the rows of temp1 and temp2 into a single dataframe, sui_dat_stack
# sui_dat_stack <- rbind(temp1, temp2)
# 
# 
# # Create a new variable 'trip_treat' in the sui_dat_stack dataframe.
# # If the 'treat' column is 1 and the 'type' column is 'gun', set 'trip_treat' to 1, otherwise set it to 0.
# sui_dat_stack$trip_treat <- if_else(sui_dat_stack$treat == 1 &
#                                       sui_dat_stack$type == "gun", 1, 0)
# 
# # Create a new variable 'trip_year_since_treat' in the sui_dat_stack dataframe.
# # If the 'type' column is 'gun', concatenate the 'year_since_treat' column with 'g', otherwise concatenate it with 'o'.
# # The resulting variable is then converted into a factor with specific levels.
# sui_dat_stack$trip_year_since_treat <- if_else(sui_dat_stack$type == "gun", 
#                                                paste0(sui_dat_stack$year_since_treat, "g"), 
#                                                paste0(sui_dat_stack$year_since_treat, "o")) %>%
#   factor(levels = c("-1g", "-5+g", "-4g", "-3g", "-2g", "0g", "1g", "2g", "3g", "4g","5g", "6+g",
#                     "-5+o", "-4o", "-3o", "-2o", "0o", "1o", "2o", "3o", "4o","5o", "6+o", "-1o"))
# # factor(levels = c("-1g",  "-4+g", "-3g", "-2g", "0g", "1g", "2g", "3g", "4+g",
# #                   "-4+o", "-3o", "-2o", "0o", "1o", "2o", "3o", "4+o", "-1o"))
# 
# 
# # Get interaction FEs
# sui_dat_stack$year_state <- paste(sui_dat_stack$year, sui_dat_stack$state)
# sui_dat_stack$type_state <- paste(sui_dat_stack$type, sui_dat_stack$state)
# sui_dat_stack$type_year <- paste(sui_dat_stack$type, sui_dat_stack$year)
# 
# # Saving
# save(sui_dat_stack, file = "./Data/sui_data_stack_cleaned.RData")
# write.csv(sui_dat_stack, file = "./Data/GPT/sui_data_stack_cleaned.csv")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 


