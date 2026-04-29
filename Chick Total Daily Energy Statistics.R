#Statistical tests for chick total daily energy published in: CHANGES ON THE MENU: INTRA-ANNUAL SHIFTS IN DIET AND FOOD PROVISIONING OF ATLANTIC PUFFIN (FRATERCULA ARCTICA) CHICKS
#Authored by: AIMEE FRASER, CHRISTINA PETALAS, RAPHAËL A. LAVOIE, KYLE H. ELLIOTT
#2026

####Format and clean data#######################################################

rm(list=ls())

library(dplyr)
library(tidyr)
library(lme4)
library(emmeans)

###GPS data
      #Read data
        trips<-read.csv(file.choose(), header = T, na.string ="", stringsAsFactors = T)
      #Make a data frame of total trips per day that each bird went on
        trips$foraging_trip <- 1 #Create a column with a 1 to track and sum number of foraging trips per bird
        daily_data <- trips %>%
          group_by(dep_id, day_number, year) %>%
          summarise(total_trips = sum(foraging_trip))

      #Add rows for days where puffins did not go on any trips
      #First remove last row with NA
        daily_data <- daily_data[-383,]
      #Get full range of days and all bird IDs      
       all_days <- seq(min(daily_data$day_number), max(daily_data$day_number))
       all_birds <- unique(daily_data$dep_id)
      #
      #Create full grid of bird-day combinations
       full_grid <- expand.grid(dep_id = all_birds, day = all_days)
      
      #Left join with your original data
        colnames(daily_data)[colnames(daily_data) == 'day_number'] <- 'day' #rename day_number to day so the join will work (col names will be the same)
      #
        daily_data <- full_grid %>%
          left_join(daily_data, by = c("dep_id", "day")) %>%
          mutate(total_trips = coalesce(total_trips, 0))
      #
      #Fill in years for the new 0 trip columns based on dep_id
      #Make a data frame with each year that corresponds to dep_id
        dep_years <- daily_data %>%
          filter(!is.na(year)) %>%
          distinct(dep_id, year)
      #Join with daily_data
      daily_data <- daily_data %>%
        dplyr::select(dep_id, day, total_trips, year) %>%  #keep these columns
        left_join(dep_years, by = "dep_id", suffix = c("", "_from_lookup")) %>%
        mutate(year = coalesce(year, year_from_lookup)) %>%
        dplyr::select(-year_from_lookup)
      
      #Remove rows after the last logged trip (to not over-inflate few late period trips)
      #Find each bird's last foraging day
        last_trip_day <- daily_data %>%
          filter(total_trips > 0) %>%
          group_by(dep_id) %>%
          summarise(last_day = max(day), .groups = "drop")
      #
      #Join and filter out days after last trip
        daily_data <- daily_data %>%
          left_join(last_trip_day, by = "dep_id") %>%
          filter(day <= last_day) %>%
          dplyr::select(-last_day)
      
      #Add a period column
        daily_data <- daily_data %>%
          mutate(period = case_when(
            day <= 4 ~ "1-4",
            day <= 8 ~ "5-8",
            day <= 12 ~ "9-12",
            day <= 16 ~ "13-16",
            day > 16 ~ ">16"
          ))
      #Check levels and rearrange 
        daily_data$period <- factor(daily_data$period, levels=c('1-4', '5-8', '9-12', '13-16','>16'))
        table(daily_data$period)
      
      #rename the data frame
       daily_trips <- daily_data

### Camera observation data
       
      #Read data
       puffins_raw<-read.csv(file.choose(), header = T, na.string ="", stringsAsFactors = T)
       
      #Drop extra rows and columns
       puffins_raw <- puffins_raw[-c(950:1108), ]
       puffins_raw <- subset(puffins_raw, select = -c(total_prey_length_px, total_prey_length_cm))
       
      #Fix typos
       levels(puffins_raw$species)
       which(puffins_raw$species=="Sandance")
       puffins_raw$species[147]<- "Sandlance"
       which(puffins_raw$species=="Fish larvae")
       puffins_raw$species[29]<- "Fish larva"
       which(puffins_raw$species=="Capelin  ")
       puffins_raw$species[841]<- "Capelin"
       table(puffins_raw$species)
       puffins_raw$species <- droplevels(puffins_raw$species)
       
      #Set up different data frames so we can try the statistical analyses excluding select observations
       library(dplyr)
       puffins_cleaned <- puffins_raw
     #Change species to unknown when confidence is below 50%
       levels(puffins_cleaned$species) <- c(levels(puffins_cleaned$species),"LowConfidence")
       puffins_cleaned <- puffins_cleaned %>%
         mutate(species = if_else(confidence < 50, "Unknown", species)) #Renames LowConfidence to Unknown
      #Remove rows with flying birds
       puffins_cleaned <- subset(puffins_cleaned, position  == "grounded")
      #Remove other island data
       puffins_cleaned <- subset(puffins_cleaned, island  == "Betchouanes")
      #Set data frames to either puffins_raw or puffins_cleaned
       puffins <- puffins_cleaned
       rm(puffins_raw, puffins_cleaned)
       
      #Calculation for prey species frequency of occurrence
      #Count totals of each species to know what proportions to use for each species when calculating the unknown length modifiers
       species_totals <- puffins %>%
         group_by(species, year) %>%
         summarise(total_count = sum(prey_number, na.rm = TRUE))
       
      #Unknown Spp Formulas
      #Constants for each year to calculate unknown prey lengths from tail end measurements based on species abundance
       unknown_2021 <- 0.6688*((species_totals$total_count[1])/(sum(species_totals$total_count[c(1, 4, 7)])))+
         0.5669*((species_totals$total_count[4])/(sum(species_totals$total_count[c(1,4,7)])))+
         0.7119*((species_totals$total_count[7])/(sum(species_totals$total_count[c(1,4,7)])))
      #
       unknown_2023 <-0.6688*((species_totals$total_count[2])/(sum(species_totals$total_count[c(2,5,8)])))+
         0.5669*((species_totals$total_count[5])/(sum(species_totals$total_count[c(2,5,8)])))+
         0.7119*((species_totals$total_count[8])/(sum(species_totals$total_count[c(2,5,8)])))
      #
       unknown_2024 <- 0.6688*((species_totals$total_count[3])/(sum(species_totals$total_count[c(3,6,9)])))+
         0.5669*((species_totals$total_count[6])/(sum(species_totals$total_count[c(3,6,9)])))+
         0.7119*((species_totals$total_count[9])/(sum(species_totals$total_count[c(3,6,9)])))
       
      #Species Length Formulas
      #Define species-specific formulas for each measurement type
      #Add different formulas to "Unknown" for each year because species proportion change between years
      #
       species_formula <- function(species, year, width, tail_end, head_length, length_est_px) {
         if (species == "Sandlance") {
           total_length <-  ifelse(!is.na(length_est_px), length_est_px, 0)+
             ifelse(!is.na(tail_end), tail_end/0.7119, 0)+
             ifelse(!is.na(width), width/0.0765, 0)
         } else if (species == "Capelin") {
           total_length <- ifelse(!is.na(length_est_px), length_est_px, 0)+
             ifelse(!is.na(tail_end), tail_end/0.6688, 0) +
             ifelse(!is.na(head_length), head_length/0.2028, 0)
         } else if (species == "Fish larva") {
           total_length <- ifelse(!is.na(length_est_px), length_est_px, 0)+
             ifelse(!is.na(tail_end), tail_end/0.5669, 0)
         } else if (species %in% c("Unknown", "LowConfidence")) {
           if (year == 2021) {
             total_length <- ifelse(!is.na(length_est_px), length_est_px, 0)+
               ifelse(!is.na(tail_end), tail_end/unknown_2021, 0)
           } else if (year == 2023) {
             total_length <- ifelse(!is.na(length_est_px), length_est_px, 0)+
               ifelse(!is.na(tail_end), tail_end/unknown_2023, 0)
           } else if (year == 2024) {
             total_length <- ifelse(!is.na(length_est_px), length_est_px, 0)+
               ifelse(!is.na(tail_end), tail_end/unknown_2024, 0)
           } else {
             total_length <- NA
           }
         } else {
           total_length <- NA
         }
         return(total_length)
       }
      #
      # Apply formulas and update the prey_length_px column
       puffins <- puffins %>%
         rowwise() %>%
         mutate(prey_length_px = species_formula(species, year, width, tail_end, head_length, length_est_px)) %>%
         ungroup()
      #
      #Change the "0" values to NA
       puffins$prey_length_px[puffins$prey_length_px == 0] <- NA
      #
      #Calculate prey length in cm
       length_cm <- function(year, prey_length_px, bill) {
         if (year == 2021) {
           total_length_cm <- prey_length_px * 4.935/ bill 
         } else if (year == 2023) {
           total_length_cm <- prey_length_px * 4.935/bill 
         } else if (year == 2024) {
           total_length_cm <- prey_length_px * 4.935/bill 
         } else {
           total_length_cm <- NA
         }
         return(total_length_cm)
       }
      #
       puffins <- puffins %>%
         rowwise() %>%
         mutate(row_prey_length_cm = length_cm(year, prey_length_px, bill)) %>%
         ungroup()
       
      #Remove NAs
       puffins <- subset(puffins, !is.na(row_prey_length_cm))
       
       
      #Add prey_number = 1 for all NAs
       puffins <- puffins %>% mutate(prey_number = ifelse(is.na(prey_number), 1, prey_number))
       
      #Prey mass calculations
      #Run mass calculations
       mass_g <- function(species, row_prey_length_cm, prey_number) {
         if (species == "Sandlance") {
           mass <- ifelse(!is.na(row_prey_length_cm), (((row_prey_length_cm * 10) ^ 2.11) * (209 * 10^-6)), 0)
         } else if (species == "Capelin") {
           mass <- ifelse(!is.na(row_prey_length_cm), (((row_prey_length_cm * 10) ^ 2.76) * (2.7 * 10^-5)), 0)
         } else if (species %in% c("Unknown", "LowConfidence", "Fish larva")) {
           mass <- ifelse(!is.na(row_prey_length_cm), (((row_prey_length_cm * 10) ^ 2.41) * (1.18 * 10^-4)), 0)
         } else {
           mass <- NA
         }
         return(mass)
       }
      #
      #Apply the formulas and add a mass_g column
       puffins <- puffins %>%
         rowwise() %>%
         mutate(longest_mass_g = mass_g(species, row_prey_length_cm, prey_number)) %>%
         ungroup()
      #
      #Apply for the prey load as a whole (accounting for prey quantity)
       row_mass_g <- function(species, row_prey_length_cm, prey_number) {
         if (species == "Sandlance") {
           mass <- ifelse(!is.na(row_prey_length_cm), (((row_prey_length_cm * 10) ^ 2.11) * (209 * 10^-6))*prey_number, 0)
         } else if (species == "Capelin") {
           mass <- ifelse(!is.na(row_prey_length_cm), (((row_prey_length_cm * 10) ^ 2.76) * (2.7 * 10^-5))*prey_number, 0)
         } else if (species %in% c("Unknown", "LowConfidence", "Fish larva")) {
           mass <- ifelse(!is.na(row_prey_length_cm), (((row_prey_length_cm * 10) ^ 2.41) * (1.18 * 10^-4))*prey_number, 0)
         } else {
           mass <- NA
         }
         return(mass)
       }
      #
      # Apply the formulas and add a mass_g column
       puffins <- puffins %>%
         rowwise() %>%
         mutate(row_mass_g = row_mass_g(species, row_prey_length_cm, prey_number)) %>%
         ungroup()
      #
      #Now we want the total mass of each prey load (with combined species)
       total_load_mass <- puffins %>%
         group_by(ID, day, year) %>%
         summarize(
           mass_g = sum(row_mass_g),
           species_combined = paste(unique(species), collapse = ", ")
         ) %>%
         ungroup()
      #
      #Remove columns we don't need
       total_load_mass <- total_load_mass [,c(1,4)]
      #
      #Add load mass data from total_load_mass to puffins
       puffins <- left_join(puffins, total_load_mass, "ID")
      #Rename newly input mass_g to load_mass_g
       colnames(puffins)[colnames(puffins) == 'mass_g'] <- 'load_mass_g'
       
      #Number of prey items per prey load calculations
      #First we sum prey_number grouped by observation ID
       prey_load_size <- puffins %>%
         group_by(ID, day,cam, year) %>%
         summarize(
           total_load_size = sum(prey_number),
           species_combined = paste(unique(species), collapse = ", ")
         ) %>%
         ungroup()
      #
      #Remove column we don't need
       prey_load_size <- prey_load_size [,c(1,5)]
      #
      #Add prey number data from prey_load_size to puffins
       puffins <- left_join(puffins, prey_load_size, "ID")
       
       
      #Prey load energy calculations
      #Apply the formulas to each species
       energy_kJ <- function(species, row_prey_length_cm, prey_number) {
         if (species == "Sandlance") {
           energy <- ifelse(!is.na(row_prey_length_cm), (((row_prey_length_cm * 10) ^ 2.11) * (209 * 10^-6)) * 5.06 * prey_number, 0)
         } else if (species == "Capelin") {
           energy <- ifelse(!is.na(row_prey_length_cm), (((row_prey_length_cm * 10) ^ 2.76) * (2.7 * 10^-5)) * 4.90 * prey_number, 0)
         } else if (species %in% c("Unknown", "LowConfidence", "Fish larva")) {
           energy <- ifelse(!is.na(row_prey_length_cm), (((row_prey_length_cm * 10) ^ 2.41) * (1.18 * 10^-4)) * 4.98 * prey_number, 0)
         } else {
           energy <- NA
         }
         return(energy)
       }
      #
      #Apply the formulas and update the prey_length_px column
       puffins <- puffins %>%
         rowwise() %>%
         mutate(prey_energy_kJ = energy_kJ(species, row_prey_length_cm, prey_number)) %>%
         ungroup()
      #
      #Rename "0" values
       puffins$prey_energy_kJ[puffins$prey_energy_kJ == 0] <- NA
      #
      #Calculate the energy per load (not energy per observation per prey species)
       prey_load_energy <- puffins %>%
         group_by(ID, day, year) %>%
         summarize(
           total_energy_content = sum(prey_energy_kJ),
           species_combined = paste(unique(species), collapse = ", ")
         ) %>%
         ungroup()
      #
      #Remove columns we don't need
       prey_load_energy <- prey_load_energy [,c(1,4)]
      #
      #Add load energy data from prey_load_energy to puffins
       puffins <- left_join(puffins, prey_load_energy, "ID")
       
       
###Add a chick-rearing period column
       puffins <- puffins %>%
         mutate(period = case_when(
           day <= 4 ~ "1-4",
           day <= 8 ~ "5-8",
           day <= 12 ~ "9-12",
           day <= 16 ~ "13-16",
           day > 16 ~ ">16"
         ))
      #Check levels and rearrange 
       table(puffins$period)
       puffins$period <- factor(puffins$period, levels=c('1-4', '5-8', '9-12', '13-16','>16'))
       
      #Set relevant variables to factor
       puffins$species <- factor(puffins$species)
       puffins$year <- factor(puffins$year)
       puffins$period <- factor(puffins$period)
       
      #Create a data frame with grouped observation data not separated by species (avoids duplicate total_load_energy and total_load_size)
      #Pool prey species IDs and remove duplicate observation IDs
       grouped_puffins <- puffins %>%
         group_by(ID, bird_ID, day, period, cam, year, total_energy_content, total_load_size) %>%
         summarize(
           species_combined = paste(unique(species), collapse = ", ")
         ) %>%
         ungroup()
       
      #Remove all excess environment elements (including GPS data frame)
       rm(list=setdiff(ls(), c("puffins", "grouped_puffins", "daily_trips")))
       
       table(grouped_puffins$period)
       
### GPS + Photo observation joined #############################################
       
      #Calculate average trips per day
       avg_trip_data <- daily_trips %>%
         group_by(year, day) %>%
         summarise(avg_trips_per_day = mean(total_trips))
       
      #Make a year_day identifier column in avg_trip_data
       avg_trip_data$identifier <- paste(avg_trip_data$year, avg_trip_data$day, sep = "_")
      #Remove day and year columns
       avg_trip_data <- subset(avg_trip_data, select = -c(day, year))
       
      #Make a year_day identifier column in avg_daily_energy
       grouped_puffins$identifier <- paste(grouped_puffins$year, grouped_puffins$day, sep = "_")
       
      #Join with trips data frame using the identifier columns
       joined_data <- grouped_puffins %>%
         left_join(avg_trip_data, by = "identifier")
       
      #Remove late 2023 (and other years that go beyond last GPS measure) until I figure out what to do with it   
       joined_data <- joined_data %>% filter(!identifier %in% c("2023_10","2023_11", "2023_12", "2023_13", "2023_14", "2023_15", "2023_16", "2023_17", "2023_18", "2023_19", "2023_20",
                                                                "2024_41"))
      #Apply formula to calculate total daily energy delivered to chicks on average
       joined_data$total_daily_energy <- with (joined_data, total_energy_content * avg_trips_per_day * 2 ) #times 2 for each parent
       
      #Remove rows for days no energy and no trips
       joined_data <- joined_data[joined_data$total_daily_energy != 0, ]
       joined_data <- joined_data[joined_data$avg_trips_per_day != 0, ]
       
### Statistical test ###########################################################
      
    #Log transform data
      shapiro.test(joined_data$total_daily_energy)
      joined_data_log <-joined_data
      joined_data_log$total_daily_energy_log <- log(joined_data$total_daily_energy + 1)
      shapiro.test(joined_data_log$total_daily_energy_log)
      #
    #Run the linear model 
      total_energy_model <- lm(total_daily_energy_log ~ period + year, data =joined_data_log)
      summary(total_energy_model)
      plot(total_energy_model)
    #Check overdispersion
      library(DHARMa)
      testDispersion(total_energy_model) #Dispersion = 0.97888, model not overdispersed
    #Overall results - Type II Wald Test to get p values 
      library(car) 
      Anova(total_energy_model, type = "II")
    #Pairwise comparisons to check results
      em_total_energy <- emmeans(total_energy_model, ~ period)
      contrast(em_total_energy, method = "pairwise")
      #Model variance
      performance::r2(total_energy_model)
      