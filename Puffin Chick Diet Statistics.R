#Statistical tests for chick diet metrics published in: CHANGES ON THE MENU: INTRA-ANNUAL SHIFTS IN DIET AND FOOD PROVISIONING OF ATLANTIC PUFFIN (FRATERCULA ARCTICA) CHICKS
#Authored by: AIMEE FRASER, CHRISTINA PETALAS, RAPHAËL A. LAVOIE, KYLE H. ELLIOTT
#2026

####Format and clean data#######################################################

  rm(list=ls())
  
  library(dplyr)
  library(glmmTMB) #GLMMs -> species_models and number_model
  library(lme4) #LMM -> all other diet metrics analysed
  library(emmeans) #for pairwise comparisons
  library(stringr)

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
  
  #Clean data
  puffins_cleaned <- puffins_raw
  #change species to unknown when confidence is below 50%
  levels(puffins_cleaned$species) <- c(levels(puffins_cleaned$species),"LowConfidence")
  puffins_cleaned <- puffins_cleaned %>%
    mutate(species = if_else(confidence < 50, "Unknown", species)) #Renames LowConfidence to Unknown
  #remove rows with flying birds (where prey are more difficult to analyze)
  puffins_cleaned <- subset(puffins_cleaned, position  == "grounded")
  #remove other island data (we only want to use data from Betchouanes)
  puffins_cleaned <- subset(puffins_cleaned, island  == "Betchouanes")
  #set data frames to either puffins_raw or puffins_cleaned
  puffins <- puffins_cleaned
  rm(puffins_raw, puffins_cleaned)


###Formulas and Calculations for diet metrics

  #Add prey_number = 1 for all NAs
  puffins <- puffins %>% mutate(prey_number = ifelse(is.na(prey_number), 1, prey_number))
  
#Calculation for prey species frequency of occurrence
  #Count totals of each species to know what proportions to use for each species when calculating the length modifiers of observations with 'unknown' species identifications.
  species_totals <- puffins %>%
    group_by(species, year) %>%
    summarise(total_count = sum(prey_number, na.rm = TRUE))

  #Unknown Spp Formulas
  #Weighted mean annual constants to calculate 'unknown' prey lengths from prey tail-end measurements based on annual species abundance
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

#Prey species length formulas
  #Species-specific formulas for each measurement type
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
  # Apply the formulas and update the prey_length_px column
  puffins <- puffins %>%
    rowwise() %>%
    mutate(prey_length_px = species_formula(species, year, width, tail_end, head_length, length_est_px)) %>%
    ungroup()
  #
  #Change "0" values to NA
  puffins$prey_length_px[puffins$prey_length_px == 0] <- NA
  #
  #Convert pixel prey lengths to cm
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
  #
  #Remove NAs
  puffins <- subset(puffins, !is.na(row_prey_length_cm))
  
#Prey mass calculations
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
  #Now calculate mass for the prey load as a whole (accounting for prey_number)
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
  #Apply the formulas and add a mass_g column
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
  #remove columns we don't need
  total_load_mass <- total_load_mass [,c(1,4)]
  #
  #Add total load mass data from total_load_mass to puffins data frame
  puffins <- left_join(puffins, total_load_mass, "ID")
  #rename newly input mass_g to load_mass_g
  colnames(puffins)[colnames(puffins) == 'mass_g'] <- 'load_mass_g'


#Number of prey items per prey load calculations (prey quantity)
  #Sum prey_number grouped by observation ID
  prey_load_size <- puffins %>%
    group_by(ID, day,cam, year) %>%
    summarize(
      total_load_size = sum(prey_number),
      species_combined = paste(unique(species), collapse = ", ")
    ) %>%
    ungroup()
  #
  #Remove columns we don't need
  prey_load_size <- prey_load_size [,c(1,5)]
  #
  #Add prey number data from prey_load_size to puffins
  puffins <- left_join(puffins, prey_load_size, "ID")


#Prey load energy calculations
  #Apply the formulas to each species with species-specific constants
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
  #Rename "0" values to NA
  puffins$prey_energy_kJ[puffins$prey_energy_kJ == 0] <- NA
  #
  #Calculate the energy per prey load (not energy per observation per prey species)
  prey_load_energy <- puffins %>%
    group_by(ID, day, year) %>%
    summarize(
      total_energy_content = sum(prey_energy_kJ),
      species_combined = paste(unique(species), collapse = ", ")
    ) %>%
    ungroup()
  #
  #remove columns we don't need
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
  puffins$period <- factor(puffins$period, levels=c('1-4', '5-8', '9-12', '13-16',">16"))
  table(puffins$period)
  
  #remove all excess environment elements
  rm(list=setdiff(ls(), c("puffins")))
  
  #Set relevant variables to factor
  puffins$species <- factor(puffins$species)
  puffins$year <- factor(puffins$year)
  puffins$period <- factor(puffins$period)

#Create a data frame with grouped observation data not separated by species (avoids duplicate total_load_energy and total_load_size)
  #pool prey species IDs and remove duplicate observation IDs
  grouped_puffins <- puffins %>%
    group_by(ID, day, bird_ID, period, cam, year, total_energy_content, total_load_size) %>%
    summarize(
      species_combined = paste(unique(species), collapse = ", ")
    ) %>%
    ungroup()
  
  #Set day as factor
  grouped_puffins$day <- factor(grouped_puffins$day)
  
  table(grouped_puffins$period)
  
  
  
###Statistical Analyses#########################################################    

#Prey Species analysis
  #We go species by species here because the GLMM family = binomial (multinomial cannot be used) so it will use information on if a species is there or not there
  #First set up the data to have presence of sandlace, capelin and larvae indicated, then group prey loads (so R doesnt log one observation with capelin and larvae as 2 observations without sandlance)
  grouped_puffins <- grouped_puffins %>%
    mutate(capelin_presence = if_else(str_detect(species_combined, "Capelin"), 1, 0))
  #
  grouped_puffins <- grouped_puffins %>%
    mutate(sandlance_presence = if_else(str_detect(species_combined, "Sandlance"), 1, 0))
  #
  grouped_puffins <- grouped_puffins %>%
    mutate(larva_presence = if_else(str_detect(species_combined, "Fish larva"), 1, 0))
  
  ###Capelin
    #Run the model
      species_capelin_model <- glmer(capelin_presence ~ period + year + (1 | bird_ID), 
                   data = grouped_puffins, family = binomial)
      summary(species_capelin_model)
    #Check overdispersion
      library(DHARMa)
      testDispersion(species_capelin_model) #Dispersion = 1.0706, model not overdispersed
    #plot to check residual interdependence
      res <- simulateResiduals(species_capelin_model)
      plot(res) #residuals look good
    #Overall results - Type II Wald Test to get p values 
      library(car) 
      Anova(species_capelin_model, type = "II")
    #Pairwise comparisons to check results
      em_capelin <- emmeans(species_capelin_model, ~ period)
      contrast(em_capelin, method = "pairwise")
    #Model variance
      VarCorr(species_capelin_model)
      performance::r2(species_capelin_model)
  
  ###Sandlance
    #Run the model
      species_sandlance_model <- glm(sandlance_presence ~ period + year, 
                                     data = grouped_puffins, family = binomial)
      summary(species_sandlance_model)
    #Check overdispersion
      library(DHARMa)
      testDispersion(species_sandlance_model) #Dispersion = 1.0033, model not overdispersed
    #plot to check residual interdependence
      res <- simulateResiduals(species_sandlance_model)
      plot(res) #looks good
    #Overall results - Type II Wald Test to get p values 
      library(car) 
      Anova(species_sandlance_model, type = "II")
    #Pairwise comparisons to check results
      em_sandlance <- emmeans(species_sandlance_model, ~ period)
      contrast(em_sandlance, method = "pairwise")
    #model variance
      performance::r2(species_sandlance_model)
    
  ###Larvae
    #Run the model
      species_larva_model <- glm(larva_presence ~ period + year, 
                                       data = grouped_puffins, family = binomial)
       summary(species_larva_model)
     #Check overdispersion
       library(DHARMa)
       testDispersion(species_larva_model) #Dispersion = 1.0166, model not overdispersed
     #Plot to check residual interdependence
       res <- simulateResiduals(species_larva_model)
       plot(res) #only one quantile in the red, with no overdispersion detected earlier, this should be ok.
    #Overall results - Type II Wald Test to get p values 
      library(car) 
      Anova(species_larva_model, type = "II")
    #Pairwise comparisons to check results
      em_larva <- emmeans(species_larva_model, ~ period)
      contrast(em_larva, method = "pairwise")
    #model variance
      performance::r2(species_larva_model)
      
###Number of prey items analysis (prey quantity)
  #Fit to a GLMM model (count data not continuous for the response variable) with family = poisson for count
    number_model <- glmer(total_load_size ~ period + year + (1 | bird_ID), 
                   data = grouped_puffins, family = poisson)
    summary(number_model)
  #Check overdispersion
    library(DHARMa)
    testDispersion(number_model) #Dispersion = 1.0083, model not overdispersed
  #Plot to check residual interdependence
    res <- simulateResiduals(number_model)
    plot(res)#only one quantile in the red, with no overdispersion detected earlier, this should be ok.
  #Overall results - Type II Wald Test to get p values 
    library(car) 
    Anova(number_model, type = "II")
  #Pairwise comparisons to check results
    em_number <- emmeans(number_model, ~ period)
    contrast(em_number, method = "pairwise")
  #Model variance
    VarCorr(number_model)
    performance::r2(number_model)
    
###Mass analysis (log-transformed)
  #Log-transform data
    shapiro.test(puffins$longest_mass_g)
    puffins_log <-puffins
    puffins_log$longest_mass_g_log <- log(puffins$longest_mass_g + 1)
    shapiro.test(puffins_log$longest_mass_g_log)
  #Fit to a LM model (log transformed)
    mass_model <- lm(longest_mass_g_log ~ period + year, data = puffins_log)
    summary(mass_model)
    AIC(mass_model)
    plot(mass_model) #looks good
  #Overall results - Type II Wald Test to get p values 
    library(car) 
    Anova(mass_model, type = "II")
  #Pairwise comparisons to check results
    em_mass <- emmeans(mass_model, ~ period)
    contrast(em_mass, method = "pairwise")
    #Model variance
    performance::r2(mass_model)
    
###Energy analysis (log transformed)
  #Log-transform data
    shapiro.test(grouped_puffins$total_energy_content)
    grouped_puffins_log <-grouped_puffins
    grouped_puffins_log$total_energy_content_log <- log(grouped_puffins$total_energy_content + 1)
    shapiro.test(grouped_puffins_log$total_energy_content_log)
  #Fit to a LM model
    energy_model <- lm(total_energy_content_log ~ period + year, data = grouped_puffins_log)
    AIC(energy_model)
    plot(energy_model) #looks good
  #Overall results - Type II Wald Test to get p values 
    library(car) 
    Anova(energy_model, type = "II")
  #Pairwise comparisons to check results
    em_energy <- emmeans(energy_model, ~ period)
    contrast(em_energy, method = "pairwise")
  #Model variance
    performance::r2(energy_model)
  

###Test to check effect of burrow ID and if it should be included in the models above
  
    #Prey types
      #Capelin
          #Model without individual ID
          m1 <- glm(capelin_presence ~ period + year, 
                                         data = grouped_puffins, family = binomial)
          #
          #Model with individual ID
          m2 <- glmer(capelin_presence ~ period + year + (1 | bird_ID), 
                                         data = grouped_puffins, family = binomial)
          #Compare and extract variance
          AIC(m1, m2) 
          VarCorr(m2)
          #bird_ID is significant and helps the model - use it
          
      #Sandlance
          #Model without individual ID
          m1 <- glm(sandlance_presence ~ period + year, 
                      data = grouped_puffins, family = binomial)
          #
          #Model with individual ID
          m2 <- glmer(sandlance_presence ~ period + year + (1 | bird_ID), 
                      data = grouped_puffins, family = binomial)
          #Compare and extract variance
          AIC(m1, m2)
        #bird_ID is not significant 
          
      #Larvae
          #Model without individual ID
          m1 <- glm(larva_presence ~ period + year, 
                      data = grouped_puffins, family = binomial)
          #
          #Model with individual ID
          m2 <- glmer(larva_presence ~ period + year + (1 | bird_ID), 
                      data = grouped_puffins, family = binomial)
          #Compare and extract variance
          AIC(m1, m2)
        #bird_ID is not significant 
          
          
  #Number of prey items per load (prey quantity)
    #Model without individual ID
    m1 <- glm(total_load_size ~ period + year, 
                data = grouped_puffins, family = poisson)
    #
    #Model with individual ID
    m2 <- glmer(total_load_size ~ period + year + (1 | bird_ID), 
                data = grouped_puffins, family = poisson)
    #Compare and extract variance
    AIC(m1, m2)
    VarCorr(m2)
    library(performance)
    icc(m1)
    icc(m2)
  #bird_ID is significant and improves the model - use it
    
  #Prey Mass
    #Model without individual ID
    m1 <- lm(longest_mass_g_log ~ period + year, data = puffins_log)
    #
    #Model with individual ID
    m2 <- lmer(longest_mass_g_log ~ period + year + (1 | bird_ID), data = puffins_log)
    #
    #Compare and extract variance
    AIC(m1, m2) 
    VarCorr(m2)
    library(performance)
    icc(m1) 
    icc(m2)
  #bird_ID is significant but harms the model - do not use
    
  #Energy of prey loads
    #Model without individual ID
    m1 <- lm(total_energy_content_log ~ period + year, data = grouped_puffins_log)
    #
    #Model with individual ID
    m2 <- lmer(total_energy_content_log ~ period + year + (1 | bird_ID), data = grouped_puffins_log)
    #
    #Compare and extract variance
    AIC(m1, m2)
    VarCorr(m2)
    library(performance)
    icc(m1)
    icc(m2)
  #bird_ID is significant but makes the model worse - do not use

  

###Mean values##################################################################

 #Prey quantity
    #Means
    prey_number_mean <- grouped_puffins %>%
      group_by(period) %>%
      summarise(Mean = mean(total_load_size, na.rm = TRUE))
    print(prey_number_mean)
    #Stdev
    prey_number_sd <- grouped_puffins %>%
      group_by(period) %>%
      summarise(Dev = sd(total_load_size, na.rm = TRUE))
    print(prey_number_sd)
  
  #Prey mass
    #Means
    mass_mean <- puffins %>%
      group_by(period) %>%
      summarise(Mean = mean(longest_mass_g, na.rm = TRUE))
    print(mass_mean)
    #Stdev
    mass_sd <- puffins %>%
      group_by(period) %>%
      summarise(Dev = sd(longest_mass_g, na.rm = TRUE))
    print(mass_sd)

    
   #Energy of prey loads
    #means
    energy_mean <- grouped_puffins %>%
      group_by(period) %>%
      summarise(Mean = mean(total_energy_content, na.rm = TRUE))
    print(as.data.frame(energy_mean), digits = 4)
    #stdev
    energy_sd <- grouped_puffins %>%
      group_by(period) %>%
      summarise(Dev = sd(total_energy_content, na.rm = TRUE))
    print(as.data.frame(energy_sd), digits = 4)


###Sample size by groups
    #Total
    table(grouped_puffins$period)

    #Subset year
     puffins_2021 <- subset(grouped_puffins, year  == "2021")
     puffins_2023 <- subset(grouped_puffins, year  == "2023")
     puffins_2024 <- subset(grouped_puffins, year  == "2024")
     
     table(puffins_2021$period)
     table(puffins_2023$period)
     table(puffins_2024$period)
     
    #Frequency of Occurrence and Numerical Abundance values
     puffins_1to4 <- subset(puffins, period  == "1-4")
     puffins_5to8 <- subset(puffins, period  == "5-8")
     puffins_9to12 <- subset(puffins, period  == "9-12")
     puffins_13to16 <- subset(puffins, period  == "13-16")
     puffins_16plus <- subset(puffins, period  == ">16")
     
     #FO
     table(puffins_1to4$species)
     table(puffins_5to8$species)
     table(puffins_9to12$species)
     table(puffins_13to16$species)
     table(puffins_16plus$species)
     
     #Numerical Abundance
     puffins_numerical_abundance <- puffins %>%
       group_by(species, period) %>%
       summarise(total_caught = sum(prey_number), .groups = "drop")
     
     #View the result
     print(puffins_numerical_abundance)
     
