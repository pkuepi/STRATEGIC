#################################################################
#                                                               #
#  Strategies to inTerrupt RAbies Transmission for the          #
#  Elimination Goal by 2030 In China (STRATEGIC) study          #
#                                                               #        
#  8.plot_function                                              #
#  This file is used for output of the plots                    #
#  Modified from the WHO rabies modelling consortium study      #
#  URL: https://github.com/katiehampson1978/rabies_PEP_access   #
#  Ref.: The Lancet Infectious Diseases 2019; 19(1): 102-11.    #
#                                                               #
#################################################################



########################################################
# 1.Functions to prepare senario names for result plot #
########################################################

####rename senarios of different strategies

scenario_prep <- function(dataframe){
  
  dataframe$scenario[dataframe$scenario=="a1"] <- "1"
  dataframe$scenario[dataframe$scenario=="a2"] <- "2"
  dataframe$scenario[dataframe$scenario=="a3"] <- "3"
  dataframe$scenario[dataframe$scenario=="a4"] <- "4"
  dataframe$scenario[dataframe$scenario=="a5"] <- "5"
  dataframe$scenario[dataframe$scenario=="a6"] <- "6"
  dataframe$scenario[dataframe$scenario=="a7"] <- "7"
  dataframe$scenario[dataframe$scenario=="a8"] <- "8"
  return(dataframe)
}

####fix metric_name  

fix_metric_name <- function(dataframe){
  
  # Read metric column as caharacter
  dataframe$metric <- as.character(dataframe$metric)
  
  # Improve names visually for plot
  dataframe$metric[dataframe$metric=="cost_per_death_averted"] <- "Cost per death averted"
  dataframe$metric[dataframe$metric=="cost_per_YLL_averted"] <- "Cost per DALY averted"
  dataframe$metric[dataframe$metric=="cost_per_DALY_averted"] <- "Cost per DALY averted"
  dataframe$metric[dataframe$metric=="deaths_averted_per_100k_vaccinated"] <- "Deaths averted per 100k vaccinated"
  dataframe$metric[dataframe$metric=="fully_vaccinated"] <- "Fully vaccinated"
  dataframe$metric[dataframe$metric=="total_cost"] <- "Total cost"
  dataframe$metric[dataframe$metric=="total_deaths"] <- "Total deaths"
  dataframe$metric[dataframe$metric=="total_deaths_averted"] <- "Total deaths averted"
  dataframe$metric[dataframe$metric=="total_vials"] <- "Total vials"
  dataframe$metric[dataframe$metric=="total_YLL"] <- "Total DALY"
  dataframe$metric[dataframe$metric=="total_DALY"] <- "Total DALY"
  dataframe$metric[dataframe$metric=="total_YLL_averted"] <- "Total DALY averted"
  dataframe$metric[dataframe$metric=="total_DALY_averted"] <- "Total DALY averted"
  dataframe$metric[dataframe$metric=="vaccinated"] <- "Vaccinated"
  
  return(dataframe)
}
####arrange_factor_levels  
arrange_factor_levels <- function(dataframe, type){
  
  if(type=="multi"){
    levels(dataframe$metric)[levels(dataframe$metric)=="total_deaths"] <- "A) Total deaths"
    levels(dataframe$metric)[levels(dataframe$metric)=="total_YLL"] <- "B) Total DALYs (x1000)"
    levels(dataframe$metric)[levels(dataframe$metric)=="cost_per_death_averted"] <- "C) Cost per death averted"
    levels(dataframe$metric)[levels(dataframe$metric)=="cost_per_YLL_averted"] <- "D) Cost per DALY averted"
    
    new_order <- c("A) Total deaths","B) Total DALYs (x1000)",  "C) Cost per death averted", "D) Cost per DALY averted")
  } else if(type=="uni"){
    levels(dataframe$metric)[levels(dataframe$metric)=="total_deaths"] <- "A) Total deaths (x1000)"
    levels(dataframe$metric)[levels(dataframe$metric)=="total_YLL"] <- "B) Total DALYs (x1000)"
    levels(dataframe$metric)[levels(dataframe$metric)=="cost_per_death_averted"] <- "C) Cost per death averted"
    levels(dataframe$metric)[levels(dataframe$metric)=="cost_per_YLL_averted"] <- "D) Cost per DALY averted"
    
    new_order <- c("A) Total deaths (x1000)","B) Total DALYs (x1000)",  "C) Cost per death averted", "D) Cost per DALY averted")
  }
  
  dataframe$metric <- as.factor(dataframe$metric)
  dataframe <- arrange(transform(dataframe,
                                 metric=factor(metric, levels=new_order)), metric)
  
  return(dataframe)
}

####################################################
#2.Functions to summarise univariate data for plot #
####################################################

summarise_univar_data <- function(horizon_dataframe, lci, uci){
  
  horizon_dataframe$scenario <- as.character(horizon_dataframe$scenario)
  horizon_dataframe$scenario[horizon_dataframe$scenario=="a2"] <- "Scenario 2 \n Improved PEP"
  horizon_dataframe$scenario[horizon_dataframe$scenario=="a3"] <- "Scenario 3a \n Dog vax"
  horizon_dataframe$scenario[horizon_dataframe$scenario=="a7"] <- "Scenario 4c \n Dog vax+IBCM"
  
  #' * TOTAL_DEATHS *
  summed <- horizon_dataframe %>%
    group_by(scenario, variable) %>%
    summarise(m = mean(total_deaths, na.rm=TRUE),
              lci=mean(total_deaths_lci, na.rm=TRUE),
              uci=mean(total_deaths_uci, na.rm=TRUE),
              metric="total_deaths")
  final_df = summed
  
  #' * TOTAL_VIALS *
  summed <- horizon_dataframe %>%
    group_by(scenario, variable) %>%
    summarise(m = mean(total_vials, na.rm=TRUE),
              lci=mean(total_vials_lci, na.rm=TRUE),
              uci=mean(total_vials_uci, na.rm=TRUE),
              metric="total_vials")
  final_df <- rbind(final_df, summed)
  
  #' * TOTAL_COST *
  summed <- horizon_dataframe %>%
    group_by(scenario, variable) %>%
    summarise(m = mean(total_cost, na.rm=TRUE),
              lci=mean(total_cost_lci, na.rm=TRUE),
              uci=mean(total_cost_uci, na.rm=TRUE),
              metric="total_cost")
  final_df <- rbind(final_df, summed)
  
  #' * TOTAL_DEATHS_AVERTED *
  summed <- horizon_dataframe %>%
    group_by(scenario, variable) %>%
    summarise(m = mean(total_deaths_averted, na.rm=TRUE),
              lci=mean(total_deaths_averted_lci, na.rm=TRUE),
              uci=mean(total_deaths_averted_uci, na.rm=TRUE),
              metric="total_deaths_averted")
  final_df <- rbind(final_df, summed)
  
  #' * TOTAL_YLL *
  summed <- horizon_dataframe %>%
    group_by(scenario, variable) %>%
    summarise(m = mean(total_YLL, na.rm=TRUE),
              lci=mean(total_YLL_lci, na.rm=TRUE),
              uci=mean(total_YLL_uci, na.rm=TRUE),
              metric="total_YLL")
  final_df <- rbind(final_df, summed)
  
  #' * TOTAL_YLL_AVERTED *
  summed <- horizon_dataframe %>%
    group_by(scenario, variable) %>%
    summarise(m = mean(total_YLL_averted, na.rm=TRUE),
              lci=mean(total_YLL_averted_lci, na.rm=TRUE),
              uci=mean(total_YLL_averted_uci, na.rm=TRUE),
              metric="total_YLL_averted")
  final_df <- rbind(final_df, summed)

  #' * VACCINATED *
  summed <- horizon_dataframe %>%
    group_by(scenario, variable) %>%
    summarise(m = mean(vaccinated, na.rm=TRUE),
              lci=mean(vaccinated_lci, na.rm=TRUE),
              uci=mean(vaccinated_uci, na.rm=TRUE),
              metric="vaccinated")
  final_df <- rbind(final_df, summed)
  
  #' * FULLY_VACCINATED *
  summed <- horizon_dataframe %>%
    group_by(scenario, variable) %>%
    summarise(m = mean(fully_vaccinated, na.rm=TRUE),
              lci=mean(fully_vaccinated_lci, na.rm=TRUE),
              uci=mean(fully_vaccinated_uci, na.rm=TRUE),
              metric="fully_vaccinated")
  final_df <- rbind(final_df, summed)
  
  #' * COST_PER_DEATH_AVERTED *
  summed <- horizon_dataframe %>%
    group_by(scenario, variable) %>%
    summarise(m = mean(cost_per_death_averted, na.rm=TRUE),
              lci=mean(cost_per_death_averted_lci, na.rm=TRUE),
              uci=mean(cost_per_death_averted_uci, na.rm=TRUE),
              metric="cost_per_death_averted")
  final_df <- rbind(final_df, summed)
  
  #' * COST_PER_YLL_AVERTED *
  summed <- horizon_dataframe %>%
    group_by(scenario, variable) %>%
    summarise(m = mean(cost_per_YLL_averted, na.rm=TRUE),
              lci=mean(cost_per_YLL_averted_lci, na.rm=TRUE),
              uci=mean(cost_per_YLL_averted_uci, na.rm=TRUE),
              metric="cost_per_YLL_averted")
  final_df <- rbind(final_df, summed)
  
  #' * DEATHS_AVERTED_PER_100K_VACCINATED *
  summed <- horizon_dataframe %>%
    group_by(scenario, variable) %>%
    summarise(m = mean(deaths_averted_per_100k_vaccinated, na.rm=TRUE),
              lci=mean(deaths_averted_per_100k_vaccinated_lci, na.rm=TRUE),
              uci=mean(deaths_averted_per_100k_vaccinated_uci, na.rm=TRUE),
              metric="deaths_averted_per_100k_vaccinated")
  final_df <- rbind(final_df, summed)
  
  final_df$metric <- as.character(final_df$metric)
  final_df$metric <- as.factor(final_df$metric)
  
  return(final_df)
  
}

######################################################
#3.Functions to summarise multivariate data for plot #
######################################################

summarise_multivariate_data <- function(horizon_dataframe, setting, lci, uci){
  
  horizon_dataframe$scenario <- as.character(horizon_dataframe$scenario)
  horizon_dataframe$scenario[horizon_dataframe$scenario=="a1"] <- "1"
  horizon_dataframe$scenario[horizon_dataframe$scenario=="a2"] <- "2"
  horizon_dataframe$scenario[horizon_dataframe$scenario=="a3"] <- "3"
  horizon_dataframe$scenario[horizon_dataframe$scenario=="a4"] <- "4"
  horizon_dataframe$scenario[horizon_dataframe$scenario=="a5"] <- "5"
  horizon_dataframe$scenario[horizon_dataframe$scenario=="a6"] <- "6"
  horizon_dataframe$scenario[horizon_dataframe$scenario=="a7"] <- "7"
  horizon_dataframe$scenario[horizon_dataframe$scenario=="a8"] <- "8"
  
  # Take only required scenarios (above)
  req_scenarios <- c("1", "2", "3", "4", "5", "6",'7','8')
  horizon_dataframe <- horizon_dataframe[which(horizon_dataframe$scenario %in% req_scenarios),]
  sort(unique(horizon_dataframe$scenario))
  
  # Set up arguments for dplyr summaries
  if(setting == "global"){arg_setting = list(quo(scenario))}
  
  #' * TOTAL_DEATHS *
  summed <- horizon_dataframe %>%
    group_by(!!! arg_setting) %>%
    summarise(m = mean(total_deaths, na.rm=TRUE),
              lci=mean(total_deaths_lci, na.rm=TRUE),
              uci=mean(total_deaths_uci, na.rm=TRUE),
              metric="total_deaths")
  final_df = summed
  
  #' * TOTAL_VIALS *
  summed <- horizon_dataframe %>%
    group_by(!!! arg_setting) %>%
    summarise(m = mean(total_vials, na.rm=TRUE),
              lci=mean(total_vials_lci, na.rm=TRUE),
              uci=mean(total_vials_uci, na.rm=TRUE),
              metric="total_vials")
  final_df <- rbind(final_df, summed)
  
  #' * TOTAL_COST *
  summed <- horizon_dataframe %>%
    group_by(!!! arg_setting) %>%
    summarise(m = mean(total_cost, na.rm=TRUE),
              lci=mean(total_cost_lci, na.rm=TRUE),
              uci=mean(total_cost_uci, na.rm=TRUE),
              metric="total_cost")
  final_df <- rbind(final_df, summed)
  
  #' * TOTAL_DEATHS_AVERTED *
  summed <- horizon_dataframe %>%
    group_by(!!! arg_setting) %>%
    summarise(m = mean(total_deaths_averted, na.rm=TRUE),
              lci=mean(total_deaths_averted_lci, na.rm=TRUE),
              uci=mean(total_deaths_averted_uci, na.rm=TRUE),
              metric="total_deaths_averted")
  final_df <- rbind(final_df, summed)
  
  #' * TOTAL_YLL *
  summed <- horizon_dataframe %>%
    group_by(!!! arg_setting) %>%
    summarise(m = mean(total_YLL, na.rm=TRUE),
              lci=mean(total_YLL_lci, na.rm=TRUE),
              uci=mean(total_YLL_uci, na.rm=TRUE),
              metric="total_YLL")
  final_df <- rbind(final_df, summed)
  
  #' * TOTAL_YLL_AVERTED *
  summed <- horizon_dataframe %>%
    group_by(!!! arg_setting) %>%
    summarise(m = mean(total_YLL_averted, na.rm=TRUE),
              lci=mean(total_YLL_averted_lci, na.rm=TRUE),
              uci=mean(total_YLL_averted_uci, na.rm=TRUE),
              metric="total_YLL_averted")
  final_df <- rbind(final_df, summed)

  #' * VACCINATED *
  summed <- horizon_dataframe %>%
    group_by(!!! arg_setting) %>%
    summarise(m = mean(vaccinated, na.rm=TRUE),
              lci=mean(vaccinated_lci, na.rm=TRUE),
              uci=mean(vaccinated_uci, na.rm=TRUE),
              metric="vaccinated")
  final_df <- rbind(final_df, summed)
  
  #' * FULLY_VACCINATED *
  summed <- horizon_dataframe %>%
    group_by(!!! arg_setting) %>%
    summarise(m = mean(fully_vaccinated, na.rm=TRUE),
              lci=mean(fully_vaccinated_lci, na.rm=TRUE),
              uci=mean(fully_vaccinated_uci, na.rm=TRUE),
              metric="fully_vaccinated")
  final_df <- rbind(final_df, summed)
  
  #' * COST_PER_DEATH_AVERTED *
  summed <- horizon_dataframe %>%
    group_by(!!! arg_setting) %>%
    summarise(m = mean(cost_per_death_averted, na.rm=TRUE),
              lci=mean(cost_per_death_averted_lci, na.rm=TRUE),
              uci=mean(cost_per_death_averted_uci, na.rm=TRUE),
              metric="cost_per_death_averted")
  final_df <- rbind(final_df, summed)
  
  #' * COST_PER_YLL_AVERTED *
  summed <- horizon_dataframe %>%
    group_by(!!! arg_setting) %>%
    summarise(m = mean(cost_per_YLL_averted, na.rm=TRUE),
              lci=mean(cost_per_YLL_averted_lci, na.rm=TRUE),
              uci=mean(cost_per_YLL_averted_uci, na.rm=TRUE),
              metric="cost_per_YLL_averted")
  final_df <- rbind(final_df, summed)
  
  #' * DEATHS_AVERTED_PER_100K_VACCINATED *
  summed <- horizon_dataframe %>%
    group_by(!!! arg_setting) %>%
    summarise(m = mean(deaths_averted_per_100k_vaccinated, na.rm=TRUE),
              lci=mean(deaths_averted_per_100k_vaccinated_lci, na.rm=TRUE),
              uci=mean(deaths_averted_per_100k_vaccinated_uci, na.rm=TRUE),
              metric="deaths_averted_per_100k_vaccinated")
  final_df <- rbind(final_df, summed)
  
  final_df$metric <- as.character(final_df$metric)
  final_df$metric <- as.factor(final_df$metric)
  
  return(final_df)
}

###################################################
#4.Functions to summarise data for projection plot#
###################################################

projected_outcomes_plot_prep <- function(region_df){
  
  # Prepare region dataframes
  region_deaths <- data.frame(Level=region_df$region, year=region_df$year, scenario=region_df$scenario,
                               variable="Deaths", mean=region_df$total_deaths,
                               lci=region_df$total_deaths_lci, uci=region_df$total_deaths_uci, stringsAsFactors=FALSE)
  region_vaccinated <- data.frame(Level=region_df$region, year=region_df$year, scenario=region_df$scenario,
                                   variable="Persons vaccinated", mean=region_df$vaccinated,
                                   lci=region_df$vaccinated_lci, uci=region_df$vaccinated_uci, stringsAsFactors=FALSE)
  region_vials <- data.frame(Level=region_df$region, year=region_df$year, scenario=region_df$scenario,
                              variable="Vaccine vials", mean=region_df$total_vials,
                              lci=region_df$total_vials_lci, uci=region_df$total_vials_uci, stringsAsFactors=FALSE)
  
  # Bind dataframes together
  plot_df <- rbind( region_deaths, region_vaccinated, region_vials)
  
  # Set level as factor
  plot_df$Level = factor(plot_df$Level, levels=c("China","Shaanxi","Tianjin","Shandong","Guangxi","Hunan","Guizhou"))
  
  # Add text to act as text for plot panels
  plot_df$text <- NA
  for(i in 1:nrow(plot_df)){
    plot_df$text[i] <- paste0(plot_df$Level[i], "-", plot_df$variable[i])
  }
  
  return(plot_df)
}

############################################################################
#5.Functions to prepare for stacked incremental cost effectiveness analysis#
############################################################################

stacked_incremental_cost_analysis_province <- function(global_df,  scenario_1, scenario_2, scenario_3=NULL,scenario_4=NULL){
  
  # Create global row
  global_row <- data.frame(Level=unique(global_df$region),
                           total_deaths_1 = global_df$total_deaths[global_df$scenario==scenario_1],
                           total_cost_1 = global_df$total_cost[global_df$scenario==scenario_1],
                           total_deaths_2 = global_df$total_deaths[global_df$scenario==scenario_2],
                           total_cost_2 = global_df$total_cost[global_df$scenario==scenario_2])
  
 
  # Join dataframes together
  bound_df <- rbind(global_row)
  
  # Calculate differences
  bound_df$deaths_diff <- bound_df$total_deaths_2-bound_df$total_deaths_1
  bound_df$cost_diff <- bound_df$total_cost_2-bound_df$total_cost_1
  
  # Calculate Incremental cost
  bound_df$ICER_deaths = bound_df$cost_diff/bound_df$deaths_diff
  
  # If a 3rd scenario is used, calculate and cbind results to df
  if(!is.null(scenario_3)){
    # Create global row
    global_row <- data.frame(Level=unique(global_df$region),
                             total_deaths_3 = global_df$total_deaths[global_df$scenario==scenario_3],
                             total_cost_3 = global_df$total_cost[global_df$scenario==scenario_3])
    
    # Join dataframes together
    extra_df <- rbind(global_row)
    bound_df <- cbind(bound_df, extra_df[,-1])
    
    # Calculate differences
    bound_df$deaths_diff_3 <- bound_df$total_deaths_3-bound_df$total_deaths_1
    bound_df$cost_diff_3 <- bound_df$total_cost_3-bound_df$total_cost_1
    
    # Calculate Incremental cost
    bound_df$ICER_deaths_3 = bound_df$cost_diff_3/bound_df$deaths_diff_3
  }
  
  if(!is.null(scenario_4)){
    
    # Create global row
    global_row <- data.frame(Level=unique(global_df$region),
                             total_deaths_4 = global_df$total_deaths[global_df$scenario==scenario_4],
                             total_cost_4 = global_df$total_cost[global_df$scenario==scenario_4])
    
    # Join dataframes together
    extra_df <- rbind(global_row)
    bound_df <- cbind(bound_df, extra_df[,-1])
    
    # Calculate differences
    bound_df$deaths_diff_4 <- bound_df$total_deaths_4-bound_df$total_deaths_1
    bound_df$cost_diff_4 <- bound_df$total_cost_4-bound_df$total_cost_1
    
    # Calculate Incremental cost
    bound_df$ICER_deaths_4 = bound_df$cost_diff_4/bound_df$deaths_diff_4
  }
  
  # Stack dataframes
  ICER_stacked_df <- data.frame(Level=bound_df$Level, deaths_diff=bound_df$deaths_diff,
                                cost_diff=bound_df$cost_diff, ICER_deaths=bound_df$ICER_deaths,
                                Scenario=paste0(scenario_1, " vs. ", scenario_2))
  
  if(!is.null(scenario_3)){
    stacked_df_p2 <- data.frame(Level=bound_df$Level, deaths_diff=bound_df$deaths_diff_3,
                                cost_diff=bound_df$cost_diff_3, ICER_deaths=bound_df$ICER_deaths_3,
                                Scenario=paste0(scenario_1, " vs. ", scenario_3))
    ICER_stacked_df <- rbind(ICER_stacked_df, stacked_df_p2)
  }
  if(!is.null(scenario_4)){
    stacked_df_p3 <- data.frame(Level=bound_df$Level, deaths_diff=bound_df$deaths_diff_4,
                                cost_diff=bound_df$cost_diff_4, ICER_deaths=bound_df$ICER_deaths_4,
                                Scenario=paste0(scenario_1, " vs. ", scenario_4))
    ICER_stacked_df <- rbind(ICER_stacked_df, stacked_df_p3)
  }
  
  return(ICER_stacked_df)
  
}

stacked_incremental_cost_analysis <- function(global_df,  scenario_1, scenario_2, scenario_3=NULL,scenario_4=NULL){
  
  # Create global row
  global_row <- data.frame(Level="all regions",
                           total_deaths_1 = global_df$total_deaths[global_df$scenario==scenario_1],
                           total_cost_1 = global_df$total_cost[global_df$scenario==scenario_1],
                           total_deaths_2 = global_df$total_deaths[global_df$scenario==scenario_2],
                           total_cost_2 = global_df$total_cost[global_df$scenario==scenario_2])
  
  
  # Join dataframes together
  bound_df <- rbind(global_row)
  
  # Calculate differences
  bound_df$deaths_diff <- bound_df$total_deaths_2-bound_df$total_deaths_1
  bound_df$cost_diff <- bound_df$total_cost_2-bound_df$total_cost_1
  
  # Calculate Incremental cost
  bound_df$ICER_deaths = bound_df$cost_diff/bound_df$deaths_diff
  
  # If a 3rd scenario is used, calculate and cbind results to df
  if(!is.null(scenario_3)){
    # Create global row
    global_row <- data.frame(Level="all regions",
                             total_deaths_3 = global_df$total_deaths[global_df$scenario==scenario_3],
                             total_cost_3 = global_df$total_cost[global_df$scenario==scenario_3])
    
    # Join dataframes together
    extra_df <- rbind(global_row)
    bound_df <- cbind(bound_df, extra_df[,-1])
    
    # Calculate differences
    bound_df$deaths_diff_3 <- bound_df$total_deaths_3-bound_df$total_deaths_1
    bound_df$cost_diff_3 <- bound_df$total_cost_3-bound_df$total_cost_1
    
    # Calculate Incremental cost
    bound_df$ICER_deaths_3 = bound_df$cost_diff_3/bound_df$deaths_diff_3
  }
  
  if(!is.null(scenario_4)){
    # Create global row
    global_row <- data.frame(Level="all regions",
                             total_deaths_4 = global_df$total_deaths[global_df$scenario==scenario_4],
                             total_cost_4 = global_df$total_cost[global_df$scenario==scenario_4])
    
    # Join dataframes together
    extra_df <- rbind(global_row)
    bound_df <- cbind(bound_df, extra_df[,-1])
    
    # Calculate differences
    bound_df$deaths_diff_4 <- bound_df$total_deaths_4-bound_df$total_deaths_1
    bound_df$cost_diff_4 <- bound_df$total_cost_4-bound_df$total_cost_1
    
    # Calculate Incremental cost
    bound_df$ICER_deaths_4 = bound_df$cost_diff_4/bound_df$deaths_diff_4
  }
  
  # Stack dataframes
  ICER_stacked_df <- data.frame(Level=bound_df$Level, deaths_diff=bound_df$deaths_diff,
                                cost_diff=bound_df$cost_diff, ICER_deaths=bound_df$ICER_deaths,
                                Scenario=paste0(scenario_1, " vs. ", scenario_2))
  
  if(!is.null(scenario_3)){
    stacked_df_p2 <- data.frame(Level=bound_df$Level, deaths_diff=bound_df$deaths_diff_3,
                                cost_diff=bound_df$cost_diff_3, ICER_deaths=bound_df$ICER_deaths_3,
                                Scenario=paste0(scenario_1, " vs. ", scenario_3))
    ICER_stacked_df <- rbind(ICER_stacked_df, stacked_df_p2)
  }
  if(!is.null(scenario_4)){
    stacked_df_p3 <- data.frame(Level=bound_df$Level, deaths_diff=bound_df$deaths_diff_4,
                                cost_diff=bound_df$cost_diff_4, ICER_deaths=bound_df$ICER_deaths_4,
                                Scenario=paste0(scenario_1, " vs. ", scenario_4))
    ICER_stacked_df <- rbind(ICER_stacked_df, stacked_df_p3)
  }
  
  return(ICER_stacked_df)
  
}

########################### End of the code ###########################