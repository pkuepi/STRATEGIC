#################################################################
#                                                               #
#  Strategies to inTerrupt RAbies Transmission for the          #
#  Elimination Goal by 2030 In China (STRATEGIC) study          #
#                                                               #        
#  7.results summary                                            #
#  This file is used for output of the analyses                 #
#  Modified from the WHO rabies modelling consortium study      #
#  URL: https://github.com/katiehampson1978/rabies_PEP_access   #
#  Ref.: The Lancet Infectious Diseases 2019; 19(1): 102-11.    #
#                                                               #
#################################################################


#####################
# 1.horizon summary #
#####################

region_horizon_iter <- function(model_output){
  sums = model_output %>%
    dplyr::group_by(region, scenario, iter) %>%
    dplyr::summarise(
                     total_deaths = sum(total_deaths),
                     total_vials = sum(total_vials),
                     total_cost = sum(total_cost),
                     total_deaths_averted = sum(total_deaths_averted),
                     total_YLL = sum(total_YLL),
                     total_YLL_averted = sum(total_YLL_averted),
                     vaccinated = sum(vaccinated),
                     fully_vaccinated = sum(fully_vaccinated),
                     RIG = sum(RIG))
  sums$cost_per_death_averted <-  sums$total_cost/sums$total_deaths_averted
  sums$cost_per_YLL_averted <-  sums$total_cost/sums$total_YLL_averted
  sums$deaths_averted_per_100k_vaccinated <-  sums$total_deaths_averted/sums$vaccinated/100000
  return(sums)
}

#########################
# 2.multivar_summary(1) #
#########################

multivar_region_summary = function(model_output, year){
  if(year==TRUE){arg_setting = list(quo(year), quo(region), quo(scenario))}
  if(year==FALSE){arg_setting = list(quo(region), quo(scenario))}
  
  sums <- model_output %>%
    dplyr::group_by(!!! arg_setting) %>% # Does not include model_output!
    dplyr::summarise(
      total_deaths_lci = quantile(total_deaths, 0.025),
      total_deaths_uci = quantile(total_deaths, 0.975),
      total_deaths = mean(total_deaths),
      
      total_vials_lci = quantile(total_vials, 0.025),
      total_vials_uci = quantile(total_vials, 0.975),
      total_vials = mean(total_vials),
      
      total_cost_lci = quantile(total_cost, 0.025),
      total_cost_uci = quantile(total_cost, 0.975),
      total_cost = mean(total_cost),
      
      total_deaths_averted_lci = quantile(total_deaths_averted, 0.025),
      total_deaths_averted_uci = quantile(total_deaths_averted, 0.975),
      total_deaths_averted = mean(total_deaths_averted),
      
      total_YLL_lci = quantile(total_YLL, 0.025),
      total_YLL_uci = quantile(total_YLL, 0.975),
      total_YLL = mean(total_YLL),
      
      total_YLL_averted_lci = quantile(total_YLL_averted, 0.025),
      total_YLL_averted_uci = quantile(total_YLL_averted, 0.975),
      total_YLL_averted = mean(total_YLL_averted),
      
      vaccinated_lci = quantile(vaccinated, 0.025),
      vaccinated_uci = quantile(vaccinated, 0.975),
      vaccinated = mean(vaccinated),
      
      fully_vaccinated_lci = quantile(fully_vaccinated, 0.025),
      fully_vaccinated_uci = quantile(fully_vaccinated, 0.975),
      fully_vaccinated = mean(fully_vaccinated),
      
      RIG_lci = quantile(RIG, 0.025),
      RIG_uci = quantile(RIG, 0.975),
      RIG = mean(RIG),
      
    )
  
  sums$cost_per_death_averted_lci <-  sums$total_cost_lci/sums$total_deaths_averted_uci
  sums$cost_per_death_averted_uci <-  sums$total_cost_uci/sums$total_deaths_averted_lci
  sums$cost_per_death_averted <-  sums$total_cost/sums$total_deaths_averted
  
  sums$cost_per_YLL_averted_lci <-  sums$total_cost_lci/sums$total_YLL_averted_uci
  sums$cost_per_YLL_averted_uci <-  sums$total_cost_uci/sums$total_YLL_averted_lci
  sums$cost_per_YLL_averted <-  sums$total_cost/sums$total_YLL_averted
  
  sums$deaths_averted_per_100k_vaccinated_lci <-  sums$total_deaths_averted_lci/sums$vaccinated_uci/100000
  sums$deaths_averted_per_100k_vaccinated_uci <-  sums$total_deaths_averted_uci/sums$vaccinated_lci/100000
  sums$deaths_averted_per_100k_vaccinated <-  sums$total_deaths_averted/sums$vaccinated/100000
  
  return(sums)
}

#########################
# 2.multivar_summary(2) #
#########################

multivar_summary = function(model_output, year, setting){
  # summary over entire time horizon or by year - globally
  if(year==TRUE & setting == "global"){arg_setting = list(quo(year), quo(scenario))}
  if(year==FALSE & setting == "global"){arg_setting = list(quo(scenario))}
  
  sums <- model_output %>%
    dplyr::group_by(!!! arg_setting) %>%
    dplyr::summarise(
      total_deaths_lci = sum(total_deaths_lci),
      total_deaths_uci = sum(total_deaths_uci),
      total_deaths = sum(total_deaths),
      
      total_vials_lci = sum(total_vials_lci),
      total_vials_uci = sum(total_vials_uci),
      total_vials = sum(total_vials),
      
      total_cost_lci = sum(total_cost_lci),
      total_cost_uci = sum(total_cost_uci),
      total_cost = sum(total_cost),
      
      total_deaths_averted_lci = sum(total_deaths_averted_lci),
      total_deaths_averted_uci = sum(total_deaths_averted_uci),
      total_deaths_averted = sum(total_deaths_averted),
      
      total_YLL_lci = sum(total_YLL_lci),
      total_YLL_uci = sum(total_YLL_uci),
      total_YLL = sum(total_YLL),
      
      total_YLL_averted_lci = sum(total_YLL_averted_lci),
      total_YLL_averted_uci = sum(total_YLL_averted_uci),
      total_YLL_averted = sum(total_YLL_averted),
      
      vaccinated_lci = sum(vaccinated_lci),
      vaccinated_uci = sum(vaccinated_uci),
      vaccinated = sum(vaccinated),
      
      fully_vaccinated_lci = sum(fully_vaccinated_lci),
      fully_vaccinated_uci = sum(fully_vaccinated_uci),
      fully_vaccinated = sum(fully_vaccinated),
      
      RIG_lci = sum(RIG_lci),
      RIG_uci = sum(RIG_uci),
      RIG = sum(RIG),
      cost_per_death_averted_lci = total_cost_lci/total_deaths_averted_uci,
      cost_per_death_averted_uci = total_cost_uci/total_deaths_averted_lci,
      cost_per_death_averted = total_cost/total_deaths_averted,
      
      cost_per_YLL_averted_lci = total_cost_lci/total_YLL_averted_uci,
      cost_per_YLL_averted_uci = total_cost_uci/total_YLL_averted_lci,
      cost_per_YLL_averted = total_cost/total_YLL_averted,
      
      deaths_averted_per_100k_vaccinated_lci = total_deaths_averted_lci/vaccinated_uci/100000,
      deaths_averted_per_100k_vaccinated_uci = total_deaths_averted_uci/vaccinated_lci/100000,
      deaths_averted_per_100k_vaccinated = total_deaths_averted/vaccinated/100000)
  return(sums)
}

####################
# 3.univar_summary #
####################

univar_region_summary = function(model_output, year){
  
  if(year==TRUE){arg_setting = list(quo(year), quo(region), quo(variable), quo(scenario))}
  if(year==FALSE){arg_setting = list(quo(region), quo(variable), quo(scenario))}
  
  sums <- model_output %>%
    dplyr::group_by(!!! arg_setting) %>% # Does not include model_output!
    dplyr::summarise(
      total_deaths_lci = quantile(total_deaths, 0.025),
      total_deaths_uci = quantile(total_deaths, 0.975),
      total_deaths = mean(total_deaths),
      
      total_vials_lci = quantile(total_vials, 0.025),
      total_vials_uci = quantile(total_vials, 0.975),
      total_vials = mean(total_vials),
      
      total_cost_lci = quantile(total_cost, 0.025),
      total_cost_uci = quantile(total_cost, 0.975),
      total_cost = mean(total_cost),
      
      total_deaths_averted_lci = quantile(total_deaths_averted, 0.025),
      total_deaths_averted_uci = quantile(total_deaths_averted, 0.975),
      total_deaths_averted = mean(total_deaths_averted),
      
      total_YLL_lci = quantile(total_YLL, 0.025),
      total_YLL_uci = quantile(total_YLL, 0.975),
      total_YLL = mean(total_YLL),
      
      total_YLL_averted_lci = quantile(total_YLL_averted, 0.025),
      total_YLL_averted_uci = quantile(total_YLL_averted, 0.975),
      total_YLL_averted = mean(total_YLL_averted),
      
      vaccinated_lci = quantile(vaccinated, 0.025),
      vaccinated_uci = quantile(vaccinated, 0.975),
      vaccinated = mean(vaccinated),
      
      fully_vaccinated_lci = quantile(fully_vaccinated, 0.025),
      fully_vaccinated_uci = quantile(fully_vaccinated, 0.975),
      fully_vaccinated = mean(fully_vaccinated),
      
      RIG_lci = quantile(RIG, 0.025),
      RIG_uci = quantile(RIG, 0.975),
      RIG = mean(RIG),
      
    )
  
  sums$cost_per_death_averted_lci <-  sums$total_cost_lci/sums$total_deaths_averted_uci
  sums$cost_per_death_averted_uci <-  sums$total_cost_uci/sums$total_deaths_averted_lci
  sums$cost_per_death_averted <-  sums$total_cost/sums$total_deaths_averted
  
  sums$cost_per_YLL_averted_lci <-  sums$total_cost_lci/sums$total_YLL_averted_uci
  sums$cost_per_YLL_averted_uci <-  sums$total_cost_uci/sums$total_YLL_averted_lci
  sums$cost_per_YLL_averted <-  sums$total_cost/sums$total_YLL_averted
  
  sums$deaths_averted_per_100k_vaccinated_lci <-  sums$total_deaths_averted_lci/sums$vaccinated_uci/100000
  sums$deaths_averted_per_100k_vaccinated_uci <-  sums$total_deaths_averted_uci/sums$vaccinated_lci/100000
  sums$deaths_averted_per_100k_vaccinated <-  sums$total_deaths_averted/sums$vaccinated/100000
  
  return(sums)
}

univar_summary = function(model_output, year, setting){
  
  # summary over entire time horizon or by year - globally
  if(year==TRUE & setting == "global"){arg_setting = list(quo(year), quo(variable), quo(scenario))}
  if(year==FALSE & setting == "global"){arg_setting = list(quo(variable), quo(scenario))}
  
  sums <- model_output %>%
    dplyr::group_by(!!! arg_setting) %>%
    dplyr::summarise(
      total_deaths_lci = sum(total_deaths_lci),
      total_deaths_uci = sum(total_deaths_uci),
      total_deaths = sum(total_deaths),
      
      total_vials_lci = sum(total_vials_lci),
      total_vials_uci = sum(total_vials_uci),
      total_vials = sum(total_vials),
      
      total_cost_lci = sum(total_cost_lci),
      total_cost_uci = sum(total_cost_uci),
      total_cost = sum(total_cost),
      
      total_deaths_averted_lci = sum(total_deaths_averted_lci),
      total_deaths_averted_uci = sum(total_deaths_averted_uci),
      total_deaths_averted = sum(total_deaths_averted),
      
      total_YLL_lci = sum(total_YLL_lci),
      total_YLL_uci = sum(total_YLL_uci),
      total_YLL = sum(total_YLL),
      
      total_YLL_averted_lci = sum(total_YLL_averted_lci),
      total_YLL_averted_uci = sum(total_YLL_averted_uci),
      total_YLL_averted = sum(total_YLL_averted),
      
      vaccinated_lci = sum(vaccinated_lci),
      vaccinated_uci = sum(vaccinated_uci),
      vaccinated = sum(vaccinated),
      
      fully_vaccinated_lci = sum(fully_vaccinated_lci),
      fully_vaccinated_uci = sum(fully_vaccinated_uci),
      fully_vaccinated = sum(fully_vaccinated),
      
      RIG_lci = sum(RIG_lci),
      RIG_uci = sum(RIG_uci),
      RIG = sum(RIG),
      
      cost_per_death_averted_lci = total_cost_lci/total_deaths_averted_uci,
      cost_per_death_averted_uci = total_cost_uci/total_deaths_averted_lci,
      cost_per_death_averted = total_cost/total_deaths_averted,
      
      cost_per_YLL_averted_lci = total_cost_lci/total_YLL_averted_uci,
      cost_per_YLL_averted_uci = total_cost_uci/total_YLL_averted_lci,
      cost_per_YLL_averted = total_cost/total_YLL_averted,
      
      deaths_averted_per_100k_vaccinated_lci = total_deaths_averted_lci/vaccinated_uci/100000,
      deaths_averted_per_100k_vaccinated_uci = total_deaths_averted_uci/vaccinated_lci/100000,
      deaths_averted_per_100k_vaccinated = total_deaths_averted/vaccinated/100000)
  
  return(sums)
}

########################### End of the code ###########################