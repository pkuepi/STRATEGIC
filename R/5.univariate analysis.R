#################################################################
#                                                               #
#  Strategies to inTerrupt RAbies Transmission for the          #
#  Elimination Goal by 2030 In China (STRATEGIC) study          #
#                                                               #        
#  5.univariate analysis                                        #
#  This file is used for one-way sensitivity analyses           #
#  Modified from the WHO rabies modelling consortium study      #
#  URL: https://github.com/katiehampson1978/rabies_PEP_access   #
#  Ref.: The Lancet Infectious Diseases 2019; 19(1): 102-11.    #
#                                                               #
#################################################################


#########################
#  univariate analysis  #
#########################
univariate_analysis <- function(ndraw, horizon, GAVI_status, DogVax_TF, VaxRegimen,
                                DALYrabies, LE, RIG_status, discount, breaks, IBCM){
  
  variables<-c("suspect_bite_incidence","non_suspect_bite_incidence","p_transmission","p_prevent")
  
  results<-vector("list", length(variables)*nrow(data))
  
  for(j in 1:length(variables)){
    for(i in 1:nrow(data)){
      set.seed(5) # set the seed before running all the draws for this region - same random numbers will be drawn for each region
      index<-nrow(data)*(j-1)+i
      results[[index]]<-cbind.data.frame(decision_tree_univ_ndraw(ndraw=ndraw, vary=variables[j], region=data$region[i],
                                                                  horizon=horizon, GAVI_status=GAVI_status, DogVax_TF=DogVax_TF,
                                                                  VaxRegimen=VaxRegimen, DALYrabies=DALYrabies, LE=LE,
                                                                  RIG_status=RIG_status, discount=discount, breaks=breaks, IBCM = IBCM),
                                         region=data$region[i], vary=variables[j])
      # Print % completion to console
      print(paste(round(100/nrow(data)*i, digits=1), "%", sep=" "))
    }
    # Print variable
    print(variables[j])
  }
  
  out<-do.call(rbind.data.frame, results)
  
  # Label data.frame names properly
  colnames(out) <- c("year","total_deaths","total_vials","total_cost","total_deaths_averted",
                     "total_YLL","total_YLL_averted","vaccinated","fully_vaccinated",
                     "p_seek_rabid", "p_seek_healthy", "p_receive", "p_receive_IBCM", "p_receive_RIG", "p_complete",
                     "TargetPopulation_rabid","TargetPopulation_healthy","RIG","Gavi_support",
                     "iter","region","variable")
  return(out)
}

############################
#  1.decision tree n_draw  #
############################

decision_tree_univ_ndraw<-function(ndraw, vary, region, horizon,
                                   GAVI_status, DogVax_TF, VaxRegimen,
                                   DALYrabies, LE, RIG_status,
                                   discount, breaks, IBCM){
  draws<-vector("list",ndraw)
  for(i in 1:ndraw){
    draws[[i]]<-cbind.data.frame(decision_tree_univ_draw(vary=vary, region=region, horizon=horizon,
                                                         GAVI_status=GAVI_status, DogVax_TF=DogVax_TF, VaxRegimen=VaxRegimen,
                                                         DALYrabies=DALYrabies, LE=LE, RIG_status=RIG_status,
                                                         discount=discount, breaks=breaks, IBCM = IBCM),
                                 iter=i)
  }
  draws<-do.call("rbind",draws)
  return(draws)
}

#
region_horizon_iter_univar <- function(model_output){
  
  sums = model_output %>%
    dplyr::group_by(region, scenario, variable, iter) %>%
    dplyr::summarise(cluster = unique(cluster),
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

##########################
#  2.decision tree draw  #
##########################

decision_tree_univ_draw<-function(vary,
                                  region, horizon,
                                  GAVI_status, DogVax_TF,
                                  VaxRegimen, DALYrabies, LE, RIG_status, discount, breaks, IBCM){
  
  # Ensure all inputs are correct:
  if(!GAVI_status%in%c("none","base","low","high")){
    stop("'GAVI_status' must take one of the following arguments: 'none','base','low','high'.")
  }
  if(!is.logical(DogVax_TF)){stop("DogVax_TF must be either true or false.")}
  
  if(!RIG_status %in% c("none","high risk","all")){stop("'RIG_status' must take one of the following arguments: 'none', 'high risk', or 'all'.")}
  
  if(!VaxRegimen %in% c("Updated TRC","IPC")){stop("'Vaccine regimen' must take one of the following arguments: 'Updated TRC', or 'IPC'.")}
  
  if(!is.logical(IBCM)){stop("IBCM must be either true or false")}
  
  # Set conditions:
  
  # 1). region
  index = which(data$region==region)
  gavi_phase=2023
  
  # 2). Bite incidence - from suspect rabid dogs
  # create human pop table for the horizon of interest: 2023 - 2035
  pop2023_col <- which(colnames(data)=="pop2023")
  pop_cols <- c(pop2023_col:(pop2023_col+(horizon-1)))
  pop <- data[index, pop_cols]
  # create dog table for the horizon of interest
  dogs <- dogs[which(dogs$region == region),] # make sure data and dogs have same regions
  dog2023_col <- which(colnames(dogs)=="dogpop2023")
  dogpop_cols <- c(dog2023_col:(dog2023_col + (horizon-1)))
  dogpop <- dogs[,dogpop_cols]
  # mean bites per rabid dog
  rabid_dog_bites = 0.38
  
  if(DogVax_TF==T) {
    Y1 = grep("2023", names(elimination_traj))
    yr_index = Y1:(Y1+horizon-1)
    trajectories = elimination_traj[, yr_index]
    
    if(vary=="suspect_bite_incidence"){ # select a random number from 0-100 to select an incidence trajectory
      elim_index <- round(runif(1, 1, 100))
      ET <- unlist(trajectories[elim_index,])
    } else { ET <- unlist(apply(trajectories, 2, mean))}
    
    # translate into a region-SPECIFIC rabid dog bite incidence
    susp_bite_inc_traj = as.numeric((dogpop * rabid_dog_bites * ET)/ pop)
    susp_bite_inc <- unlist(susp_bite_inc_traj)
    
  } else {
    if(vary=="suspect_bite_incidence"){ 
      susp_bite_inc <-  (dogpop[1:horizon] * rabid_dog_bites * rnorm(horizon, mean = rabies$inc, sd = rabies$sd))/pop
    } else { susp_bite_inc <- (dogpop[1:horizon] * rabid_dog_bites * rep(rnorm(1, mean = rabies$inc, sd = rabies$sd), horizon))/pop[1:horizon]}
  }
  susp_bite_inc <- unlist(susp_bite_inc)
  
  # 3). Bite incidence - from non-suspect dogs
  mean <- data$bite_inc_non_susp[index]
  if(vary=="non_suspect_bite_incidence"){
    non_susp_bite_inc <- rep(rtriangle(n=1, a=mean-0.2*mean, b=mean+0.2*mean), horizon)
  } else {non_susp_bite_inc <- rep(mean, horizon)}
  
  # 4). Probability of developing rabies if bitten by rabid dog - Differences due to age and bite site ignored - separate analysis justifies this
  if(vary=="p_transmission"){
    prob_transmission <- mixture_model(1) #### Varying p_transmission
  }else{
    prob_transmission <- (transmission$prob_head * transmission$prob_death_given_head) +
      (transmission$prob_arm * transmission$prob_death_given_arm) +
      (transmission$prob_trunk * transmission$prob_death_given_trunk) +
      (transmission$prob_leg * transmission$prob_death_given_leg)
  }
  if(GAVI_status=="base") {
    
    # 5). Health seeking
    prob_seek_pep <- gavi_phasing(gavi_phase, p_SQ = data$p_seek[index], p_step$base, p_increment, p_seek_cap$base)
    
    # 6). Access to PEP
    prob_receive_pep <- gavi_phasing(phase = gavi_phase, p_SQ = data$p_receive[index], step = p_step$base, increment = p_increment, cap = p_receive_cap$base) 
    
    # 7). Rate of PEP completion
    prob_complete_pep <- gavi_phasing(phase = gavi_phase, p_SQ = data$p_complete[index], step = p_step$base, increment = p_increment, cap = p_complete_cap$base)
    
  }else if(GAVI_status=="low") { warning("Univariate analysis only for base GAVI support.")
  } else if(GAVI_status=="high"){ warning("Univariate analysis only for base GAVI support.")
  }else{
    prob_seek_pep <- data$p_seek[index]
    prob_receive_pep <- data$p_receive[index]
    prob_complete_pep <- data$p_complete[index]
  }
  # 8). Access to RIG
  RIG_risk = data$prop_III[index] * data$p_rig[index]
  if(RIG_status=="none"){ prob_RIG <-0
  }else{prob_RIG <- RIG_risk}
  
  # 9). Probability that PEP prevents rabies
  # Depends on completion and RIG, and may differ for IM vs ID - see PEP simulations & Tz data analysis
  if(vary=="p_prevent"){
    bin_confint_no_rigrisk <- Hmisc::binconf(params$p_prevent_given_complete_n, n=params$p_prevent_given_complete_n)
    bin_confint_rigrisk <- 1 
    
    prob_prevent_rabies_given_complete_pep <- rtriangle(n=1, a=bin_confint_no_rigrisk[2], b=bin_confint_no_rigrisk[3])
    
    prob_prevent_rabies_given_complete_pep_RIG <- ((1-RIG_risk) * prob_prevent_rabies_given_complete_pep) + (RIG_risk * bin_confint_rigrisk) # No marginal survival benefit of RIG
    
    prob_prevent_rabies_given_imperfect_pep <- rbinom(n=1, size=params$p_prevent_given_imperfect_n,
                                                      prob=params$p_prevent_given_imperfect)/params$p_prevent_given_imperfect_n
    
    
  } else {
    bin_confint_no_rigrisk <- Hmisc::binconf(params$p_prevent_given_complete_n, n=params$p_prevent_given_complete_n)
    bin_confint_rigrisk <- 1
    
    prob_prevent_rabies_given_complete_pep <- mean(rtriangle(n=1000000, a=bin_confint_no_rigrisk[2], b=bin_confint_no_rigrisk[3]))
    prob_prevent_rabies_given_complete_pep_RIG <- ((1-RIG_risk) * prob_prevent_rabies_given_complete_pep) + (RIG_risk * bin_confint_rigrisk)
    prob_prevent_rabies_given_imperfect_pep <- params$p_prevent_given_imperfect
    
  }
  
  # 10). Vaccine usage
  pep = PEP_ts(gavi_phase, VaxRegimen, region)
  
  # 11). Population/ demographic data for region
  population_columns <- data[index,grep("pop20", colnames(data))]
  offst <- which(colnames(population_columns)=="pop2023")
  population <- population_columns[1:horizon+(offst-1)]
  
  ################################
  # Calculate decision tree:     #
  ################################
  
  # Measures in the decision tree
  susp_bites = susp_bite_inc * population
  # print(susp_bite_inc)
  non_susp_bites = non_susp_bite_inc * population
  
  # IBCM
  IBCM_endemic = 0.5
  IBCM_elim = 0.1
  
  # No IBCM means patients of bites by healthy and rabid dogs are treated the same
  if(IBCM == F) {
    patients = (susp_bites * prob_seek_pep) + (non_susp_bites * prob_seek_pep)
    prob_receive_pep_IBCM = prob_receive_pep
  } else {
    # IBCM means patients of bites by healthy dogs are given a risk assessment
    # If rabies is still endemic - treat just 50% of healthy bite patients
    # If rabies eliminated - treat just 10% of healthy bite patients
    prob_receive_pep_IBCM = rep(prob_receive_pep[1], length(prob_receive_pep))
    prob_receive_pep_IBCM[which(gavi_intro(gavi_phase)=="support")] <- IBCM_endemic
    prob_receive_pep_IBCM[which(susp_bite_inc==0)] <- IBCM_elim
    patients = (susp_bites * prob_seek_pep) + (non_susp_bites * prob_seek_pep * prob_receive_pep_IBCM)
  }
  
  # Human rabies deaths & deaths averted
  if(RIG_risk == 0){
    human_rabies_deaths = susp_bites * prob_transmission *
      (prob_seek_pep * prob_receive_pep * prob_complete_pep * (1-prob_prevent_rabies_given_complete_pep) + #### Those who complete PEP
         prob_seek_pep * prob_receive_pep * (1-prob_complete_pep)  * (1-prob_prevent_rabies_given_imperfect_pep) + #### Those who do not complete PEP
         prob_seek_pep * (1-prob_receive_pep) + #### Those who seek PEP but don't receive it
         (1-prob_seek_pep)) #### Those who don't seek PEP
    
    human_rabies_deaths_averted = susp_bites * prob_transmission *
      (prob_seek_pep * prob_receive_pep * prob_complete_pep * prob_prevent_rabies_given_complete_pep +
         prob_seek_pep * prob_receive_pep * (1-prob_complete_pep) * prob_prevent_rabies_given_imperfect_pep)
    
  } else {
    prob_prevent_rabies_given_complete_pep <- rep(prob_prevent_rabies_given_complete_pep_RIG, length(gavi_intro(gavi_phase))) # create p_prevent variable adjusted for RIG support
    
    human_rabies_deaths = susp_bites * prob_transmission * #### Top half of tree - those infected
      (prob_seek_pep * prob_receive_pep * prob_complete_pep * (1-prob_prevent_rabies_given_complete_pep) + #### Those who complete PEP
         prob_seek_pep * prob_receive_pep * (1-prob_complete_pep)  * (1-prob_prevent_rabies_given_imperfect_pep) + #### Those who do not complete PEP
         prob_seek_pep * (1-prob_receive_pep) + #### Those who seek PEP but don't receive it
         (1-prob_seek_pep)) #### Those who don't seek PEP
    
    human_rabies_deaths_averted = susp_bites * prob_transmission *
      (prob_seek_pep * prob_receive_pep * prob_complete_pep * prob_prevent_rabies_given_complete_pep +
         prob_seek_pep * prob_receive_pep * (1-prob_complete_pep) * prob_prevent_rabies_given_imperfect_pep)
  }
  # Vaccine vials used per year
  vials_per_year = (patients * prob_receive_pep * prob_complete_pep * pep$vials_complete) +
    (patients * prob_receive_pep * (1-prob_complete_pep) * pep$vials_imperfect)
  
  # Costs of PEP - discounted
  future <- (1:horizon)+2
  cost_PEP_per_year = patients * (prob_receive_pep * (prob_complete_pep * pep$costs_complete +
                                                        (1-prob_complete_pep) * pep$costs_imperfect)) * exp(-discount*future)
  
  # Courses of RIG per year - and costs discounted
  p_RIG <- rep(prob_RIG, length(gavi_intro(gavi_phase))) # create p_RIG adjusted for RIG support
  courses_RIG_per_year = patients * prob_receive_pep * p_RIG
  cost_RIG = data$price_per_vial_RIG[index]
  cost_RIG_per_year = courses_RIG_per_year * cost_RIG * exp(-discount*future)
  
  if(DogVax_TF==T){
    cost_dog_per_year = c(11.33*6.8996*0.4*dogpop[1] * exp(-discount*future[1]),
                          11.33*6.8996*0.7*dogpop[2:horizon] * exp(-discount*future[2:horizon]))
    }else{
    cost_dog_per_year = 11.33*6.8996*0.4*dogpop[1:horizon] * exp(-discount*future)
  }
  if(IBCM == F){
    cost_IBCM=rep(0,horizon)
  }else{
    cost_IBCM=c(100000*6.8996*exp(-3),rep(0,horizon-1))
  }
  
  # TOTAL COSTS
  cost_per_year = cost_PEP_per_year + cost_RIG_per_year + cost_dog_per_year + cost_IBCM
  
  # DALYs for rabies (discounted YLL)
  YLL_rabies_case = YLLcalc(DALYtable=DALYrabies, LTvalues=LE, cause="rabies", discount=discount, C=1, Beta=0, alpha=0, breaks=breaks) * exp(-discount*future) # check not discounting twice!
  YLL_rabies = human_rabies_deaths * YLL_rabies_case
  YLL_averted = human_rabies_deaths_averted * YLL_rabies_case
  
  # Vaccinated and fully vaccinated persons
  vaccinated = patients * prob_receive_pep
  fully_vaccinated = vaccinated * prob_complete_pep
  
  years<-2023:2070
  
  ## p_seek:suspect and p_seek:healthy.
  p_seek_rabid <- prob_seek_pep
  p_seek_healthy <- prob_seek_pep
  
  ## p_receive.
  p_receive <- prob_receive_pep
  p_receive_RIG <- p_RIG
  
  ## p_complete.
  p_complete <- prob_complete_pep
  
  ## "Target population".
  population <- as.data.frame(t(population))
  rownames(population) <- NULL
  
  # target1: genuinely rabid exposed persons
  exposure_inc <- susp_bite_inc
  TargetPopulation_rabid <- as.numeric(unlist(exposure_inc * t(population)))
  
  # target2: persons bitten by healthy animals
  healthy_exposure_inc <- non_susp_bite_inc
  TargetPopulation_healthy <- as.numeric(unlist(healthy_exposure_inc * population))
  
  # RETURN RESULTS
  return(cbind.data.frame(year=years[1:horizon],
                          human_rabies_deaths = as.numeric(human_rabies_deaths),
                          vials_per_year = as.numeric(vials_per_year),
                          cost_per_year = as.numeric(cost_per_year),
                          human_rabies_deaths_averted = as.numeric(human_rabies_deaths_averted),
                          YLL_rabies = as.numeric(YLL_rabies),
                          YLL_averted = as.numeric(YLL_averted),
                          vaccinated = as.numeric(vaccinated),
                          fully_vaccinated = as.numeric(fully_vaccinated),
                          p_seek_rabid = as.numeric(p_seek_rabid),
                          p_seek_healthy = as.numeric(p_seek_healthy),
                          p_receive = as.numeric(p_receive),
                          p_receive_IBCM = prob_receive_pep_IBCM,
                          p_receive_RIG = p_RIG,
                          p_complete = as.numeric(p_complete),
                          TargetPopulation_rabid = as.numeric(TargetPopulation_rabid),
                          TargetPopulation_healthy = as.numeric(TargetPopulation_healthy),
                          RIG = as.numeric(courses_RIG_per_year),
                          gavi_support = gavi_intro(gavi_phase)
                          #, pep_vials = pep$vials_complete # Check on vials per patient!
  )
  )
}

########################### End of the code ###########################