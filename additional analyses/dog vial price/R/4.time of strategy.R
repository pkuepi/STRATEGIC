#################################################################
#                                                               #
#  Strategies to inTerrupt RAbies Transmission for the          #
#  Elimination Goal by 2030 In China (STRATEGIC) study          #
#                                                               #        
#  4.time of strategy                                           #
#  This file is used to simulate the time of strategies         #
#  Modified from the WHO rabies modelling consortium study      #
#  URL: https://github.com/katiehampson1978/rabies_PEP_access   #
#  Ref.: The Lancet Infectious Diseases 2019; 19(1): 102-11.    #
#                                                               #
#################################################################


##### time when the strategy introduced
gavi_intro = function(phase){
  support = rep("none", length(2023:2035))
  if(phase < 2036){support[which(2023:2035 >= phase)] <- "support"}
  support
}


PEP_ts = function(phase, VaxRegimen, region){
  support = gavi_intro(phase) # Time series of the new strategy
  regimen = rep(VaxRegimen, length(support)) # Start with SQ regimen
  pep_outputs = data.frame(vials_complete = rep(NA, length(support)), # Create a data frame for pep outputs
                           vials_imperfect = rep(NA, length(support)),
                           costs_complete = rep(NA, length(support)),
                           costs_imperfect = rep(NA, length(support)))
  for(i in 1:length(support)){ # Run through years of the new strategy to create PEP outputs
    pep_outputs[i,] = PEPoutputs(regimen[i], region, as.character(support[i]))
  }
  pep_outputs
}

PEPoutputs = function(regimen, region, GAVI_status){
  # Index for identifying PEP scenario
  complete <- which(vacc$throughput=="Med" & vacc$Setting=="all" & vacc$regimen== regimen & vacc$completeness=="Complete")
  imperfect <- which(vacc$throughput=="Med" & vacc$Setting=="all" & vacc$regimen== regimen & vacc$completeness=="Incomplete")
  
  n_vials_complete <- vacc$vial[complete]
  n_vials_imperfect <- vacc$vial[imperfect]
  
  cost_complete <- n_vials_complete*data$price_per_vial[which(data$region==region)] +
    data$cost_first_visit[which(data$region==region)] +
    data$cost_followup_visit[which(data$region==region)]*(vacc$nvisit[complete]-1)
  
  cost_imperfect <- n_vials_imperfect*data$price_per_vial[which(data$region==region)] +
    data$cost_first_visit[which(data$region==region)] +
    data$cost_followup_visit[which(data$region==region)]*(vacc$nvisit[imperfect]-1)
  
  output <- data.frame(vials_complete = n_vials_complete,
                       vials_imperfect = n_vials_imperfect,
                       costs_complete = cost_complete,
                       costs_imperfect= cost_imperfect)
  return(output)
}

########################### End of the code ###########################