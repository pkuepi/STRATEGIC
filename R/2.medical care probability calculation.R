#################################################################
#                                                               #
#  Strategies to inTerrupt RAbies Transmission for the          #
#  Elimination Goal by 2030 In China (STRATEGIC) study          #
#                                                               #        
#  2.medical care probability calculation                       #
#  This file is used to calculate the probabilities of          #
#  seeking, receiving, completing medical care                  #
#  Modified from the WHO rabies modelling consortium study      #
#  URL: https://github.com/katiehampson1978/rabies_PEP_access   #
#  Ref.: The Lancet Infectious Diseases 2019; 19(1): 102-11.    #
#                                                               #
#################################################################


##### probabilities of seeking, receiving, completing medical care 
p_seek_cap  <- data.frame(base = 0.9, low = 0.85, high = 0.95) # Maximum p_seek for rabid bite victims for each scenario
p_receive_cap <- data.frame(base = 0.99, low = 0.98, high = 1) # Maximum p_receive for bite victims for each scenario
p_complete_cap <- data.frame(base = 0.975, low = 0.95, high = 1) # Maximum p_complete for bite victims for each scenario
p_step <- data.frame(base = 0.01, low = 0.01, high = 0.01) # 1% increase on introduction of Gavi support
p_increment <- 0.01 # 1% per annum as per other Gavi vaccines
gavi_phasing = function(phase, p_SQ, step, increment, cap){
  p = rep(p_SQ, length(2023:2035)) # create vector of status quo probabilities
  if(p_SQ > cap){p} # make sure that the new strategy is NEVER below the status quo!
  else
    if(p_SQ+step > cap){p[which(2023:2035 >= phase)] <- cap } # make sure that the new strategy is NEVER below the status quo!
  else{p[which(2023:2035 == phase)] <- p_SQ + step # step change with year of the new strategy
  p_traj = c(seq(from = p_SQ + step, to = cap, by = increment), rep(cap, 15)) # subsequent trajectories of probabilities
  yrs = which(2023:2035 >= phase) # fit the trajectory into the simulation time horizon
  p[yrs] = p_traj[1:length(yrs)]
  }
  p
}

########################### End of the code ###########################