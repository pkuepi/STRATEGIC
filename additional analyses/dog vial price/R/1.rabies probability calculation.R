#################################################################
#                                                               #
#  Strategies to inTerrupt RAbies Transmission for the          #
#  Elimination Goal by 2030 In China (STRATEGIC) study          #
#                                                               #        
#  1.rabies probability calculation                             #
#  This file is used to calculate the probability of infection  #
#  Modified from the WHO rabies modelling consortium study      #
#  URL: https://github.com/katiehampson1978/rabies_PEP_access   #
#  Ref.: The Lancet Infectious Diseases 2019; 19(1): 102-11.    #
#                                                               #
#################################################################


transmission <- read.csv("data/biological_param_p_infect.csv")
mixture_model<-function(n){ 
  N = transmission$n_rabid.bites # N BITES
  
  # BITES by site
  n_head = transmission$n_head.bites
  n_arm = transmission$n_arm.bites
  n_trunk = transmission$n_trunk.bites
  n_leg = transmission$n_leg.bites
  
  # DEATHS by site
  death_head = transmission$n_deaths.head
  death_arm = transmission$n_deaths.arm
  death_trunk = transmission$n_deaths.trunk
  death_leg = transmission$n_deaths.leg
  
  # No PEP by site
  n_no.pep.head = transmission$n_no.pep.head
  n_no.pep.arm = transmission$n_no.pep.arm
  n_no.pep.trunk = transmission$n_no.pep.trunk
  n_no.pep.leg = transmission$n_no.pep.leg
  
  # random draw of probability of bite site
  prob_head <- rbinom(n=n, size=N, prob=n_head/N)/N
  prob_arm <- rbinom(n=n, size=N, prob=n_arm/N)/N
  prob_trunk <- rbinom(n=n, size=N, prob=n_trunk/N)/N
  prob_leg <- rbinom(n=n, size=N, prob=n_leg/N)/N
  
  # random draw of probability of death given bite site
  prob_death_given_head <- rbinom(n=n, size=n_no.pep.head, prob=death_head/n_no.pep.head)/n_no.pep.head
  prob_death_given_arm <- rbinom(n=n, size=n_no.pep.arm, prob=death_arm/n_no.pep.arm)/n_no.pep.arm
  prob_death_given_trunk <- rbinom(n=n, size=n_no.pep.trunk, prob=death_trunk/n_no.pep.trunk)/n_no.pep.trunk
  prob_death_given_leg <-  rbinom(n=n, size=n_no.pep.leg, prob=death_leg/n_no.pep.leg)/n_no.pep.leg
  
  prob_death_given_rabid_bite <- prob_head*prob_death_given_head +
    prob_arm*prob_death_given_arm+
    prob_trunk*prob_death_given_trunk+
    prob_leg*prob_death_given_leg
  
  return(prob_death_given_rabid_bite)
}

########################### End of the code ###########################