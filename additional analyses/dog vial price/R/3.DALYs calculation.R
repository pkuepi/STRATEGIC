#######################################################################
#                                                                     #
#  Strategies to inTerrupt RAbies Transmission for the                #
#  Elimination Goal by 2030 In China (STRATEGIC) study                #
#                                                                     #        
#  3.DALYs calculation                                                #
#  This file is used to calculate the disability-adjusted life-years  #
#  Modified from the WHO rabies modelling consortium study            #
#  URL: https://github.com/katiehampson1978/rabies_PEP_access         #
#  Ref.: The Lancet Infectious Diseases 2019; 19(1): 102-11.          #
#                                                                     #
#######################################################################


YLLcalc <- function(DALYtable, LTvalues, cause, discount = 0.03, C = 0.1658,
                    Beta = 0.04, alpha = 1, breaks = "5y"){
  if(breaks == "annual"){interval = rep(1, length(LTvalues))}
  else {interval=c(1, 4, rep(5,17))}
  LE <- rep(LTvalues, interval)
  curr_age <- seq(from = 0, by = 1, length.out = length(LE))
  death_age <- round(curr_age + LE)
  cost <- sapply(1:(length(curr_age)), function(i){
    ages <- curr_age[i]:death_age[i]
    future <- (1:length(ages))-1
    disc_year <- exp(-discount*(future))
    weighting = ifelse(alpha==0, 1, C*ages^alpha*exp(-Beta*ages))
    return(sum(disc_year*weighting))
  })
  death_pc_prop <- DALYtable$death_pc/sum(DALYtable$death_pc)
  burden_year <- rep(death_pc_prop/interval, interval)
  return(sum(cost*burden_year))
}

########################### End of the code ###########################