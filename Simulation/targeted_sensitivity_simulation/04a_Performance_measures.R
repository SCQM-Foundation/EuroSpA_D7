
# Title: Estimating performance measures ----------------------------------
# From additional targeted simulations


# Set-up ------------------------------------------------------------------

# rm(list = ls()) #For interactive use
# .rs.restartR() #For interactive use. Or type "R --vanilla" in terminal for a completely fresh R process

library(dplyr)
library(rsimsum) 


# Loading estimates data --------------------------------------------------

load(file="/home/sharedfolders/study_D7/simulation_results/sim_res_add_targeted_scen.Rdata")


# Calculating performance measures ----------------------------------------

pm_dat <- simsum(data=est_df, estvarname="est_logHR", se="se", true="logHR", methodvar="model", ref="SM1", by="DGM", x=TRUE) 

pm_dat$summ <- pm_dat$summ %>% 
  mutate(CI_LL = est - 1.96*mcse) %>%
  mutate(CI_UL = est + 1.96*mcse) %>% 
  mutate(DGM = as.numeric(DGM)) %>% 
  
  left_join(y=sim_scen %>% select(-iteration) %>% distinct) %>% 
  
  as.data.frame()


# Saving performance measures data ----------------------------------------

save(pm_dat, file="Simulation/simulation_results/pm_data_add_targeted_scen.Rdata")

