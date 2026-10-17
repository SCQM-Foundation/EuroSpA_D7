
# Title: Additional targeted simulations ----------------------------------

# Notes:
# 1) For a hazard ratio of 1, the four scenarios (exponential, Weibull, 500, 1000) with a standard dosing interval of 56 days for drug A and 7 days for drug B are repeated with within-drug variability in the dosing interval.
# 2) The distribution of dosing intervals was roughly guided by the frequencies of dosing intervals observed for infliximab and etanercept in the registry data (with diagnoses combined): 
# Infliximab: 70% 56 days (standard), 20% 42 days, and 10% 14 days
# Etanercept: 90% 7 days (standard) and 10% 14 days


# Set-up ------------------------------------------------------------------

# rm(list = ls()) #For interactive use
# .rs.restartR() #For interactive use. Or type "R --vanilla" in terminal for a completely fresh R process

library(dplyr)
library(survival)
library(tram)

load(file="Simulation/simulation_results/Sim_param.Rdata")


# Simulation parameters ---------------------------------------------------

nsim <- 2000

dist <- c("exp", "wei")
exp_rate_A <- exp_rate_B <- -log(0.8)
weibull_shape_A <- weibull_shape_B <- log(-log(2)/log(0.8))/log(8) 
weibull_scale_A <- weibull_scale_B <- 1/(-log(0.8))^(1/weibull_shape_A) 

sample_size <- c(100, 500, 1000)

drugA_intervals <- c(56, 42, 14)
drugA_int_prob <- c(0.7, 0.2, 0.1)
drugB_intervals <- c(7, 14)
drugB_int_prob <- c(0.9, 0.1)

last_received <- c("A", "B", "C", "D", "E")
decision <- c("F", "G", "H")
first_missed <- c("I", "J")

# > Simulation DGMs and repetitions: --------------------------------------
sim_scen <- data.frame(DGM = rep(1:6, each=nsim), dist = rep(dist, each=3*nsim), nobs = rep(rep(sample_size, times=2), each=nsim), standard_dosing_intervalA = rep(56, 6*nsim), standard_dosing_intervalB = rep(7, 6*nsim), HR = rep(1, 6*nsim), iteration = rep(1:nsim, times=6))


# Simulation --------------------------------------------------------------

sim_dat <- vector(mode="list", length=nrow(sim_scen))

tdata_state <- vector(mode="list", length=nrow(sim_scen))
tdata_state[[1]] <- tdg_state

tRNG_state <- vector(mode="list", length=nrow(sim_scen))

start_time <- Sys.time()

for (i in 1:nrow(sim_scen)) {
  
  .Random.seed <- tdata_state[[i]]
  
  sim_dat[[i]] <- data.frame(model=rep(NA_character_, 5), est_logHR=rep(NA_real_, 5), var=rep(NA_real_, 5), pval=rep(NA_real_, 5))
  
  # > Data generation: ------------------------------------------------
  
  df <- data.frame(Registry=as.factor(rep(rep(LETTERS[1:10], each=sim_scen$nobs[i]/10), times=2)), FU=FU, Drug=as.factor(rep(c("A", "B"), each=sim_scen$nobs[i])))
  
  if(sim_scen$dist[i] == "exp") {
    df[df$Drug == "A", "Disc_time"] <- ceiling(rexp(n=sim_scen$nobs[i], rate=exp_rate_A)*365)
    df[df$Drug == "B", "Disc_time"] <- ceiling(rexp(n=sim_scen$nobs[i], rate=exp_rate_B)*365)
  } else {
    df[df$Drug == "A", "Disc_time"] <- ceiling(rweibull(n=sim_scen$nobs[i], shape=weibull_shape_A, scale=weibull_scale_A)*365)
    df[df$Drug == "B", "Disc_time"] <- ceiling(rweibull(n=sim_scen$nobs[i], shape=weibull_shape_B, scale=weibull_scale_B)*365)
  }
  
  df[df$Drug == "A", "Dosing"] <- sample(x=drugA_intervals, size=sim_scen$nobs[i], replace=TRUE, prob=drugA_int_prob)
  df[df$Drug == "B", "Dosing"] <- sample(x=drugB_intervals, size=sim_scen$nobs[i], replace=TRUE, prob=drugB_int_prob)
  
  df <- df %>% 
    mutate(Last_received = Dosing * Disc_time %/% Dosing) %>% 
    mutate(First_missed = Dosing * (Disc_time %/% Dosing + 1)) %>% 
    
    mutate(decision_event = as.numeric(Disc_time <= FU)) %>% 
    mutate(time_left = ifelse(test=decision_event == 0, yes=FU, no=Last_received)) %>% 
    mutate(time_right = ifelse(test=decision_event == 0, yes=Inf, no=pmin(First_missed, FU))) %>% 
    
    as.data.frame()
  
  if (i < nrow(sim_scen)) tdata_state[[i+1]] <- .Random.seed
  
  # > Analyses: -------------------------------------------------------
  # >> SM1: -----------------------------------------------------------
  # Everyone reports the decision to discontinue if before end of follow-up.
  
  fit_SM1 <- survival::coxph(Surv(time=pmin(Disc_time, FU), event=decision_event) ~ Drug, data=df)
  
  sim_dat[[i]]$model[1] <- "SM1"
  sim_dat[[i]]$est_logHR[1] <- fit_SM1$coefficients
  sim_dat[[i]]$var[1] <- fit_SM1$var
  sim_dat[[i]]$pval[1] <- summary(fit_SM1)$coefficients[5]
  
  if (is.null(fit_SM1$coefficients) | is.null(fit_SM1$var) | is.null(summary(fit_SM1)$coefficients[5])) cat(i, ": no SM1 fit") 
  
  # >> SM2: -----------------------------------------------------------
  # Either the date of the last dose received or the date of the first dose missed is reported if the decision to discontinue happened before end of follow-up. From that, information on an interval for the decision is derived.
  
  tRNG_state[[i]][["before_tram"]] <- list(fit_state[[1]][[i]], fit_state[[2]][[i]])
  tRNG_state[[i]][["after_tram"]] <- fit <- vector(mode="list", length=2)
  
  for (j in 1:2) {
    .Random.seed <- fit_state[[j]][[i]]
    
    fit[[j]] <- tram::Coxph(Surv(time=time_left, time2=time_right, type="interval2") ~ Drug , data=df, log_first=T)
    
    tRNG_state[[i]][["after_tram"]][[j]] <- .Random.seed
  }
  
  fit_SM2 <- if (as.numeric(logLik(fit[[1]])) >= as.numeric(logLik(fit[[2]]))) fit[[1]] else fit[[2]]
  
  sim_dat[[i]]$model[2] <- "SM2"
  sim_dat[[i]]$est_logHR[2] <- fit_SM2$coef["DrugB"]
  sim_dat[[i]]$var[2] <- summary(fit_SM2)$test$vcov
  sim_dat[[i]]$pval[2] <- summary(fit_SM2)$test$test$pvalues
  
  if (is.null(fit_SM2$coef["DrugB"]) | is.null(summary(fit_SM2)$test$vcov) | is.null(summary(fit_SM2)$test$test$pvalues)) cat(i, ": no SM2 fit") 
  
  # >> SM3: -----------------------------------------------------------
  # Everyone reports the date of the last dose received if the decision to discontinue happened before end of follow-up. Time to the last dose received is erroneously right-censored at end of follow-up if the decision is taken less than one dosing interval after end of follow-up.
  
  fit_SM3 <- survival::coxph(Surv(time=ifelse(test=decision_event == 0, yes=FU, no=Last_received), event=decision_event) ~ Drug, data=df)
  
  sim_dat[[i]]$model[3] <- "SM3"
  sim_dat[[i]]$est_logHR[3] <- fit_SM3$coefficients
  sim_dat[[i]]$var[3] <- fit_SM3$var
  sim_dat[[i]]$pval[3] <- summary(fit_SM3)$coefficients[5]
  
  if (is.null(fit_SM3$coefficients) | is.null(fit_SM3$var) | is.null(summary(fit_SM3)$coefficients[5])) cat(i, ": no SM3 fit") 
  
  # >> SM4: -----------------------------------------------------------
  # Everyone reports the date of the first dose missed if the decision to discontinue happened before end of follow-up. In case the decision is taken less than one dosing interval before the end of follow-up, the time to the first dose missed is correctly right-censored at end of follow-up.
  
  fit_SM4 <- survival::coxph(Surv(time=pmin(First_missed, FU), event=First_missed <= FU) ~ Drug, data=df)
  
  sim_dat[[i]]$model[4] <- "SM4"
  sim_dat[[i]]$est_logHR[4] <- fit_SM4$coefficients
  sim_dat[[i]]$var[4] <- fit_SM4$var
  sim_dat[[i]]$pval[4] <- summary(fit_SM4)$coefficients[5]
  
  if (is.null(fit_SM4$coefficients) | is.null(fit_SM4$var) | is.null(summary(fit_SM4)$coefficients[5])) cat(i, ": no SM4 fit")
  
  # >> SM5: -----------------------------------------------------------
  # Mixture of reported events.
  
  df_SM5 <- df %>% 
    mutate(type = case_when(
      Registry %in% decision ~ "decision",
      Registry %in% last_received ~ "last_received",
      Registry %in% first_missed ~ "first_missed"
    )) %>% 
    mutate(event_SM5 = case_when(
      type %in% c("decision", "last_received") ~ decision_event,
      type == "first_missed" ~ as.numeric(First_missed <= FU)
    )) %>% 
    mutate(time_SM5 = case_when(
      type == "decision" ~ pmin(Disc_time, FU),
      type == "last_received" ~ ifelse(test=decision_event == 0, yes=FU, no=Last_received),
      type == "first_missed" ~ pmin(First_missed, FU)
    )) %>% 
    as.data.frame()
  
  fit_SM5 <- survival::coxph(Surv(time=time_SM5, event=event_SM5) ~ Drug + strata(Registry), data=df_SM5)
  
  sim_dat[[i]]$model[5] <- "SM5"
  sim_dat[[i]]$est_logHR[5] <- fit_SM5$coefficients
  sim_dat[[i]]$var[5] <- fit_SM5$var
  sim_dat[[i]]$pval[5] <- summary(fit_SM5)$coefficients[5]
  
  if (is.null(fit_SM5$coefficients) | is.null(fit_SM5$var) | is.null(summary(fit_SM5)$coefficients[5])) cat(i, ": no SM5 fit") 
  
}

end_time <- Sys.time()
run_time <- end_time - start_time

# > Estimates data long format: -------------------------------------------
est_df <- bind_rows(sim_dat, .id="repetition") %>%
  mutate(se = sqrt(var)) %>% 
  mutate(DGM = rep(sim_scen$DGM, each=5)) %>%
  mutate(dist = rep(sim_scen$dist, each=5)) %>% 
  mutate(nobs = rep(sim_scen$nobs, each=5)) %>% 
  mutate(standard_dosing_intervalA = rep(sim_scen$standard_dosing_intervalA, each=5)) %>% 
  mutate(standard_dosing_intervalB = rep(sim_scen$standard_dosing_intervalB, each=5)) %>% 
  mutate(HR = rep(sim_scen$HR, each=5)) %>% 
  mutate(logHR = log(HR)) %>% 
  as.data.frame()


# Saving output -----------------------------------------------------------

Info_session <- sessionInfo()
Info_sys <- Sys.info()

save(Info_session, Info_sys, sim_dat, est_df, sim_scen, nsim, run_time, tdata_state, tRNG_state, file="/home/sharedfolders/study_D7/simulation_results/sim_res_add_targeted_scen.Rdata")


