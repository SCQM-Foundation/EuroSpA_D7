# ******************************************************************************
# Description: Run modules for whole process
# ****************************************************************************** 
# Note: 02_Simulation.R must be run for each combination of distribution and sample size.

# Clean environment (for interactive use)
rm(list = ls())

# Clean unused memory (for interactive use, run this always after a big run)
# .rs.restartR()

# Source functions, get packages etc
source("Simulation/settings.R")

# Load simulation parameters
load(file="Simulation/simulation_results/Sim_param.Rdata")

# For this run
sample_size <- 100 # Set either to 100, 500, or 1000
dist <- "weibull" # Set to either "weibull" or "exponential"

# Store object - Need to decide later how to deal with big data
M1_Surv_Cox <- data.frame(Est = NA, Var = NA, pval=NA, time=NA, iter=NA)
M2_Surv_Tram <- data.frame(Est = NA, Var = NA, pval=NA, time=NA, iter=NA)
M3_Surv_Cox <- data.frame(Est = NA, Var = NA, pval=NA, time=NA, iter=NA)
M4_Surv_Cox <- data.frame(Est = NA, Var = NA, pval=NA, time=NA, iter=NA)
M5_Surv_Cox <- data.frame(Est = NA, Var = NA, pval=NA, time=NA, iter=NA)

i <-1

data_state <- vector(mode="list", length=nrow(Sim_table))
data_state[[1]] <- dg_state[[paste(dist, sample_size, sep="_")]]

RNG_state <- vector(mode="list", length=nrow(Sim_table))

start_all <- Sys.time()

for (i in 1:nrow(Sim_table)){
  .Random.seed <- data_state[[i]]
  
  RNG_state[[i]] <- vector(mode="list", length=2)
  names(RNG_state[[i]]) <- c("before_tram", "after_tram")

  print(paste0("Starting iteration: ", i))
  start_data <- Sys.time()
  
  if(dist == "exponential"){
    df <- data.table(Disc_drugA = ceiling(rexp(n = sample_size,
                                               rate = Sim_table[i, HA])*365),
                     Disc_drugB = ceiling(rexp(n = sample_size,
                                               rate = Sim_table[i, HB])*365))  
  } else if(dist == "weibull"){
    df <- data.table(Disc_drugA = ceiling(rweibull(n = sample_size, shape=Sim_table[i, k_A], scale=Sim_table[i, lambda_A])*365),
                     Disc_drugB = ceiling(rweibull(n = sample_size, shape=Sim_table[i, k_B], scale=Sim_table[i, lambda_B])*365))
  } else{
    print("only exponential & weibull distributions are supported for now")
  }
  
  df[, DrugA_length := Sim_table[i]$DrugA]
  df[, DrugB_length := Sim_table[i]$DrugB]
  
  # Last dose received
  df[, Last_doseA := DrugA_length * Disc_drugA %/% DrugA_length]
  df[, Last_doseB := DrugB_length * Disc_drugB %/% DrugB_length]
  
  # First dose missed
  df[, First_missedA := Last_doseA + DrugA_length] 
  df[, First_missedB := Last_doseB + DrugB_length]
  
  # Drug
  df[, DrugA := "A"]
  df[, DrugB := "B"]
  
  df[, Registry := rep(LETTERS[1:10], each = sample_size/10)]
  
  tmp_dfA <- df[, .(Disc_drugA, DrugA, Last_doseA, First_missedA, Registry)]
  setnames(tmp_dfA, c("Disc_drugA", "DrugA", "Last_doseA", "First_missedA"),
                    c("Disc_drug", "Drug", "Last_dose", "First_missed"))
  tmp_dfB <- df[, .(Disc_drugB, DrugB, Last_doseB, First_missedB, Registry)]
  setnames(tmp_dfB, c("Disc_drugB", "DrugB", "Last_doseB", "First_missedB"),
           c("Disc_drug", "Drug","Last_dose", "First_missed"))

  df_long <- rbind(tmp_dfA, tmp_dfB)
  
  df_long[,  Drug := as.factor(Drug)]
  
  # Note: For consistency reasons, replace any last dose received times of zero with 0.5.
  df_long[Last_dose == 0, Last_dose := 0.5]
  
  end_data <- Sys.time()
  if (i < nrow(Sim_table)) data_state[[i+1]] <- .Random.seed
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  #                                      Methods 
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # M1 Surv Cox ----------------------------------------------------------------
  # Event of interest is discontinuation decision and is observed exactly.
  
  start_M1_surv_cox <- Sys.time()
  
  fit_survival <- survival::coxph(Surv(time = pmin(Disc_drug, FU), event=Disc_drug <= FU) ~Drug,
                                  data = df_long)
  
  M1_Surv_Cox[i, "Est"] <- fit_survival$coefficients
  M1_Surv_Cox[i, "Var"] <- fit_survival$var
  M1_Surv_Cox[i, "pval"] <- summary(fit_survival)$coefficients[5]
  
  end_M1_surv_cox <- Sys.time()
  M1_Surv_Cox[i, "time"] <- difftime(end_M1_surv_cox, start_M1_surv_cox, units="secs")
  M1_Surv_Cox[i, "iter"] <- Sim_table[i, iter]
  
  print("M1_Surv_Cox complete")
  
  # M2 tram Surv ---------------------------------------------------------------
  # Event of interest is discontinuation decision but it is interval-censored unless dosing interval is <= one day.
  # Note: To guard against occasional RNG-dependent inferior optimization solutions in tram::Coxph for interval-censored data, models were fitted twice and the fit with the higher log-likelihood was retained.
  # Note: Reproducibility: Small numerical differences in results obtained with tram::Coxph() may occur across software environments despite using identical data and random seeds. These differences may arise from changes in tram or its dependencies, including optimization packages such as alabama and nloptr. Exact reproduction therefore requires the same package versions.
  
  start_M2_Surv_Tram <- Sys.time()
  tmp_df_long <- copy(df_long)
  tmp_df_long[, Event := as.numeric(Disc_drug <= FU)]
  tmp_df_long[Event == 0, time := FU]
  tmp_df_long[Event == 1, time := Last_dose]
  tmp_df_long[Event == 0, time2 := Inf]
  tmp_df_long[Event == 1, time2 := pmin(First_missed, FU)] 
  tmp_df_long[, Registry := as.factor(Registry)]
  if(Sim_table[i, DrugA] == 1){
    tmp_df_long[Drug == "A" & Event == 1, time2 := time]
  }
  if(Sim_table[i, DrugB] == 1){
    tmp_df_long[Drug == "B" & Event == 1, time2 := time]
  }
  
  RNG_state[[i]][["before_tram"]] <- list(fit_state[[1]][[i]], fit_state[[2]][[i]])
  RNG_state[[i]][["after_tram"]] <- fit <- vector(mode="list", length=2)
  
  for (j in 1:2) {
    .Random.seed <- fit_state[[j]][[i]]
    
    fit[[j]] <- tram::Coxph(Surv(time=time, time2=time2, type="interval2") ~ Drug, data=tmp_df_long, log_first = T)
    
    RNG_state[[i]][["after_tram"]][[j]] <- .Random.seed
  }
  
  fit_tram_surv <- if (as.numeric(logLik(fit[[1]])) >= as.numeric(logLik(fit[[2]]))) fit[[1]] else fit[[2]]
  
  M2_Surv_Tram[i, "Est"] <- fit_tram_surv$coef["DrugB"]
  M2_Surv_Tram[i, "Var"] <- summary(fit_tram_surv)$test$vcov
  M2_Surv_Tram[i, "pval"] <- summary(fit_tram_surv)$test$test$pvalues
  
  end_M2_Surv_Tram <- Sys.time()
  M2_Surv_Tram[i, "time"] <- difftime(end_M2_Surv_Tram, start_M2_Surv_Tram, units="secs")
  M2_Surv_Tram[i, "iter"] <- Sim_table[i, iter]
  
  print("M2_Surv_Tram complete")
    
  # M3 Surv Cox-----------------------------------------------------------------
  # Event of interest is the last dose received and is observed exactly.
  
  start_M3_surv_cox <- Sys.time()
  
  fit_survival <- survival::coxph(Surv(time = ifelse(test=Disc_drug <= FU, yes=Last_dose, no=FU), event=Disc_drug <= FU) ~
                                    Drug,
                                  data = df_long)
  
  M3_Surv_Cox[i, "Est"] <- fit_survival$coefficients
  M3_Surv_Cox[i, "Var"] <- fit_survival$var
  M3_Surv_Cox[i, "pval"] <- summary(fit_survival)$coefficients[5]
  
  end_M3_surv_cox <- Sys.time()
  M3_Surv_Cox[i, "time"] <- difftime(end_M3_surv_cox, start_M3_surv_cox, units="secs")
  M3_Surv_Cox[i, "iter"] <- Sim_table[i, iter]
  print("M3_Surv_Cox complete")
  
  # M4 Surv Cox ----------------------------------------------------------------
  # Event of interest is the first dose missed and is observed exactly. 
  
  start_M4_surv_cox <- Sys.time()
  
  fit_survival <- survival::coxph(Surv(time = pmin(First_missed, FU), event=First_missed <= FU) ~
                                    Drug,
                                  data = df_long)
  
  M4_Surv_Cox[i, "Est"] <- fit_survival$coefficients
  M4_Surv_Cox[i, "Var"] <- fit_survival$var
  M4_Surv_Cox[i, "pval"] <- summary(fit_survival)$coefficients[5]
  
  end_M4_surv_cox <- Sys.time()
  M4_Surv_Cox[i, "time"] <- difftime(end_M4_surv_cox, start_M4_surv_cox, units="secs")
  M4_Surv_Cox[i, "iter"] <- Sim_table[i, iter]
  print("M4_Surv_Cox complete")
  
  # M5 Surv Cox ----------------------------------------------------------------
  # Event of interest depends on registry and is always observed exactly.
  
  start_M5_Surv_surv <- Sys.time()

  tmp_df_long <- copy(df_long)
  last_dose <- c("A", "B", "C", "D", "E")
  exact_disc <- c("F", "G", "H")
  first_missed <- c("I", "J")
  tmp_df_long[Registry %in% exact_disc, type := "exact_disc"]
  tmp_df_long[Registry %in% last_dose, type := "last_dose"]
  tmp_df_long[Registry %in% first_missed, type := "first_missed"]
  tmp_df_long[type == "exact_disc", Event := Disc_drug <= FU]
  tmp_df_long[type == "last_dose", Event := Disc_drug <= FU]
  tmp_df_long[type == "first_missed", Event := First_missed <= FU]
  tmp_df_long[type == "exact_disc", EFU := pmin(Disc_drug, FU)]
  tmp_df_long[type == "last_dose", EFU := ifelse(test=Disc_drug <= FU, yes=Last_dose, no=FU)]
  tmp_df_long[type == "first_missed", EFU := pmin(First_missed, FU)] 
  
  fit_cox_surv <- survival::coxph(Surv(time = EFU, event=Event) ~ Drug + strata(Registry),
                               data = tmp_df_long)
  
  M5_Surv_Cox[i, "Est"] <- fit_cox_surv$coef["DrugB"]
  M5_Surv_Cox[i, "Var"] <- fit_cox_surv$var
  M5_Surv_Cox[i, "pval"] <- summary(fit_cox_surv)$coefficients[5]
  
  end_M5_Surv_surv <- Sys.time()
  M5_Surv_Cox[i, "time"] <- difftime(end_M5_Surv_surv, start_M5_Surv_surv, units="secs")
  M5_Surv_Cox[i, "iter"] <- Sim_table[i, iter]
  print("M5_Surv_Cox complete")
  
  
  # End methods ----------------------------------------------------------------
}

end_all <- Sys.time()
time_elapsed <- end_all - start_all

setDT(M1_Surv_Cox)
M1_Surv_Cox[, method := "M1_Surv_Cox"]
M1_Surv_Cox <- cbind(M1_Surv_Cox, Sim_table[,-"iter"])

setDT(M2_Surv_Tram)
M2_Surv_Tram[, method := "M2_Surv_Tram"]
M2_Surv_Tram <- cbind(M2_Surv_Tram, Sim_table[,-"iter"])

setDT(M3_Surv_Cox)
M3_Surv_Cox[, method := "M3_Surv_Cox"]
M3_Surv_Cox <- cbind(M3_Surv_Cox, Sim_table[,-"iter"])

setDT(M4_Surv_Cox)
M4_Surv_Cox[, method := "M4_Surv_Cox"]
M4_Surv_Cox <- cbind(M4_Surv_Cox, Sim_table[,-"iter"])

setDT(M5_Surv_Cox)
M5_Surv_Cox[, method := "M5_Surv_Cox"]
M5_Surv_Cox <- cbind(M5_Surv_Cox, Sim_table[,-"iter"])

sim_res <- rbind(M1_Surv_Cox,M2_Surv_Tram,M3_Surv_Cox,M4_Surv_Cox,M5_Surv_Cox)

Info_session <- sessionInfo()
Info_sys <- Sys.info()

save(Info_session, Info_sys, time_elapsed, sim_res, data_state, RNG_state, file = paste0("/home/sharedfolders/study_D7/simulation_results/", paste0("sim_res_", dist, sample_size, "_2000_iter_number_revised.Rdata")))


