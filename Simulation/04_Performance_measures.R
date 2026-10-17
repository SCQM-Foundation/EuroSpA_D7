# ******************************************************************************
# Description: Calculate performance measures
# ****************************************************************************** 

rm(list = ls()) # Clean environment (for interactive use)

# Source functions, get packages etc
source("Simulation/settings.R")

# load simulated data
load(file = "Simulation/simulation_results/all_data_revised.Rdata")

# Bias -------------------------------------------------------------------------
bias <- all_data[, .(bias_estimate = sum(Est - log(HR)) / .N,
                     bias_MCSE = sqrt(sum((Est - mean(Est))**2)/(.N*(.N-1)))),
                 by = c("string_para", "method", "distribution", "sample_size")][
                       ,lower_CI_bias := bias_estimate -1.96*bias_MCSE][
                         ,higher_CI_bias := bias_estimate +1.96*bias_MCSE]
bias[, HR:= unlist(lapply(strsplit(string_para,"_"), "[", 3))]
bias[, setting:= paste(unlist(lapply(strsplit(string_para,"_"), "[", 1)),
                       unlist(lapply(strsplit(string_para,"_"), "[", 2)), sep = "_")]
save(bias, file = "Simulation/simulation_results/bias_revised.Rdata")

# EmpSE ------------------------------------------------------------------------
empse <- all_data[, .(Empse_est = sqrt(sum((Est-mean(Est))**2)/(.N-1)),
                          reps = .N),
                      by = c("string_para", "method", "distribution", "sample_size")][
                        , Empse_mcse := Empse_est / sqrt(2*(reps-1))][
                          , lower_CI_empse := Empse_est - 1.96*Empse_mcse][
                            , higher_CI_empse := Empse_est + 1.96*Empse_mcse]
empse[, HR:= unlist(lapply(strsplit(string_para,"_"), "[", 3))]
empse[, setting:= paste(unlist(lapply(strsplit(string_para,"_"), "[", 1)),
                       unlist(lapply(strsplit(string_para,"_"), "[", 2)), sep = "_")]

save(empse, file = "Simulation/simulation_results/empse_revised.Rdata")

# Relative error ---------------------------------------------------------------
relerror  <- all_data[, .(Empse_est = sqrt(sum((Est-mean(Est))**2)/(.N-1)),
                              reps = .N,
                              avg_modSE_estimate = sqrt(mean(Var))),
                          by = c("string_para", "method", "distribution", "sample_size")][
                            , rlt_error_modSE_estimate := ((avg_modSE_estimate / Empse_est)-1)*100]


d1 <- copy(all_data)
d1 <- merge(d1, unique(relerror[, .(string_para, method, distribution, sample_size, avg_modSE_estimate)]),
            all.x = T, by = c("string_para", "method", "distribution", "sample_size"))
d1[, avg_modVar_estimate := avg_modSE_estimate**2]
key <- c("string_para", "method", "distribution", "sample_size")
d1 <- d1[, .(var_var = sum((Var - avg_modVar_estimate)**2) / (.N-1) ), by = key]
relerror <- merge(relerror, d1[, .(var_var, string_para, method, distribution, sample_size)], by = key, all.x = T)
relerror[, avg_modSE_MCSE :=  sqrt(var_var / (4*reps*avg_modSE_estimate**2))]
relerror[, rlt_error_modSE_MCSE := 100*(avg_modSE_estimate / Empse_est)*(sqrt( (var_var/(4*reps*avg_modSE_estimate**4)) + (1/(2*(reps-1)))))]
relerror[, lower_CI_rlt_error_modse := rlt_error_modSE_estimate - 1.96*rlt_error_modSE_MCSE]
relerror[, higher_CI_rlt_error_modse := rlt_error_modSE_estimate + 1.96*rlt_error_modSE_MCSE]
relerror[, HR:= unlist(lapply(strsplit(string_para,"_"), "[", 3))]
relerror[, setting:= paste(unlist(lapply(strsplit(string_para,"_"), "[", 1)),
                        unlist(lapply(strsplit(string_para,"_"), "[", 2)), sep = "_")]

save(relerror, file = "Simulation/simulation_results/relerror_revised.Rdata")

# Coverage ---------------------------------------------------------------------
wald_df <- copy(all_data)
wald_df[, ci_high := Est + 1.96*sqrt(Var)]
wald_df[, ci_low := Est - 1.96*sqrt(Var)]
wald_df[log(HR) >= ci_low & log(HR) <= ci_high, Flag := 1]
coverage <- wald_df[, .(coverage_estimate = sum(Flag, na.rm = T)/.N,
                        reps = .N), by = c("string_para", "method", "distribution", "sample_size")]
coverage[, coverage_MCSE := sqrt((coverage_estimate*(1-coverage_estimate))/reps)]
coverage[, lower_CI_coverage := coverage_estimate -1.96*coverage_MCSE]
coverage[, higher_CI_coverage := coverage_estimate +1.96*coverage_MCSE]
coverage[, HR:= unlist(lapply(strsplit(string_para,"_"), "[", 3))]
coverage[, setting:= paste(unlist(lapply(strsplit(string_para,"_"), "[", 1)),
                           unlist(lapply(strsplit(string_para,"_"), "[", 2)), sep = "_")]

save(coverage, file = "Simulation/simulation_results/coverage_revised.Rdata")

# Precision gain from M1 -------------------------------------------------------
succesful_reps <- reshape(empse[, .(string_para, method, distribution, sample_size, reps)],
                          idvar = c("string_para", "distribution", "sample_size"), timevar = "method", direction = "wide")

df_prec <- reshape(empse[, .(string_para, method, distribution, sample_size, Empse_est)], 
                   idvar = c("string_para", "distribution", "sample_size"), timevar = "method", direction = "wide") 

df_prec[, ':='(Precision_M2 = (((Empse_est.M1_Surv_Cox/Empse_est.M2_Surv_Tram)**2)-1)*100,
               Precision_M3 = (((Empse_est.M1_Surv_Cox/Empse_est.M3_Surv_Cox)**2)-1)*100,
               Precision_M4 = (((Empse_est.M1_Surv_Cox/Empse_est.M4_Surv_Cox)**2)-1)*100,
               Precision_M5 = (((Empse_est.M1_Surv_Cox/Empse_est.M5_Surv_Cox)**2)-1)*100,
               MCSE_M2_1 = 200*(Empse_est.M1_Surv_Cox/Empse_est.M2_Surv_Tram)**2,
               MCSE_M3_1 = 200*(Empse_est.M1_Surv_Cox/Empse_est.M3_Surv_Cox)**2,
               MCSE_M4_1 = 200*(Empse_est.M1_Surv_Cox/Empse_est.M4_Surv_Cox)**2,
               MCSE_M5_1 = 200*(Empse_est.M1_Surv_Cox/Empse_est.M5_Surv_Cox)**2)]

cor_df <- reshape(all_data[, .(Est, method, string_para, distribution, sample_size, iter)],
                  idvar = c("string_para", "distribution", "sample_size", "iter"), timevar = "method", direction = "wide")

corr_state <- cor_df[, .(cor_stateM2 =sqrt((1- cor(Est.M1_Surv_Cox, Est.M2_Surv_Tram)**2)/(.N-1)),
                              cor_stateM3 =sqrt((1- cor(Est.M1_Surv_Cox, Est.M3_Surv_Cox)**2)/(.N-1)),
                              cor_stateM4 =sqrt((1- cor(Est.M1_Surv_Cox, Est.M4_Surv_Cox)**2)/(.N-1)),
                              cor_stateM5 =sqrt((1- cor(Est.M1_Surv_Cox, Est.M5_Surv_Cox)**2)/(.N-1))), by = c("string_para", "distribution", "sample_size")]

mcse_pres <- merge(df_prec, corr_state, by = c("string_para", "distribution", "sample_size"), all.x = T)
mcse_pres[, mcse_presicion2 := MCSE_M2_1*cor_stateM2]
mcse_pres[, mcse_presicion3 := MCSE_M3_1*cor_stateM3]
mcse_pres[, mcse_presicion4 := MCSE_M4_1*cor_stateM4]
mcse_pres[, mcse_presicion5 := MCSE_M5_1*cor_stateM5]

# Estimate
precision_est_df_long <- gather(mcse_pres[, .(string_para, distribution, sample_size, Precision_M2, Precision_M3, Precision_M4, Precision_M5)], method, estimate, Precision_M2:Precision_M5)
setDT(precision_est_df_long)
precision_est_df_long[, method := gsub("Precision_", "", method)]

# MCSE
precision_mcse_df_long <- gather(mcse_pres[, .(string_para, distribution, sample_size, mcse_presicion2, mcse_presicion3, mcse_presicion4, mcse_presicion5)], method, mcse, mcse_presicion2:mcse_presicion5)
setDT(precision_mcse_df_long)
precision_mcse_df_long[, method := paste0("M", gsub("mcse_presicion", "", method))]

# Join estimate with MCSE
precision <- merge(precision_est_df_long, precision_mcse_df_long)
precision[, lower_CI := estimate -1.96*mcse]
precision[, higher_CI := estimate +1.96*mcse]
precision[, HR:= unlist(lapply(strsplit(string_para,"_"), "[", 3))]
precision[, setting:= paste(unlist(lapply(strsplit(string_para,"_"), "[", 1)),
                           unlist(lapply(strsplit(string_para,"_"), "[", 2)), sep = "_")]

precision[, method := case_when(
  method == "M1" ~ "M1_Surv_Cox",
  method == "M2" ~ "M2_Surv_Tram",
  method == "M3" ~ "M3_Surv_Cox",
  method == "M4" ~ "M4_Surv_Cox",
  method == "M5" ~ "M5_Surv_Cox",
)]


save(precision, file = "Simulation/simulation_results/precision_revised.Rdata")

# Power ------------------------------------------------------------------------
df_pow <- copy(all_data)
df_pow[, logHR := log(HR)]
df_pow[, flag_pow := ifelse(pval < 0.05, 1, 0)]

# Estimate
df_pow <- df_pow[, .(power=sum(flag_pow)/.N, nsim_success = .N) , by = c("string_para", "method", "distribution", "sample_size", "HR")]

# MCSE
df_pow[, mcse_power := sqrt(  (power*(1-power))/nsim_success)]
df_pow[, lower_CI := power -1.96*mcse_power]
df_pow[, higher_CI := power +1.96*mcse_power]
df_pow[, HR:= unlist(lapply(strsplit(string_para,"_"), "[", 3))]
df_pow[, setting:= paste(unlist(lapply(strsplit(string_para,"_"), "[", 1)),
                            unlist(lapply(strsplit(string_para,"_"), "[", 2)), sep = "_")]
save(df_pow, file = "Simulation/simulation_results/df_pow_revised.Rdata")
