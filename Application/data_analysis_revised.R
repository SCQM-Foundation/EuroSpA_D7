
# Title: Analysing registry data ------------------------------------------

# Initializing ------------------------------------------------------------

# rm(list = ls()) #For interactive use
# .rs.restartR() #For interactive use. Or type "R --vanilla" in terminal for a completely fresh R process

source(file="Application/initialise.R") 

library(dplyr)
library(survival)
library(tram)
library(nnet)
library(lubridate)
library(ggplot2)
library(cobalt) 
library(purrr) 
library(ggbeeswarm)

diagnoses <- c("psa", "spa")
contrasts <- c("infliximab vs etanercept", "adalimumab vs certolizumab")
reg_sets <- c("Any registry", "LDR / FDM registries", "Single registry")


# Loading prepared registry data ------------------------------------------

load(file="Application/data/prepared_data.RData")


# Loading differences from simulations ------------------------------------

load(file="Simulation/simulation_results/estimate_differences.Rdata") #Note: loads estimate differences wrt SM2

load(file="Simulation/simulation_results/estimate_differences_add_targeted_scen.Rdata") #Note: loads estimate differences wrt SM2


# Figure file path --------------------------------------------------------

fig_path <- "Application/figures/"


# Setting figure font -----------------------------------------------------

theme_set(theme_bw(base_family = "Arial")) #for use with ggplot


# Analysis data sets ------------------------------------------------------
# Note: Below code creates eight analysis data sets (combo of diagnosis, treatment contrast, and registry set) in d7_final_list.
# Note: All character covariates need to be factors for tram::Coxph() to work.

d7_final_set <- d7_final %>% 
  filter(is.na(withdrawal_reason)) %>% #Note: to disregard patients excluded from the analyses because of missing or incomplete information.
  as.data.frame()

table(d7_final_set$country, d7_final_set$bio_drug_name, d7_final_set$diagnosis_group) #Note: provides information on TNFi data by diagnosis and registry to define selset. Some registries have zero numbers for certain TNFi or very low numbers.

selset <- data.frame(diagnosis=rep(diagnoses, each=2), contrast=rep(contrasts, times=2)) %>% 
  mutate(full_registry_set = case_when(
    diagnosis == "psa" & contrast == "infliximab vs etanercept" ~ "CH, ES, IS, NO, PT, RO, SI",
    diagnosis == "psa" & contrast == "adalimumab vs certolizumab" ~ "CH, ES, NO, PT, SI",
    diagnosis == "spa" & contrast == "infliximab vs etanercept" ~ "CH, ES, NO, PT, RO",
    diagnosis == "spa" & contrast == "adalimumab vs certolizumab" ~ "CH, ES, NO, PT, RO, UK"
  )) %>% 
  as.data.frame()
# Note: EE (an FDM registry as IS) did not meet the data requirement for inclusion (at least 5 patients for both treatments) for any analysis data set.

d7_final_list <- bal_list <- vector(mode="list", length=length(diagnoses))
names(d7_final_list) <- names(bal_list) <- diagnoses

for (diag in diagnoses) {
  
  d7_final_list[[diag]] <- bal_list[[diag]] <- vector(mode="list", length=length(contrasts))
  names(d7_final_list[[diag]]) <- names(bal_list[[diag]]) <- contrasts
  
  for (c in contrasts) {
    
    d7_final_list[[diag]][[c]] <- bal_list[[diag]][[c]] <- vector(mode="list", length=length(reg_sets))
    names(d7_final_list[[diag]][[c]]) <- names(bal_list[[diag]][[c]]) <- reg_sets
    
    for (s in reg_sets) {
      
      d7_final_list[[diag]][[c]][[s]] <- d7_final_set %>% 
        filter(diagnosis_group == diag & bio_drug_name %in% strsplit(c, "\\s*vs\\s*")[[1]] & country %in% strsplit(selset %>% filter(diagnosis == diag & contrast == c) %>% pull(full_registry_set), ",\\s*")[[1]]) %>% 
        filter(if (s == "Any registry") TRUE else if (s == "LDR / FDM registries") bio_drug_stop_date_definition != "decision" else country == "PT") %>% 
        mutate(bio_drug_name = factor(bio_drug_name, levels=rev(strsplit(c, "\\s*vs\\s*")[[1]]))) %>% 
        mutate(gender = factor(gender)) %>% 
        mutate(calendar_year = factor(calendar_year)) %>% 
        mutate(country = factor(country)) %>% 
        mutate(concom_csdmard_baseline = factor(concom_csdmard_baseline)) %>% 
        as.data.frame()
      
    }
  }
}


# IPTW approach for confounding -------------------------------------------
# Note: IPTW are determined separately for each of the eight analysis data sets contained in d7_final_list

bal_df <- trt_freq <- data.frame()

for (diag in diagnoses) {
  cat(diag, ": \n")
  for (c in contrasts) {
    cat(c, ": \n")
    for (s in reg_sets) {
      cat(s, ": \n")
      
      if (s != "Single registry") {
      ps_fit <- glm(formula=bio_drug_name ~ country + gender + age_drug_baseline + calendar_year + concom_csdmard_baseline, data=d7_final_list[[diag]][[c]][[s]], family=binomial())
      } else {
        ps_fit <- glm(formula=bio_drug_name ~ gender + age_drug_baseline + calendar_year + concom_csdmard_baseline, data=d7_final_list[[diag]][[c]][[s]], family=binomial())
      }
      
      ps_1 <- predict(object=ps_fit, type="response")
      cat("Summary of estimated propensities: ", summary(ps_1), "\n")
      
      ps_obs <- ifelse(test=d7_final_list[[diag]][[c]][[s]]$bio_drug_name == levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[2], yes=ps_1, no=1 - ps_1)
      names(ps_obs) <- ifelse(test=d7_final_list[[diag]][[c]][[s]]$bio_drug_name == levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[2], yes=levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[2], no=levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[1])
      cat("Summary of actual propensities (acc. to treatment received): ", levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[1], ": ", summary(ps_obs[names(ps_obs) == levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[1]]), "\n")
      cat("Proportion of unique propensities outside the range observed for ", levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[1], ": ", sum(unique(ps_1) > max(1-ps_obs[names(ps_obs) == levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[1]]) | unique(ps_1) < min(1-ps_obs[names(ps_obs) == levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[1]])) / length(unique(ps_1)), "\n")
      cat("Summary of actual propensities (acc. to treatment received): ", levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[2], ": ", summary(ps_obs[names(ps_obs) == levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[2]]), "\n")
      cat("Proportion of unique propensities outside the range observed for ", levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[2], ": ", sum(unique(ps_1) > max(ps_obs[names(ps_obs) == levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[2]]) | unique(ps_1) < min(ps_obs[names(ps_obs) == levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[2]])) / length(unique(ps_1)), "\n")
      
      p_treat_obs <- prop.table(x=d7_final_list[[diag]][[c]][[s]] %>% pull(bio_drug_name) %>% table)
      cat("Marginal treatment frequencies: ", levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[1], ": ", p_treat_obs[levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[1]], "\n")
      cat("Marginal treatment frequencies: ", levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[2], ": ", p_treat_obs[levels(d7_final_list[[diag]][[c]][[s]]$bio_drug_name)[2]], "\n")
      
      trt_freq <- rbind(
        trt_freq,
        as.data.frame(p_treat_obs) %>% rename(bio_drug_name='.') %>% mutate(diagnosis=diag, contrast=c, regset=s) %>% as.data.frame()
      )
      
      treat_index <- as.numeric(d7_final_list[[diag]][[c]][[s]] %>% pull(bio_drug_name))
      w_stab <- p_treat_obs[treat_index] / ps_obs
      cat("Mean of stabilized weights: ", mean(w_stab), "\n") #Should be close to one.
      
      cat("Order of rows preserved: ", all.equal(current=model.frame(ps_fit), target=d7_final_list[[diag]][[c]][[s]] %>% select(names(model.frame(ps_fit))), check.attributes=F), "\n")
      
      d7_final_list[[diag]][[c]][[s]] <- d7_final_list[[diag]][[c]][[s]] %>% 
        mutate(w_stab = as.vector(w_stab)) %>% 
        as.data.frame()
      
      bal_list[[diag]][[c]][[s]] <- cobalt::bal.tab(x=reformulate(termlabels=c(if (s != "Single registry") "country", "gender", "age_drug_baseline", "calendar_year", "concom_csdmard_baseline"), response="bio_drug_name"), data=d7_final_list[[diag]][[c]][[s]], weights=d7_final_list[[diag]][[c]][[s]]$w_stab, method="weighting", estimand="ATE", s.d.denom="pooled", binary="std", abs=TRUE, quick=F, un=T)
      
      bal_df <- rbind(
        bal_df,
        bal_list[[diag]][[c]][[s]]$Balance %>% mutate(covariate=row.names(.), diagnosis=diag, contrast=c, regset=s) %>% as.data.frame()
      )
      
      print(cobalt::love.plot(bal_list[[diag]][[c]][[s]], abs=TRUE, threshold=0.1, title=paste(diag, c, s, sep="_")))
      
      rm(ps_fit, ps_1, treat_index, ps_obs, p_treat_obs, w_stab)
    }
  }
}

unlist(lapply(X=bal_list, FUN=function(d) lapply(X=d, FUN=function(c) lapply(X=c, FUN=function(s) max(s[["Balance"]][["Diff.Adj"]])))))

lapply(X=d7_final_list, FUN=function(d) lapply(X=d, FUN=function(c) lapply(X=c, FUN=function(s) summary(s$w_stab))))

lapply(X=d7_final_list, FUN=function(d) lapply(X=d, FUN=function(c) lapply(X=c, FUN=function(s) sum(s$w_stab > 10))))

lapply(X=d7_final_list, FUN=function(d) lapply(X=d, FUN=function(c) lapply(X=c, FUN=function(s) table(s %>% filter(w_stab > 10) %>% pull(bio_drug_name)))))

lapply(X=d7_final_list, FUN=function(d) lapply(X=d, FUN=function(c) lapply(X=c, FUN=function(s) table(s %>% filter(w_stab > 10) %>% pull(country)))))

lapply(X=d7_final_list, FUN=function(d) lapply(X=d, FUN=function(c) lapply(X=c, FUN=function(s) by(s$w_stab, s$bio_drug_name, summary))))


# KM plots by event and diagnosis -----------------------------------------

for (diag in diagnoses) {
  for (e in c("last", "first")) {
    ds <- bind_rows(
      d7_final_list[[diag]][["infliximab vs etanercept"]][["LDR / FDM registries"]],
      d7_final_list[[diag]][["adalimumab vs certolizumab"]][["LDR / FDM registries"]]
      )
    
    cairo_pdf(file=paste0(fig_path, "KM_plots_", diag, "_", e, ".pdf"), family="Arial")
    par(las=1, mar=c(10, 5, 2, 2)) 
    
    if (e == "last") {
      sf <- survfit(formula=Surv(time=last_received_time_left, event=!is.na(bio_drug_stop_date)) ~ bio_drug_name, data=ds)
    } else {
      sf <- survfit(formula=Surv(time=first_missed_time_left, event=!(is.na(bio_drug_stop_date) | (!is.na(bio_drug_stop_date) & date_first_dose_missed > date_of_last_visit))) ~ bio_drug_name, data=ds)
    } 
    
    risk_times <- c(0, 365, 2*365, 3*365, 4*365, 5*365, 6*365, 7*365)
    s <- summary(sf, times=risk_times, extend=TRUE)
    nrisk <- matrix(s$n.risk, nrow=length(sf$strata), byrow=TRUE)
    
    xticks <- seq(from=0, to=as.numeric(ceiling_date(as.Date(max(ds$data_cut_date)), unit="months") - as.Date(study_start)), by=182)
    xticks <- c(xticks, last(xticks)+182)
    
    plot(sf, col=1:5, xaxt="n", xlab="Time since TNFi start (months)", ylab="Estimated probability of retaining TNFi", mark.time=T, mark="/")
    axis(side=1, at=xticks, labels=seq(from=0, by=6, length.out=length(xticks)))
    
    usr <- par("usr")
    ypos <- usr[3] - c(0.23, 0.27, 0.31, 0.35) * diff(usr[3:4])
    
    for (i in seq_len(nrow(nrisk))) {
      text(x=risk_times, y=ypos[i], labels=nrisk[i, ], xpd=NA, col=i, cex=0.75)
    }
    
    labels <- sub(pattern="bio_drug_name=", replacement="", x=names(sf$strata))
    labels <- sub(pattern="certolizumab", replacement="CZP", x=labels)
    
    text(x=usr[1], y=ypos, labels=labels, pos=2, xpd=NA, col=1:4, cex=0.7)
    
    text(x=usr[1], y=usr[3] - 0.19 * diff(usr[3:4]), labels="Number at risk:", pos=2, xpd=NA, cex=0.7)
    
    legend(x="topright", legend=labels, col=1:4, lty=1, bty="n")
    
    dev.off()
  }
}


# AM2 - decision, interval-censored approach ------------------------------

fit_AM2_IPTW_untrimmed <- list()
fit_AM2_ca <- list()
fit_AM2 <- list()

for (diag in diagnoses) {
  for (c in contrasts) {
    for (s in reg_sets) {
      
      set.seed(1722008)
      
      fit_AM2_IPTW_untrimmed[[diag]][[c]][[s]] <- tram::Coxph(formula=reformulate(termlabels="bio_drug_name", response=paste0("Surv(time=decision_time_left, time2=decision_time_right, type='interval2')", if (s != "Single registry") " | 0 + country")), data=d7_final_list[[diag]][[c]][[s]], weights=w_stab, cluster=patient_id, log_first=TRUE)
      
      set.seed(1722008)
      
      fit_AM2_ca[[diag]][[c]][[s]] <- tram::Coxph(formula=reformulate(termlabels=c("bio_drug_name", "gender", "age_drug_baseline", "calendar_year", "concom_csdmard_baseline"), response=paste0("Surv(time=decision_time_left, time2=decision_time_right, type='interval2')", if (s != "Single registry") " | 0 + country")), data=d7_final_list[[diag]][[c]][[s]], log_first=TRUE)
      
      set.seed(1722008)
      
      fit_AM2[[diag]][[c]][[s]] <- tram::Coxph(formula=reformulate(termlabels="bio_drug_name", response=paste0("Surv(time=decision_time_left, time2=decision_time_right, type='interval2')", if (s != "Single registry") " | 0 + country")), data=d7_final_list[[diag]][[c]][[s]], log_first=TRUE)
      
    }
  }
}


# AM3 - last dose received, exact approach --------------------------------

fit_AM3_survival_IPTW_untrimmed <- list()
fit_AM3_survival_ca <- list()
fit_AM3_survival <- list()

for (diag in diagnoses) {
  for (c in contrasts) {
    for (s in reg_sets) {
      
      if (s == "Any registry") {
        fit_AM3_survival_IPTW_untrimmed[[diag]][[c]][[s]] <- fit_AM3_survival_ca[[diag]][[c]][[s]] <- fit_AM3_survival[[diag]][[c]][[s]] <- NULL
      } else {
        
        fit_AM3_survival_IPTW_untrimmed[[diag]][[c]][[s]] <- survival::coxph(formula=Surv(time=last_received_time_left, event=as.numeric(!is.na(bio_drug_stop_date))) ~ bio_drug_name + strata(country) + cluster(patient_id), data=d7_final_list[[diag]][[c]][[s]], weights=w_stab)
        
        fit_AM3_survival_ca[[diag]][[c]][[s]] <- survival::coxph(formula=Surv(time=last_received_time_left, event=as.numeric(!is.na(bio_drug_stop_date))) ~ bio_drug_name + gender + age_drug_baseline + calendar_year + concom_csdmard_baseline + strata(country), data=d7_final_list[[diag]][[c]][[s]])
        
        fit_AM3_survival[[diag]][[c]][[s]] <- survival::coxph(formula=Surv(time=last_received_time_left, event=as.numeric(!is.na(bio_drug_stop_date))) ~ bio_drug_name + strata(country), data=d7_final_list[[diag]][[c]][[s]])
        
      }
    }
  }
}


# AM4 - first dose missed, exact approach ----------------------------------

fit_AM4_survival_IPTW_untrimmed <- list()
fit_AM4_survival_ca <- list()
fit_AM4_survival <- list()

for (diag in diagnoses) {
  for (c in contrasts) {
    for (s in reg_sets) {
      
      d7_final_list[[diag]][[c]][[s]] %>% filter(!is.na(bio_drug_stop_date) & date_first_dose_missed > date_of_last_visit) %>% group_by(bio_drug_stop_date_definition, bio_drug_name) %>% summarise(n=n()) %>% print
      
      if (s == "Any registry") {
        fit_AM4_survival_IPTW_untrimmed[[diag]][[c]][[s]] <- fit_AM4_survival_ca[[diag]][[c]][[s]] <- fit_AM4_survival[[diag]][[c]][[s]] <- NULL
      } else {
        
        fit_AM4_survival_IPTW_untrimmed[[diag]][[c]][[s]] <- survival::coxph(formula=Surv(time=first_missed_time_left, event=as.numeric(!(is.na(bio_drug_stop_date) | (!is.na(bio_drug_stop_date) & date_first_dose_missed > date_of_last_visit)))) ~ bio_drug_name + strata(country) + cluster(patient_id), data=d7_final_list[[diag]][[c]][[s]], weights=w_stab)
        
        fit_AM4_survival_ca[[diag]][[c]][[s]] <- survival::coxph(formula=Surv(time=first_missed_time_left, event=as.numeric(!(is.na(bio_drug_stop_date) | (!is.na(bio_drug_stop_date) & date_first_dose_missed > date_of_last_visit)))) ~ bio_drug_name + gender + age_drug_baseline + calendar_year + concom_csdmard_baseline + strata(country), data=d7_final_list[[diag]][[c]][[s]])
        
        fit_AM4_survival[[diag]][[c]][[s]] <- survival::coxph(formula=Surv(time=first_missed_time_left, event=as.numeric(!(is.na(bio_drug_stop_date) | (!is.na(bio_drug_stop_date) & date_first_dose_missed > date_of_last_visit)))) ~ bio_drug_name + strata(country), data=d7_final_list[[diag]][[c]][[s]])
        
      }
    }
  }
}


# AM5 - mixture, exact approach --------------------------------------------

fit_AM5_survival_IPTW_untrimmed <- list()
fit_AM5_survival_ca <- list()
fit_AM5_survival <- list()

for (diag in diagnoses) {
  for (c in contrasts) {
    for (s in reg_sets) {
      
      if (s == "Single registry") {
        fit_AM5_survival_IPTW_untrimmed[[diag]][[c]][[s]] <- fit_AM5_survival_ca[[diag]][[c]][[s]] <- fit_AM5_survival[[diag]][[c]][[s]] <- NULL
      } else {
        
        df <- d7_final_list[[diag]][[c]][[s]] %>% 
          mutate(event_AM5 = case_when(
            bio_drug_stop_date_definition == "first dose missed" ~ as.numeric(!(is.na(bio_drug_stop_date) | (!is.na(bio_drug_stop_date) & date_first_dose_missed > date_of_last_visit))), #Note: actually, there are no !is.na(bio_drug_stop_date) with date_first_dose_missed > date_of_last_visit.
            TRUE ~ as.numeric(!is.na(bio_drug_stop_date))
          )) %>% 
          as.data.frame()
        
        fit_AM5_survival_IPTW_untrimmed[[diag]][[c]][[s]] <- survival::coxph(formula=Surv(time=mixed_time_left, event=event_AM5) ~ bio_drug_name + strata(country) + cluster(patient_id), data=df, weights=w_stab)
        
        fit_AM5_survival_ca[[diag]][[c]][[s]] <- survival::coxph(formula=Surv(time=mixed_time_left, event=event_AM5) ~ bio_drug_name + gender + age_drug_baseline + calendar_year + concom_csdmard_baseline + strata(country), data=df)
        
        fit_AM5_survival[[diag]][[c]][[s]] <- survival::coxph(formula=Surv(time=mixed_time_left, event=event_AM5) ~ bio_drug_name + strata(country), data=df)
        
        rm(df)
      }
    }
  }
}


# Plotting CIs ------------------------------------------------------------

forest_tables <- function(fit, name.fit, contrast, reg.set){ 
  # Applied to each fit object of interest
  # Note re vcov.tram(): Argument with_baseline must be TRUE to obtain robust covariances (needed in case of weighting). This is not clear from the help file. For vcov.coxph with_baseline is not an argument and is simply ignored => we can use the same call for both model types.
  # For tram::Coxph, vcov(fit, with_baseline=TRUE) currently returns a cluster-robust sandwich covariance matrix when the model was fit with the cluster argument. This is based on the implementation in tram version 0.8-3 and may change in future versions. If tram is updated, this behaviour should be re-checked.
  
  est <- coef(fit)[1]
  var_mat <- vcov(fit, with_baseline=T)
  var <- var_mat[rownames(var_mat)[grepl(pattern="bio_drug_name", x=rownames(var_mat))], colnames(var_mat)[grepl(pattern="bio_drug_name", x=colnames(var_mat))]]
  
  tmp_ci <- data.frame(LCL=numeric(1), UCL=numeric(1), est=numeric(1))
  
  tmp_ci[1, ] <- c(exp(est + c(-1, 1)*1.96*sqrt(var)), exp(est))
  
  tmp_ci <- tmp_ci %>% 
    mutate(contrast_of_interest = contrast) %>%
    mutate(reg_set = reg.set) %>% 
    mutate(model = substr(x=name.fit, start=5, stop=7)) %>% 
    as.data.frame()
  
  return(tmp_ci)
}

# > Unweighted, marginal: -------------------------------------------------
analysis_list <- list(fit_AM2, fit_AM3_survival, fit_AM4_survival, fit_AM5_survival)
names(analysis_list) <- c("fit_AM2", "fit_AM3_survival", "fit_AM4_survival", "fit_AM5_survival")

ci_plot_tab <- data.frame()
for (i in 1:length(analysis_list)) {
  for (diag in diagnoses) {
    for (c in contrasts) {
      for (s in reg_sets) {
        if (!is.null(analysis_list[[i]][[diag]][[c]][[s]]) & s != "Single registry") {
          ci_plot_tab <- rbind(
            ci_plot_tab, 
            forest_tables(fit=analysis_list[[i]][[diag]][[c]][[s]], name.fit=names(analysis_list)[i], contrast=c, reg.set=s) %>% 
              mutate(indication = if (diag == "psa") "PsA" else "axSpA") %>% 
              mutate(method = "Unweighted, marginal") %>% 
              as.data.frame()
          )
        }
      }
    }
  }
}

# > IPTW, marginal: ------------------------------------------------------
analysis_list_IPTW_untrimmed <- list(fit_AM2_IPTW_untrimmed, fit_AM3_survival_IPTW_untrimmed, fit_AM4_survival_IPTW_untrimmed, fit_AM5_survival_IPTW_untrimmed)
names(analysis_list_IPTW_untrimmed) <- c("fit_AM2_IPTW_untrimmed", "fit_AM3_survival_IPTW_untrimmed", "fit_AM4_survival_IPTW_untrimmed", "fit_AM5_survival_IPTW_untrimmed")

for (i in 1:length(analysis_list_IPTW_untrimmed)) {
  for (diag in diagnoses) {
    for (c in contrasts) {
      for (s in reg_sets) {
        if (!is.null(analysis_list_IPTW_untrimmed[[i]][[diag]][[c]][[s]]) & s != "Single registry") {
          ci_plot_tab <- rbind(
            ci_plot_tab, 
            forest_tables(fit=analysis_list_IPTW_untrimmed[[i]][[diag]][[c]][[s]], name.fit=names(analysis_list_IPTW_untrimmed)[i], contrast=c, reg.set=s) %>% 
              mutate(indication = if (diag == "psa") "PsA" else "axSpA") %>% 
              mutate(method = "IPTW, marginal") %>% 
              as.data.frame()
          )
        }
      }
    }
  }
}

# > Covariate adjustment: -------------------------------------------------
analysis_list_ca <- list(fit_AM2_ca, fit_AM3_survival_ca, fit_AM4_survival_ca, fit_AM5_survival_ca)
names(analysis_list_ca) <- c("fit_AM2_ca", "fit_AM3_survival_ca", "fit_AM4_survival_ca", "fit_AM5_survival_ca")

for (i in 1:length(analysis_list_ca)) {
  for (diag in diagnoses) {
    for (c in contrasts) {
      for (s in reg_sets) {
        if (!is.null(analysis_list_ca[[i]][[diag]][[c]][[s]]) & s != "Single registry") {
          ci_plot_tab <- rbind(
            ci_plot_tab, 
            forest_tables(fit=analysis_list_ca[[i]][[diag]][[c]][[s]], name.fit=names(analysis_list_ca)[i], contrast=c, reg.set=s) %>% 
              mutate(indication = if (diag == "psa") "PsA" else "axSpA") %>% 
              mutate(method = "Unweighted, covariate-adjusted") %>% 
              as.data.frame()
          )
        }
      }
    }
  }
}


# > Combined plot: --------------------------------------------------------

ci_plot_tab <- ci_plot_tab %>% 
  mutate(method = factor(method, levels=c("Unweighted, marginal", "IPTW, marginal", "Unweighted, covariate-adjusted"))) %>% 
  mutate(contrast_of_interest = if_else(condition=contrast_of_interest == "adalimumab vs certolizumab", true="adalimumab vs CZP", false=contrast_of_interest)) %>% 
  as.data.frame()

p <- ggplot(data=ci_plot_tab, aes(x=contrast_of_interest, y=est, shape=model, fill=model)) +
  geom_linerange(aes(ymin=LCL, ymax=UCL), position=position_dodge(width=0.8)) +
  geom_point(size=3, position=position_dodge(width=0.8)) +
  facet_grid(rows=vars(reg_set, method), cols=vars(indication)) +
  geom_hline(yintercept=1, linetype="dotted", linewidth=0.7) +
  coord_flip() +
  theme_bw() +
  labs(x="Comparison", y="Hazard ratio", shape=NULL, fill=NULL) +
  theme(text=element_text(size=16), legend.position="bottom") +
  scale_shape_manual(values=c(AM2=22, AM3=23, AM4=24, AM5=25)) +
  scale_fill_manual(values=c(AM2="grey88", AM3="grey68", AM4="grey47", AM5="black")) 

ggsave(filename=paste0(fig_path, "/CI_plot_combined.pdf"), plot=p, height=18, width=12, family="Arial", device=cairo_pdf)


# Assessing inflation of standard errors under IPTW -----------------------
# Relative to the inflation expected from the reduction in effective sample size (ESS) implied by the stabilized weights.

obs_tab <- purrr::imap(bal_list, function(lvl1, diag) {
  purrr::imap(lvl1, function(lvl2, contrast) {
    purrr::imap(lvl2, function(df, r) {
      if (r != "Single registry") {
        t(df[["Observations"]]) %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var="bio_drug_name") %>%
          dplyr::mutate(indication = if (diag == "psa") "PsA" else "axSpA", reg_set = r, contrast_of_interest = contrast, .before=1)
      }
    }) %>%
      purrr::list_rbind()
  }) %>%
    purrr::list_rbind()
}) %>%
  purrr::list_rbind()

obs_tab <- obs_tab %>%
  tidyr::pivot_wider(names_from=bio_drug_name, values_from=c("Unadjusted", "Adjusted")) %>%
  as.data.frame()

assessing_IPTW_fit <- left_join(
  x=ci_plot_tab %>% filter(method == "IPTW, marginal") %>% mutate(UCL_LCL_ratio_IPTW = UCL/LCL) %>% select(-c(LCL, UCL, est, method)), 
  y=ci_plot_tab %>% filter(method == "Unweighted, marginal") %>% mutate(UCL_LCL_ratio_uw = UCL/LCL) %>% select(-c(LCL, UCL, est, method))) %>% 
  mutate(se_IPTW = log(UCL_LCL_ratio_IPTW)/(2*1.96)) %>% 
  mutate(se_uw = log(UCL_LCL_ratio_uw)/(2*1.96)) %>%
  mutate(r_se = se_IPTW / se_uw) %>% 
  
  left_join(y=obs_tab %>% mutate(contrast_of_interest = if_else(condition=contrast_of_interest == "infliximab vs etanercept", true=contrast_of_interest, false="adalimumab vs CZP"))) %>% 
  
  mutate(across(.cols=starts_with(match="Unadjusted"), .fns=~ get(sub(pattern="^Unadjusted", replacement="Adjusted", x=cur_column())) / .x, .names="ratio_{.col}")) %>%
  
  mutate(r_se_expected = case_when( #Note: this is an ad-hoc diagnostic based on the two-sample mean case with pooled variance estimation.
    contrast_of_interest == "infliximab vs etanercept" ~ sqrt((Adjusted_infliximab^-1 + Adjusted_etanercept^-1) / (Unadjusted_infliximab^-1 + Unadjusted_etanercept^-1)),
    contrast_of_interest == "adalimumab vs CZP" ~ sqrt((Adjusted_adalimumab^-1 + Adjusted_certolizumab^-1) / (Unadjusted_adalimumab^-1 + Unadjusted_certolizumab^-1))
  )) %>% 
  
  as.data.frame()

View(assessing_IPTW_fit %>% filter(r_se > r_se_expected)) 
sum(assessing_IPTW_fit$r_se > assessing_IPTW_fit$r_se_expected)


# Differences wrt reference outcome model (AM2) ---------------------------

appl_distances <- list()

for (diag in diagnoses) {
  for (c in contrasts) {
    for (r in reg_sets) {
      
      if (r != "Single registry") {
        appl_distances[[diag]][[r]][[c]] <- data.frame()
        
        for (a in c("analysis_list", "analysis_list_ca", "analysis_list_IPTW_untrimmed")) {
          
          res <- get(a)
          
          appl_distances[[diag]][[r]][[c]] <- rbind(
            appl_distances[[diag]][[r]][[c]],
            rbind(
              data.frame(outcome="AM2", est=coef(res[grepl("AM2", names(res))][[1]][[diag]][[c]][[r]])[1], row.names=NULL),
              if (r == "LDR / FDM registries") data.frame(outcome="AM3", est=coef(res[grepl("AM3", names(res))][[1]][[diag]][[c]][[r]])[1], row.names=NULL) else NULL,
              if (r == "LDR / FDM registries") data.frame(outcome="AM4", est=coef(res[grepl("AM4", names(res))][[1]][[diag]][[c]][[r]])[1], row.names=NULL) else NULL,
              data.frame(outcome="AM5", est=coef(res[grepl("AM5", names(res))][[1]][[diag]][[c]][[r]])[1], row.names=NULL)
            ) %>%
              tidyr::pivot_wider(names_from=outcome, values_from=est) %>% 
              mutate(AM3_AM2 = if (r == "LDR / FDM registries") AM3 - AM2 else NULL,
                     AM4_AM2 = if (r == "LDR / FDM registries") AM4 - AM2 else NULL,
                     AM5_AM2 = AM5 - AM2) %>% 
              mutate(analytic_method = if (a == "analysis_list") "Unweighted, marginal" else if (a == "analysis_list_ca") "Unweighted, covariate-adjusted" else "IPTW, marginal") %>% 
              as.data.frame()
            
          )
        }
      }
    }
  }
}


# Differences wrt SM2 from simulation -------------------------------------
sim_plot_distances <- est_diff %>% 
  mutate(setting = case_when(
    string_para %in% c("14_14_1", "14_14_2") & sample_size == 100 ~ "100_14_14",
    string_para %in% c("14_14_1", "14_14_2") & sample_size == 500 ~ "500_14_14",
    string_para %in% c("56_7_1", "56_7_2", "7_56_2") & sample_size == 100 ~ "100_7_56",
    string_para %in% c("56_7_1", "56_7_2", "7_56_2") & sample_size == 500 ~ "500_7_56",
  )) %>% 
  select(setting, SM3_SM2, SM4_SM2, SM5_SM2) %>% 
  as.data.frame()

sim_plot_distances <- rbind(
  sim_plot_distances, 
  est_diff_targeted %>% 
    mutate(setting = case_when(
      nobs == 100 ~ "100_7_56_variable", 
      nobs == 500 ~ "500_7_56_variable"
    )) %>% 
    select(setting, SM3_SM2, SM4_SM2, SM5_SM2)) %>% 
  mutate(setting = factor(setting, levels=unique(setting))) %>%
  as.data.frame()


# Compatibility checks ----------------------------------------------------

# > AM3 vs AM2: -----------------------------------------------------------
appl_subset_32_IE <- purrr::imap_dfr(.x=appl_distances, .f=function(diag_list, diag) {
  diag_list[["LDR / FDM registries"]][["infliximab vs etanercept"]] %>%
    mutate(diagnosis=if (diag == "psa") "PsA" else "axSpA", registry_set="LDR / FDM registries", contrast="infliximab vs etanercept") %>%
    select(diagnosis, registry_set, contrast, analytic_method, AM3_AM2) %>% 
    as.data.frame()
}
)

appl_subset_32_IE <- bind_rows(
  lapply(c("100_7_56", "500_7_56", "100_7_56_variable", "500_7_56_variable"), function(z) {
    appl_subset_32_IE %>% mutate(setting = z) %>% as.data.frame()
  })
)

appl_subset_32_AC <- purrr::imap_dfr(.x=appl_distances, .f=function(diag_list, diag) {
  diag_list[["LDR / FDM registries"]][["adalimumab vs certolizumab"]] %>%
    mutate(diagnosis=if (diag == "psa") "PsA" else "axSpA", registry_set="LDR / FDM registries", contrast="adalimumab vs certolizumab") %>%
    select(diagnosis, registry_set, contrast, analytic_method, AM3_AM2) %>% 
    as.data.frame()
}
)

appl_subset_32_AC <- bind_rows(
  lapply(c("100_14_14", "500_14_14"), function(z) {
    appl_subset_32_AC %>% mutate(setting = z) %>% as.data.frame()
  })
)

p <- ggplot(sim_plot_distances, aes(x=setting, y=SM3_SM2)) +
  coord_cartesian(ylim=c(-ceiling(1000*max(abs(sim_plot_distances$SM3_SM2)))/1000, ceiling(1000*max(abs(sim_plot_distances$SM3_SM2)))/1000)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey") + 
  ggbeeswarm::geom_quasirandom(width=0.3, size=0.8, alpha=0.35, shape=16, colour="grey40") +
  geom_point(data=appl_subset_32_IE, aes(x=setting, y=AM3_AM2, shape=diagnosis, fill=analytic_method), size=3, stroke=0.8, position=position_jitter(width=0.08, height=0)) +
  geom_point(data=appl_subset_32_AC, aes(x=setting, y=AM3_AM2, shape=diagnosis, fill=analytic_method), size=3, stroke=0.8, position=position_jitter(width=0.08, height=0)) +
  labs(x="Setting", y="SM3-SM2 difference", shape="Diagnosis", fill="Analytical method") +
  scale_shape_manual(values=c("PsA"=21, "axSpA"=24)) +
  scale_fill_manual(values=c("Unweighted, marginal"="white", "Unweighted, covariate-adjusted"="grey70", "IPTW, marginal"="grey30")) +
  guides(shape=guide_legend(order=1, override.aes=list(size=4, colour="black", fill="white")), fill=guide_legend(order=2, override.aes=list(shape=21, size=4, colour="black"))) +
  theme_bw()

ggsave(filename=paste0(fig_path, "/compatibility_differences_32.pdf"), plot=p, height=10, width=12, family="Arial", device=cairo_pdf)

# > AM4 vs AM2: -----------------------------------------------------------
appl_subset_42_IE <- purrr::imap_dfr(.x=appl_distances, .f=function(diag_list, diag) {
  diag_list[["LDR / FDM registries"]][["infliximab vs etanercept"]] %>%
    mutate(diagnosis=if (diag == "psa") "PsA" else "axSpA", registry_set="LDR / FDM registries", contrast="infliximab vs etanercept") %>%
    select(diagnosis, registry_set, contrast, analytic_method, AM4_AM2) %>% 
    as.data.frame()
}
)

appl_subset_42_IE <- bind_rows(
  lapply(c("100_7_56", "500_7_56", "100_7_56_variable", "500_7_56_variable"), function(z) {
    appl_subset_42_IE %>% mutate(setting = z) %>% as.data.frame()
  })
)

appl_subset_42_AC <- purrr::imap_dfr(.x=appl_distances, .f=function(diag_list, diag) {
  diag_list[["LDR / FDM registries"]][["adalimumab vs certolizumab"]] %>%
    mutate(diagnosis=if (diag == "psa") "PsA" else "axSpA", registry_set="LDR / FDM registries", contrast="adalimumab vs certolizumab") %>%
    select(diagnosis, registry_set, contrast, analytic_method, AM4_AM2) %>% 
    as.data.frame()
}
)

appl_subset_42_AC <- bind_rows(
  lapply(c("100_14_14", "500_14_14"), function(z) {
    appl_subset_42_AC %>% mutate(setting = z) %>% as.data.frame()
  })
)

p <- ggplot(sim_plot_distances, aes(x=setting, y=SM4_SM2)) +
  coord_cartesian(ylim=c(-ceiling(100*max(abs(sim_plot_distances$SM4_SM2)))/100, ceiling(100*max(abs(sim_plot_distances$SM4_SM2)))/100)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey") + 
  ggbeeswarm::geom_quasirandom(width=0.3, size=0.8, alpha=0.35, shape=16, colour="grey40") +
  geom_point(data=appl_subset_42_IE, aes(x=setting, y=AM4_AM2, shape=diagnosis, fill=analytic_method), size=3, stroke=0.8, position=position_jitter(width=0.08, height=0)) +
  geom_point(data=appl_subset_42_AC, aes(x=setting, y=AM4_AM2, shape=diagnosis, fill=analytic_method), size=3, stroke=0.8, position=position_jitter(width=0.08, height=0)) +
  labs(x="Setting", y="SM4-SM2 difference", shape="Diagnosis", fill="Analytical method") +
  scale_shape_manual(values=c("PsA"=21, "axSpA"=24)) +
  scale_fill_manual(values=c("Unweighted, marginal"="white", "Unweighted, covariate-adjusted"="grey70", "IPTW, marginal"="grey30")) +
  guides(shape=guide_legend(order=1, override.aes=list(size=4, colour="black", fill="white")), fill=guide_legend(order=2, override.aes=list(shape=21, size=4, colour="black"))) +
  theme_bw()

ggsave(filename=paste0(fig_path, "/compatibility_differences_42.pdf"), plot=p, height=10, width=12, family="Arial", device=cairo_pdf)

# > AM5 vs AM2: -----------------------------------------------------------
appl_subset_52_IE <- purrr::imap_dfr(.x=appl_distances, .f=function(diag_list, diag) {
  purrr::imap_dfr(.x=diag_list, .f=function(reg_list, registry) {
    reg_list[["infliximab vs etanercept"]] %>%
      mutate(diagnosis=if (diag == "psa") "PsA" else "axSpA", registry_set=registry, contrast="infliximab vs etanercept") %>%
      select(diagnosis, registry_set, contrast, analytic_method, AM5_AM2)
  }
  )
}
)

appl_subset_52_IE <- bind_rows(
  lapply(c("100_7_56", "500_7_56", "100_7_56_variable", "500_7_56_variable"), function(z) {
    appl_subset_52_IE %>% mutate(setting = z) %>% as.data.frame()
  })
)

appl_subset_52_AC <- purrr::imap_dfr(.x=appl_distances, .f=function(diag_list, diag) {
  purrr::imap_dfr(.x=diag_list, .f=function(reg_list, registry) {
    reg_list[["adalimumab vs certolizumab"]] %>%
      mutate(diagnosis=if (diag == "psa") "PsA" else "axSpA", registry_set=registry, contrast="adalimumab vs certolizumab") %>%
      select(diagnosis, registry_set, contrast, analytic_method, AM5_AM2)
  }
  )
}
)

appl_subset_52_AC <- bind_rows(
  lapply(c("100_14_14", "500_14_14"), function(z) {
    appl_subset_52_AC %>% mutate(setting = z) %>% as.data.frame()
  })
)

p <- ggplot(sim_plot_distances, aes(x=setting, y=SM5_SM2)) +
  coord_cartesian(ylim=c(-ceiling(100*max(abs(sim_plot_distances$SM5_SM2)))/100, ceiling(100*max(abs(sim_plot_distances$SM5_SM2)))/100)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey") + 
  ggbeeswarm::geom_quasirandom(width=0.3, size=0.8, alpha=0.35, shape=16, colour="grey40") +
  geom_point(data=appl_subset_52_IE, aes(x=setting, y=AM5_AM2, shape=diagnosis, colour=registry_set, fill=analytic_method), size=3, stroke=0.8, position=position_jitter(width=0.08, height=0)) +
  geom_point(data=appl_subset_52_AC, aes(x=setting, y=AM5_AM2, shape=diagnosis, colour=registry_set, fill=analytic_method), size=3, stroke=0.8, position=position_jitter(width=0.08, height=0)) +
  labs(x="Setting", y="SM5-SM2 difference", shape="Diagnosis", colour="Registry set", fill="Analytical method") +
  scale_shape_manual(values=c("PsA"=21, "axSpA"=24)) +
  scale_colour_manual(values=c("Any registry"="black", "LDR / FDM registries"="grey40")) +
  scale_fill_manual(values=c("Unweighted, marginal"="white", "Unweighted, covariate-adjusted"="grey70", "IPTW, marginal"="grey30")) +
  guides(shape=guide_legend(order=1, override.aes=list(size=4, colour="black", fill="white")), colour=guide_legend(order=2, override.aes=list(shape=21, size=4, fill="white")), fill=guide_legend(order=3, override.aes=list(shape=21, size=4, colour="black"))) +
  theme_bw()

ggsave(filename=paste0(fig_path, "/compatibility_differences_52.pdf"), plot=p, height=10, width=12, family="Arial", device=cairo_pdf)


# Saving objects for use in .Rmd file -------------------------------------

save(selset, d7_final_set, d7_final_list, bal_list, bal_df, analysis_list_ca, analysis_list_IPTW_untrimmed, analysis_list, assessing_IPTW_fit, appl_distances, file="Application/data/analysis_output.RData")


# Creating requirements file ----------------------------------------------
requirements(script_name="data_analysis_revised")

