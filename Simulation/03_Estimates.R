# ******************************************************************************
# Description: Combine and plot estimate data
# ****************************************************************************** 

rm(list = ls()) # Clean environment (for interactive use)

# Get packages
library(data.table)
library(dplyr)
library(ggplot2)

# Exponential, 2000 iterations, 500 sample size
load("/home/sharedfolders/study_D7/simulation_results/sim_res_exponential500_2000_iter_number_revised.Rdata")
exp_500 <- copy(sim_res); rm(sim_res, RNG_state, data_state, Info_session, Info_sys, time_elapsed)
exp_500[, distribution := "Exponential"]
exp_500[, sample_size := 500]
exp_500[, c("HA", "HB", "k_A", "lambda_A", "k_B", "lambda_B") := NULL]

# Exponential, 2000 iterations, 1000 sample size
load("/home/sharedfolders/study_D7/simulation_results/sim_res_exponential1000_2000_iter_number_revised.Rdata")
exp_1000 <- copy(sim_res); rm(sim_res, RNG_state, data_state, Info_session, Info_sys, time_elapsed)
exp_1000[, distribution := "Exponential"]
exp_1000[, sample_size := 1000]
exp_1000[, c("HA", "HB", "k_A", "lambda_A", "k_B", "lambda_B") := NULL]

# Weibull, 2000 iterations, 500 sample size
load("/home/sharedfolders/study_D7/simulation_results/sim_res_weibull500_2000_iter_number_revised.Rdata")
weib_500 <- copy(sim_res); rm(sim_res, RNG_state, data_state, Info_session, Info_sys, time_elapsed)
weib_500[, distribution := "Weibull"]
weib_500[, sample_size := 500]
weib_500[, c("HA", "HB", "k_A", "lambda_A", "k_B", "lambda_B") := NULL]

# Weibull, 2000 iterations, 1000 sample size
load("/home/sharedfolders/study_D7/simulation_results/sim_res_weibull1000_2000_iter_number_revised.Rdata")
weib_1000 <- copy(sim_res); rm(sim_res, RNG_state, data_state, Info_session, Info_sys, time_elapsed)
weib_1000[, distribution := "Weibull"]
weib_1000[, sample_size := 1000]
weib_1000[, c("HA", "HB", "k_A", "lambda_A", "k_B", "lambda_B") := NULL]

# Exponential, 2000 iterations, 100 sample size
load("/home/sharedfolders/study_D7/simulation_results/sim_res_exponential100_2000_iter_number_revised.Rdata")
exp_100 <- copy(sim_res); rm(sim_res, RNG_state, data_state, Info_session, Info_sys, time_elapsed)
exp_100[, distribution := "Exponential"]
exp_100[, sample_size := 100]
exp_100[, c("HA", "HB", "k_A", "lambda_A", "k_B", "lambda_B") := NULL]

# Weibull, 2000 iterations, 100 sample size
load("/home/sharedfolders/study_D7/simulation_results/sim_res_weibull100_2000_iter_number_revised.Rdata")
weib_100 <- copy(sim_res); rm(sim_res, RNG_state, data_state, Info_session, Info_sys, time_elapsed)
weib_100[, distribution := "Weibull"]
weib_100[, sample_size := 100]
weib_100[, c("HA", "HB", "k_A", "lambda_A", "k_B", "lambda_B") := NULL]

# Rbind all data into a long format table
all_data <- rbind(exp_100, exp_500, exp_1000, weib_100, weib_500, weib_1000)

# Export such that if there is a need to run 03_Estimates.R all_data.Rdata is up to date.
save(all_data, file = "Simulation/simulation_results/all_data_revised.Rdata")

# Quick checks
summary(all_data[, Est])
all_data[is.na(Est), ]
summary(all_data[, Var])
all_data[is.na(Var), ]


# Pairs plots ------------------------

format_diff <- function(x) {
  if (x == 0) {
    bquote("Max |diff| = " * "0")
  } else {
    if (x < 0.001) {
      exponent <- floor(log10(x))
      mantissa <- x / 10^exponent
      
      if (sprintf("%.1f", mantissa) == "1.0") {
        bquote("Max |diff| = " * 10^.(exponent))
      } else {
        bquote("Max |diff| = " * .(sprintf("%.1f", mantissa)) %*% 10^.(exponent))
      }
    } else {
      bquote("Max |diff| = " * .(sprintf("%.3f", x)))
    }
  }
}

lowerFn <- function(x, y, show.max=TRUE, ...) {
  points(x, y, col=adjustcolor("black", alpha.f=0.5), pch=16, cex=0.6)
  abline(a=0, b=1, col="red")
  
  if (show.max) {
    max_abs_diff <- max(abs(x - y), na.rm=TRUE)
    usr <- par("usr")
    text(x=usr[2], y=usr[3], labels=format_diff(max_abs_diff), adj=c(1.05, -0.5), cex=0.7)
  }
}

upperFn <- function(x, y, ...) {
  usr <- par("usr")
  on.exit(par(usr=usr))
  par(usr=c(0, 1, 0, 1))
  
  r <- cor(x, y, use="complete.obs")
  text(0.5, 0.5, labels=formatC(r, digits=6, format="f"), cex=1.5)
}

diagFn <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr=usr))
  
  x <- x[is.finite(x)]
  d <- density(x)
  y <- d$y / max(d$y) #Note: pairs() supplies the x-range; rescale density to panel coordinates
  
  par(usr=c(usr[1:2], 0, 1.5))
  lines(d$x, y)
}

Loop_tab <- unique(all_data[, .(sample_size, distribution, HR)])

# Setting figure font
theme_set(theme_bw(base_family = "Arial")) #for use with ggplot (requires ~/fonts/Arial.ttf)

i <- 1
# Pair plots
for(i in 1:nrow(Loop_tab)){
  
  # Long to wide format
  M1 <- all_data[method== "M1_Surv_Cox" & HR == Loop_tab$HR[i]  & sample_size == Loop_tab$sample_size[i] & distribution == Loop_tab$distribution[i], .(Est, Var, string_para, sample_size, iter)]
  setnames(M1, c("Est", "Var"), c("Est_M1", "Var_M1"))
  
  M2 <- all_data[method== "M2_Surv_Tram" & HR == Loop_tab$HR[i]  & sample_size == Loop_tab$sample_size[i] & distribution == Loop_tab$distribution[i], .(Est, Var, string_para, sample_size, iter)]
  setnames(M2, c("Est", "Var"), c("Est_M2", "Var_M2"))
  
  M3 <- all_data[method== "M3_Surv_Cox" & HR == Loop_tab$HR[i] & sample_size == Loop_tab$sample_size[i] & distribution == Loop_tab$distribution[i], .(Est, Var, string_para, sample_size, iter)]
  setnames(M3, c("Est", "Var"), c("Est_M3", "Var_M3"))
  
  M4 <- all_data[method== "M4_Surv_Cox"& HR == Loop_tab$HR[i] & sample_size == Loop_tab$sample_size[i] & distribution == Loop_tab$distribution[i], .(Est, Var, string_para, sample_size, iter)]
  setnames(M4, c("Est", "Var"), c("Est_M4", "Var_M4"))
  
  M5 <- all_data[method== "M5_Surv_Cox" & HR == Loop_tab$HR[i] & sample_size == Loop_tab$sample_size[i] & distribution == Loop_tab$distribution[i], .(Est, Var, string_para, sample_size, iter)]
  setnames(M5, c("Est", "Var"), c("Est_M5", "Var_M5"))
  
  results_all <- left_join(M1, M2) %>% left_join(M3)%>% left_join(M4)%>% left_join(M5)
  results_all$intervals <- paste0(unlist(lapply(strsplit(results_all$string_para, "_"), "[", 1)),"_",unlist(lapply(strsplit(results_all$string_para, "_"), "[", 2)))
  results_all$HR <- unlist(lapply(strsplit(results_all$string_para, "_"), "[", 3))
  setting <- paste0("Sample size = ", Loop_tab$sample_size[i],", Distribution = ", Loop_tab$distribution[i], ", HR = ", Loop_tab$HR[i])
  
  results_all <- results_all[order(factor(intervals, levels = c("1_1", "7_7", "14_14",
                                                 "56_56", "1_7", "7_1",
                                                 "7_14","14_7",  "1_14", "14_1",
                                                 "14_56","56_14", 
                                                 "7_56", "56_7", 
                                                 "1_56", "56_1"
                                                 ))), .SD]
  setorder(results_all, HR)
  
  # Estimates
  cairo_pdf(file = paste0("/home/sharedfolders/study_D7/simulation_results/figures/Paired Estimates ", setting, "_revised.pdf"), onefile = TRUE, family="Arial")
  j <- 1
  for (j in 1:length(unique(results_all$string_para))){ 
    interv <- unique(results_all[string_para == unique(results_all$string_para)[j],
                intervals])
    dat_tmp <- results_all[string_para == unique(results_all$string_para)[j],
                           .(Est_M1, Est_M2,Est_M3,Est_M4,Est_M5)]
    colnames(dat_tmp) <- gsub(pattern = "Est_", "S", x = colnames(dat_tmp))
    
    par(oma=c(5, 1, 3, 1))
    
    pairs(dat_tmp, lower.panel=lowerFn, upper.panel=upperFn, diag.panel=diagFn, main=paste0(setting,", Intervals = ", interv), gap=0.5)
    mtext("Estimates", side=1, outer=TRUE, line=3.5)
    
    rm(dat_tmp)
  }
  dev.off()
  
  # Estimates variance
  cairo_pdf(file = paste0("/home/sharedfolders/study_D7/simulation_results/figures/Paired Variances ", setting, "_revised.pdf"), onefile = TRUE, family="Arial")

  for (j in 1:length(unique(results_all$string_para))){
    interv <- unique(results_all[string_para == unique(results_all$string_para)[j],
                                 intervals])
    dat_tmp <- results_all[string_para == unique(results_all$string_para)[j],
                           .(Var_M1, Var_M2, Var_M3, Var_M4, Var_M5)]
    colnames(dat_tmp) <- gsub(pattern = "Var_", "S", x = colnames(dat_tmp))
    dat_tmp <- dat_tmp %>% mutate(across(.cols=everything(), .fns=sqrt))
    
    par(oma=c(5, 1, 3, 1))
    
    pairs(dat_tmp, lower.panel=lowerFn, upper.panel=upperFn, diag.panel=diagFn, main=paste0(setting,", Intervals = ", interv), gap=0.5)
    mtext("Square root of variance estimates", side=1, outer=TRUE, line=3.5)
    
    rm(dat_tmp)
  }
  dev.off()
  
}


# qqplot ------------------------
all_data[, setting := paste0(DrugA,"_", DrugB)]
all_data[, setting := factor(setting, c("1_1", "7_7", "14_14",
                                        "56_56", "1_7", "7_1",
                                        "14_7", "7_14", "1_14", "14_1",
                                        "56_14", "14_56",
                                        "7_56", "56_7", 
                                         "1_56", "56_1"))]

all_data[, method := case_when(
  method == "M1_Surv_Cox" ~"SM1",
  method == "M2_Surv_Tram" ~"SM2",
  method == "M3_Surv_Cox" ~"SM3",
  method == "M4_Surv_Cox" ~"SM4",
  method == "M5_Surv_Cox" ~"SM5"
)]
i <- 1


for (i in 1:nrow(Loop_tab)){
  
  setting <- paste0("Sample size = ", Loop_tab$sample_size[i],", Distribution = ", Loop_tab$distribution[i], ", HR = ", Loop_tab$HR[i])
  cairo_pdf(file = paste0("/home/sharedfolders/study_D7/simulation_results/figures/QQ ", setting, "_revised.pdf"), height = 8,width = 12, onefile = TRUE, family="Arial")
  p <- ggplot(all_data[distribution == Loop_tab$distribution[i] & HR == Loop_tab$HR[i]  & sample_size == Loop_tab$sample_size[i]], aes(sample = Est, color = factor(method))) +
  stat_qq_line()+
  stat_qq(alpha = 0.8, size= 1) +
  facet_grid(~setting, scales ="free")+
  theme_bw()+
  xlab("Theoretical")+
  ylab("Sample")+
  theme(legend.position="bottom") +
  ggtitle(paste0("QQ-Plot ", setting))+
  theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size = 12))+
  guides(color=guide_legend(title = "Outcome"))
  print(p)
  dev.off()
  Sys.sleep(1)
}


# Zip plots --------------------------------------------------------------------
zip <- copy(all_data)
zip[, ci_high := Est + 1.96*sqrt(Var)]
zip[, ci_low := Est - 1.96*sqrt(Var)]
zip[, modSE := sqrt(Var)]
zip[, z := abs((Est - log(HR))/ modSE )]

Loop_tab_zip <- Loop_tab %>% select(sample_size, distribution) %>% distinct

i <- 1
j <- 1
for (i in 1:nrow(Loop_tab_zip)){
  M <- zip[sample_size == Loop_tab_zip$sample_size[i] & distribution == Loop_tab_zip$distribution[i], .(Est, Var, HR, string_para, setting, sample_size, iter, ci_low, ci_high, modSE, z, method)]
  setting <- paste0("Sample size = ", Loop_tab_zip$sample_size[i],", Distribution = ", Loop_tab_zip$distribution[i])
  cairo_pdf(file = paste0("/home/sharedfolders/study_D7/simulation_results/figures/Zip ", setting, "_revised.pdf"),height = 10,width = 15, onefile = TRUE, family="Arial")
  for (j in 1:length(unique(zip$method))){
    M_sub <- M[method == unique(zip$method)[j]]
    setorder(M_sub, z)
    key <- c("string_para")
    M_sub[, nro := 1:.N, by = key]
    M_sub$flag_color <- (log(M_sub$HR) < M_sub$ci_low) | (log(M_sub$HR) > M_sub$ci_high)
    M_sub$color <- ifelse(M_sub$flag_color , "purple", "cadetblue")
    M_sub$percent <- round( (M_sub$nro/2000)*100, 2)
    
    M_text <- unique(M_sub$method)
    title <- paste0(setting,", ",M_text)
    M_sub[, setting_plot:= paste0(setting, "_", HR)] 
    M_sub[, setting_plot := factor(setting_plot, c("1_1_1", "7_7_1", "14_14_1", "56_56_1", "7_1_1", "14_7_1", "14_1_1", "56_14_1", "56_7_1" , "56_1_1",
                                                   "1_1_2", "7_7_2", "14_14_2", "56_56_2","7_1_2", "1_7_2", "7_14_2",
                                                   "14_7_2", "1_14_2", "14_1_2", "56_14_2", "14_56_2", "56_7_2", "7_56_2",
                                                   "1_56_2", "56_1_2"))]
  
    p <- ggplot(M_sub, aes(x=Est, y = as.numeric(percent))) + 
      geom_point(color = "white")+
      geom_linerange(aes(xmin=ci_low, xmax=ci_high,color = color))+
      scale_color_identity()+
      geom_hline(yintercept =  95, color  = "red", linewidth = 0.1)+
      facet_wrap(setting_plot~., nrow =3)+theme_bw() +
      ylab("Fractional centile of |z|")+
      xlab("95% confidence intervals")+
      geom_vline(aes(xintercept = log(HR)), color  = "red", linewidth =0.1)+
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))+
      theme(text = element_text(size = 14))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(p)
    print(paste0("i is at ***", i, "*** and j is at ***", j, "***"))
  }
  dev.off()
}

# Estimates plots --------------------------------------------------------------------
est_data <- copy(all_data)
est_data[, spec_mean:= mean(Est), by = c("method", "HR", "setting", "sample_size", "distribution")]
est_data[, spec_mean_sqrtvar:= mean(sqrt(Var)), by = c("method", "HR", "setting", "sample_size", "distribution")]
i <- 1

for(i in 1:nrow(Loop_tab)){
  setting <- paste0("Sample size = ", Loop_tab$sample_size[i],", Distribution = ", Loop_tab$distribution[i], ", HR = ", Loop_tab$HR[i])
  cairo_pdf(file = paste0("/home/sharedfolders/study_D7/simulation_results/figures/Estimates ", setting, "_revised.pdf"), height = 10,width = 15, onefile = TRUE, family="Arial")
  
  tmp <- est_data[sample_size == Loop_tab$sample_size[i] & distribution == Loop_tab$distribution[i] & HR == Loop_tab$HR[i]]
  key <- c("setting", "method")
  tmp[, nro := 1:.N, by =key]
  p1 <- ggplot(data=tmp, aes(x =Est, y = nro))+
    geom_point(alpha=0.2)+
    geom_vline(aes(xintercept = log(Loop_tab$HR[i])), color  = "red", linewidth =0.5)+
    geom_vline(aes(xintercept = spec_mean), color  = "yellow", linewidth =0.5)+
    facet_grid(setting~method)+
    theme_bw()+ 
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) + 
    ylab("")+
    xlab("Estimates")+
    ggtitle(setting)+
    theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size = 12))
  print(p1)
  dev.off()
  Sys.sleep(1)
}

# Variance plots --------------------------------------------------------------------
i <- 1

for(i in 1:nrow(Loop_tab)){
  setting <- paste0("Sample size = ", Loop_tab$sample_size[i],", Distribution = ", Loop_tab$distribution[i], ", HR = ", Loop_tab$HR[i])
  cairo_pdf(file = paste0("/home/sharedfolders/study_D7/simulation_results/figures/Estimates Variances ", setting, "_revised.pdf"), height = 10,width = 15, onefile = TRUE, family="Arial")
  
  tmp <- est_data[sample_size == Loop_tab$sample_size[i] & distribution == Loop_tab$distribution[i] & HR == Loop_tab$HR[i]]
  key <- c("setting", "method")
  tmp[, nro := 1:.N, by =key]
  p1 <- ggplot(data=tmp, aes(x =sqrt(Var), y = nro))+
    geom_point(alpha=0.2)+
    geom_vline(aes(xintercept = spec_mean_sqrtvar), color  = "yellow", linewidth =0.5)+
    facet_grid(setting~method)+
    theme_bw()+ 
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) + 
    ylab("")+
    xlab("Square root of variance estimates")+
    ggtitle(setting)+
    theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size = 12))
  print(p1)
  dev.off()
  Sys.sleep(1)
}


# Deriving differences to SM2 estimates -----------------------------------
# For certain DGMs

est_diff <- all_data %>%
  filter(sample_size != 1000 & ((DrugA == 14 & DrugB == 14) | (DrugA == 56 & DrugB== 7) | (DrugA == 7 & DrugB == 56))) %>% 
  mutate(Est = case_when( 
    (DrugA == 14 & DrugB == 14) | (DrugA == 7 & DrugB == 56) ~ Est,
    DrugA == 56 & DrugB== 7 ~ -Est #Note: to get 56 vs 7 days.
  )) %>% 
  tidyr::pivot_wider(id_cols=c(iter, DrugA, DrugB, HR, string_para, distribution, sample_size), names_from=method, values_from=Est) %>% 
  mutate(across(c(SM3, SM4, SM5), ~ .x - SM2, .names="{.col}_SM2")) %>% 
  as.data.frame()

save(est_diff, file="Simulation/simulation_results/estimate_differences.Rdata")
