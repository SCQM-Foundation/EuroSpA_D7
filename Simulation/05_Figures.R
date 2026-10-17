# ******************************************************************************
# Description: Create plots of performance measures
# ****************************************************************************** 

rm(list = ls()) # Clean environment (for interactive use)

# Plots for MS -----

# Get packages -----
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load summarised data -----
load("Simulation/simulation_results/bias_revised.Rdata")
load("Simulation/simulation_results/empse_revised.Rdata")
load("Simulation/simulation_results/relerror_revised.Rdata")
load("Simulation/simulation_results/coverage_revised.Rdata")
load("Simulation/simulation_results/df_pow_revised.Rdata")
load("Simulation/simulation_results/precision_revised.Rdata")

# Setting figure font
theme_set(theme_bw(base_family = "Arial")) #for use with ggplot (requires ~/fonts/Arial.ttf)

# Function to produce plots
plots <- function(data, type = "MS", x, low, high, title, vinter=NULL, is.relative = F, jit){
  data_selected <- copy(data)
  if(type == "MS"){
    data_selected <- data_selected[HR == 1]
    data_selected[, setting := factor(setting, c("1_1", "7_7", "14_14",
                                                 "56_56", "7_1", "14_7",
                                                 "14_1", "56_14", "56_7",
                                                 "56_1"))]  
  } else if (type == "SUP"){
    data_selected <- data_selected[HR == 2]
    data_selected[, setting := factor(setting, c("1_1", "7_7", "14_14", "56_56",
                                                 "7_1", "1_7", "14_7", "7_14",
                                                 "14_1", "1_14", "56_14", "14_56",
                                                 "56_7", "7_56", "56_1", "1_56"))]
  }
  
  data_selected[, method := case_when(
    method == "M1_Surv_Cox" ~ "SM1",
    method == "M2_Surv_Tram" ~ "SM2",
    method == "M3_Surv_Cox" ~ "SM3",
    method == "M4_Surv_Cox" ~ "SM4",
    method == "M5_Surv_Cox" ~ "SM5",
  )]
  
  if(!is.relative){
    p <- ggplot(data_selected, aes(x = setting, y = get(x), shape = method, fill= method))+
      geom_linerange(aes(ymin = get(low), ymax = get(high)),  colour = "black", position = position_dodge(width = jit))
  } else{
    p <- ggplot(data_selected, aes(x = setting, y = get(x)*100, shape = method, fill=method))+
      geom_linerange(aes(ymin = get(low)*100, ymax = get(high)*100),  colour = "black",position = position_dodge(width = jit))
  }
  
  p <- p + geom_point(size = 3,position = position_dodge(width = jit))+
    facet_grid(sample_size~distribution)+
    theme_bw()+
    ylab(title)+
    theme(legend.title=element_blank())+
    theme(legend.position="bottom")+
    xlab("Setting")+ 
    theme(text = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_shape_manual(values=c(SM1=21,SM2 = 22,SM3 = 23,SM4= 24,SM5= 25))+
    scale_fill_manual(values=c(SM1="white",SM2 = "grey88",SM3 = "grey68",SM4= "grey47",SM5= "black"))
  
  if(!is.null(vinter)){
    p <- p+ geom_hline(yintercept = vinter, linetype = "dotted", color = "black", linewidth = 0.7)
  }
    
  if(grepl("bias", x)) {
    p <- p + 
      geom_hline(yintercept=c(-0.05, 0.05), linetype="dotted", color="grey68", linewidth = 0.7) +
      coord_flip(ylim=c(-0.05, 0.05))
  } else (
    p <- p + coord_flip()
  )
  
  return(p)
}


# Save here
path <- "Simulation/figures/"
# Manuscript

plots(data = bias, type = "MS", x = "bias_estimate", low = "lower_CI_bias", high = "higher_CI_bias", title = "Bias with respect to log(HR)", vinter = 0, is.relative = F, jit = 1)
ggsave(paste0(path, "biasHR1_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

plots(data = empse, type = "MS", x = "Empse_est", low = "lower_CI_empse", high = "higher_CI_empse", title = "Empirical Standard Error", is.relative = F, jit = 1)
ggsave(paste0(path, "empseHR1_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

plots(data = coverage, type = "MS", x = "coverage_estimate", low = "lower_CI_coverage", high = "higher_CI_coverage", title = "Coverage %", vinter = 95,is.relative = T, jit = 1)
ggsave(paste0(path, "coverageHR1_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

plots(data = relerror, type = "MS", x = "rlt_error_modSE_estimate", low = "lower_CI_rlt_error_modse", high = "higher_CI_rlt_error_modse", title = "Relative error %", vinter = 0,is.relative = F, jit = 1) # relative calculation already in the table
ggsave(paste0(path, "relerrorHR1_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

plots(data = df_pow, type = "MS", x = "power", low = "lower_CI", high = "higher_CI", title = "Type I error %", vinter = 5, is.relative = T, jit = 0.7) 
ggsave(paste0(path, "type1error_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

plots(data = precision, type = "MS", x = "estimate", low = "lower_CI", high = "higher_CI", title = "Precision relative to decision, exact %", vinter = 0,is.relative = F, jit = 0.7) 
ggsave(paste0(path, "precisionHR1_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

# Supplementary ---------
plots(data = bias, type = "SUP", x = "bias_estimate", low = "lower_CI_bias", high = "higher_CI_bias", title = "Bias with respect to log(HR)", vinter = 0, is.relative = F, jit = 0.7)
ggsave(paste0(path, "biasHR2_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

plots(data = empse, type = "SUP", x = "Empse_est", low = "lower_CI_empse", high = "higher_CI_empse", title = "Empirical Standard Error", is.relative = F, jit = 0.7)
ggsave(paste0(path, "empseHR2_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

plots(data = coverage, type = "SUP", x = "coverage_estimate", low = "lower_CI_coverage", high = "higher_CI_coverage", title = "Coverage %", vinter = 95,is.relative = T, jit = 0.7)
ggsave(paste0(path, "coverageHR2_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

plots(data = relerror, type = "SUP", x = "rlt_error_modSE_estimate", low = "lower_CI_rlt_error_modse", high = "higher_CI_rlt_error_modse", title = "Relative error %", vinter = 0,is.relative = F, jit = 0.7) # relative calcualtion already in the table
ggsave(paste0(path, "relerrorHR2_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

plots(data = df_pow, type = "SUP", x = "power", low = "lower_CI", high = "higher_CI", title = "Power %", is.relative = T, jit = 0.7) 
ggsave(paste0(path, "powerHR2_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

plots(data = precision, type = "SUP", x = "estimate", low = "lower_CI", high = "higher_CI", title = "Precision relative to decision, exact %", vinter = 0,is.relative = F, jit = 0.7) 
ggsave(paste0(path, "precisionHR2_revised.pdf"), height = 10, width = 12, family="Arial", device=cairo_pdf)

