
# Title: Plotting estimated performance measures --------------------------
# From additional targeted simulations


# Set-up ------------------------------------------------------------------

# rm(list = ls()) #For interactive use
# .rs.restartR() #For interactive use. Or type "R --vanilla" in terminal for a completely fresh R process

library(dplyr)
library(ggplot2)


# Loading pm data ---------------------------------------------------------

# > Additional targeted simulations: --------------------------------------
load(file="Simulation/simulation_results/pm_data_add_targeted_scen.Rdata")

# > Original simulations: -------------------------------------------------
for (i in c("bias", "empse", "relerror", "coverage", "df_pow", "precision")) {
  
  load(file=paste0("Simulation/simulation_results/", i, "_revised.Rdata"))
  df <- get(i) %>% 
    filter(string_para == "56_7_1") %>% 
    as.data.frame()
  
  assign(x=i, value=df)
  rm(df)
  
}


# Plotting instructions ---------------------------------------------------

theme_set(theme_bw(base_family = "Arial")) #for use with ggplot (requires ~/fonts/Arial.ttf)

plots <- function(data_new, data_old, x_new, x_old, lcl_old, ucl_old, multiplier=1, title, vinter=NULL, jit, filename){
  df_new <- data_new %>% 
    filter(stat == x_new) %>% 
    mutate(distribution = ifelse(test=dist == "exp", yes="Exponential", no="Weibull")) %>% 
    mutate(dosing = "variable") %>% 
    select(distribution, nobs, dosing, model, est, CI_LL, CI_UL) %>% 
    as.data.frame()
  
  if (!is.null(data_old)) {
  df_old <- data_old %>% 
    mutate(model = case_when(
      method == "M1_Surv_Cox" ~ "SM1",
      method == "M2_Surv_Tram" ~ "SM2",
      method == "M3_Surv_Cox" ~ "SM3",
      method == "M4_Surv_Cox" ~ "SM4",
      method == "M5_Surv_Cox" ~ "SM5"
    )) %>% 
    mutate(dosing = "fixed") %>% 
    select(distribution, nobs=sample_size, dosing, model, est=all_of(x_old), CI_LL=all_of(lcl_old), CI_UL=all_of(ucl_old)) %>% 
    as.data.frame()
}
  
  df <- if (!is.null(data_old)) rbind(df_new, df_old) else df_new
  
  if (x_new == "relprec") {
    df <- df %>% filter(model != "SM1") %>% as.data.frame()
  }
  
    p <- ggplot(df, aes(x=dosing, y=est*multiplier, shape=model, fill=model)) +
      geom_linerange(aes(ymin=CI_LL*multiplier, ymax=CI_UL*multiplier), colour="black", position=position_dodge(width=jit))
  
  p <- p + geom_point(size=3, position=position_dodge(width=jit)) +
    facet_grid(nobs~distribution) +
    theme_bw() +
    ylab(title) +
    theme(legend.title=element_blank()) +
    theme(legend.position="bottom") +
    xlab("Dosing interval")+ 
    theme(text=element_text(size=16),
          axis.text.x=element_text(angle=45, hjust=1)) +
    scale_shape_manual(values=c(SM1=21, SM1_meta=21, SM2=22, SM3=23, SM4=24, SM5=25, SM5_meta=25)) +
    scale_fill_manual(values=c(SM1="white", SM1_meta="white", SM2="grey88", SM3="grey68", SM4="grey47", SM5="black", SM5_meta="white"))
  
  if (!is.null(vinter)) {
    p <- p + geom_hline(yintercept=vinter, linetype="dotted", color="black", linewidth=0.7)
  }
  
  if (grepl("bias", x_new)) {
    p <- p + 
      geom_hline(yintercept=c(-0.05, 0.05), linetype="dotted", color="grey68", linewidth=0.7) +
      coord_flip(ylim=c(-0.05, 0.05))
  } else (
    p <- p + coord_flip()
  )
  
  ggsave(filename=filename, height=10, width=12, family="Arial", device=cairo_pdf)
}

path <- "Simulation/figures/"


# Bias plot ---------------------------------------------------------------

plots(data_new=pm_dat$summ, data_old=bias, x_new="bias", x_old="bias_estimate", lcl_old="lower_CI_bias", ucl_old="higher_CI_bias", title="Bias with respect to log(HR)", vinter=0, jit=1, filename=paste0(path, "bias_add_targeted_scen.pdf"))


# Empirical variance plot -------------------------------------------------

plots(data_new=pm_dat$summ, data_old=empse, x_new="empse", x_old="Empse_est", lcl_old="lower_CI_empse", ucl_old="higher_CI_empse", title="Empirical Standard Error", jit=1, filename=paste0(path, "empse_add_targeted_scen.pdf"))


# Model variance plot -------------------------------------------------

plots(data_new=pm_dat$summ, data_old=NULL, x_new="modelse", title="Model-based Standard Error", jit=1, filename=paste0(path, "modse_add_targeted_scen.pdf"))


# Coverage plot -----------------------------------------------------------

plots(data_new=pm_dat$summ, data_old=coverage, x_new="cover", x_old="coverage_estimate", lcl_old="lower_CI_coverage", ucl_old="higher_CI_coverage", multiplier=100, title="Coverage %", vinter=95, jit=1, filename=paste0(path, "coverage_add_targeted_scen.pdf"))


# Relative error plot -----------------------------------------------------

plots(data_new=pm_dat$summ, data_old=relerror, x_new="relerror", x_old="rlt_error_modSE_estimate", lcl_old="lower_CI_rlt_error_modse", ucl_old="higher_CI_rlt_error_modse", title="Relative error %", vinter=0, jit=1, filename=paste0(path, "relerror_add_targeted_scen.pdf"))


# Type I error plot -------------------------------------------------------

plots(data_new=pm_dat$summ, data_old=df_pow, x_new="power", x_old="power", lcl_old="lower_CI", ucl_old="higher_CI", multiplier=100, title="Type I error %", vinter=5, jit=0.7, filename=paste0(path, "type1error_add_targeted_scen.pdf")) 


# Relative precision plot -------------------------------------------------

plots(data_new=pm_dat$summ, data_old=precision, x_new="relprec", x_old="estimate", lcl_old="lower_CI", ucl_old="higher_CI", title="Precision relative to decision, exact %", vinter=0, jit=0.7, filename=paste0(path, "precision_add_targeted_scen.pdf")) 

