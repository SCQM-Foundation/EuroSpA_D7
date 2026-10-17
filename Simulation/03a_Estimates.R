
# Title: Assessing estimates data -----------------------------------------
# From additional targeted simulations


# Set-up ------------------------------------------------------------------

# rm(list = ls()) #For interactive use
# .rs.restartR() #For interactive use. Or type "R --vanilla" in terminal for a completely fresh R process

library(tidyr)
library(dplyr)
library(ggplot2)


# Loading estimates data --------------------------------------------------

load(file="/home/sharedfolders/study_D7/simulation_results/sim_res_add_targeted_scen.Rdata")


# Quick checks ------------------------------------------------------------

# > Number of failed simulation repetitions: ------------------------------
est_df %>% group_by(DGM) %>% summarize(n_failed_est = sum(is.na(est_logHR)), n_failed_var = sum(is.na(var)), n_failed_pval = sum(is.na(pval)))

  
# Pairs plots of estimates ------------------------------------------------

theme_set(theme_bw(base_family = "Arial")) #for use with ggplot (requires ~/fonts/Arial.ttf)

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

cairo_pdf(file="/home/sharedfolders/study_D7/simulation_results/figures/Paired_estimates_add_targeted_scen.pdf", onefile=TRUE, family="Arial")

for(i in 1:6) {
  
  df <- est_df %>%
    filter(DGM == i) %>%
    select(repetition, model, est_logHR) %>%
    pivot_wider(names_from=model, values_from=est_logHR) %>%
    select(-repetition) %>%
    as.data.frame()
  
  par(oma = c(5, 1, 3, 1))
  
  pairs(df, lower.panel=lowerFn, upper.panel=upperFn, diag.panel=diagFn, main=paste0("DGM ", i), gap=0.5)
  mtext("Estimates", side=1, outer=TRUE, line=3.5)
  
}

dev.off()
rm(df)


# Pairs plots of variance estimates ---------------------------------------

cairo_pdf(file="/home/sharedfolders/study_D7/simulation_results/figures/Paired_variance_estimates_add_targeted_scen.pdf", onefile=TRUE, family="Arial")

for(i in 1:6){
  
  df <- est_df %>%
    filter(DGM == i) %>%
    select(repetition, model, se) %>%
    pivot_wider(names_from=model, values_from=se) %>%
    select(-repetition) %>%
    as.data.frame()
  
  par(oma = c(5, 1, 3, 1))
  
  pairs(df, lower.panel=lowerFn, upper.panel=upperFn, diag.panel=diagFn, main=paste0("DGM ", i), gap=0.5)
  mtext("Square root of variance estimates", side=1, outer=TRUE, line=3.5)
  
}

dev.off()
rm(df)


# QQ-plot -----------------------------------------------------------------

cairo_pdf(file="/home/sharedfolders/study_D7/simulation_results/figures/QQ_plot_add_targeted_scen.pdf", height=8, width=12, family="Arial")

p1 <- ggplot(est_df, aes(sample=est_logHR, color=factor(model))) +
  stat_qq_line() +
  stat_qq(alpha=0.8, size=1) +
  facet_grid(~DGM, scales="free") +
  theme_bw() +
  xlab("Theoretical") +
  ylab("Sample") +
  theme(legend.position="bottom") +
  ggtitle("QQ-Plot ") +
  theme(plot.title=element_text(hjust=0.5)) + theme(text=element_text(size=12)) + guides(color=guide_legend(title="Outcome"))
print(p1)

dev.off()
rm(p1)


# Zip plots ---------------------------------------------------------------

cairo_pdf(file="/home/sharedfolders/study_D7/simulation_results/figures/Zip_add_targeted_scen.pdf", height=10, width=15, onefile=TRUE, family="Arial")

for (i in unique(est_df$model)) {
  
  df <- est_df %>% 
    filter(model == i) %>% 
    mutate(LCL = est_logHR - 1.96*se) %>% 
    mutate(UCL = est_logHR + 1.96*se) %>% 
    mutate(z = abs((est_logHR - logHR) / se)) %>% 
    arrange(DGM, z) %>%
    group_by(DGM) %>% 
    mutate(percent = round(row_number()/nsim*100, digits=2)) %>% 
    mutate(flag_color = logHR < LCL | logHR > UCL) %>% 
    mutate(color = ifelse(test=flag_color, yes="purple", no="cadetblue")) %>% 
    as.data.frame()
  
  title <- paste0(df$model)
  p1 <- ggplot(df, aes(x=est_logHR, y=as.numeric(percent))) + 
    geom_point(color="white") +
    geom_linerange(aes(xmin=LCL, xmax=UCL, color=color)) +
    scale_color_identity() +
    geom_hline(yintercept=95, color="red", linewidth=0.1) +
    facet_grid(~DGM, scales="free") +
    ylab("Fractional centile of |z|") +
    xlab("95% confidence intervals") +
    geom_vline(aes(xintercept=logHR), color ="red", linewidth=0.1) +
    ggtitle(title) + theme(plot.title=element_text(hjust=0.5)) +
    theme(text=element_text(size=14)) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  print(p1)
  
  
}

dev.off()
rm(df, p1)


# Estimates plot ----------------------------------------------------------

cairo_pdf(file="/home/sharedfolders/study_D7/simulation_results/figures/Estimates_add_targeted_scen.pdf", height=10, width=15, family="Arial")

df <- est_df %>% 
  group_by(DGM, model) %>% 
  mutate(mean_group = mean(est_logHR)) %>% 
  mutate(n=row_number()) %>% 
  as.data.frame()

p1 <- ggplot(data=df, aes(x=est_logHR, y=n)) +
  geom_point(alpha=0.2) +
  geom_vline(aes(xintercept=logHR), color="red", linewidth=0.5) +
  geom_vline(aes(xintercept=mean_group), color="yellow", linewidth=0.5) +
  facet_grid(DGM~model) +
  theme_bw()+ 
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank()) + 
  ylab("") +
  xlab("Estimates")
print(p1)

dev.off()
rm(df, p1)


# Variances plot ----------------------------------------------------------

cairo_pdf(file="/home/sharedfolders/study_D7/simulation_results/figures/Variances_add_targeted_scen.pdf", height=10, width=15, family="Arial")

df <- est_df %>% 
  group_by(DGM, model) %>% 
  mutate(mean_group = mean(se)) %>% 
  mutate(n=row_number()) %>% 
  as.data.frame()

p1 <- ggplot(data=df, aes(x=se, y=n)) +
  geom_point(alpha=0.2) +
  geom_vline(aes(xintercept=mean_group), color="yellow", linewidth=0.5) +
  facet_grid(DGM~model) +
  theme_bw()+ 
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank()) + 
  ylab("") +
  xlab("Square root of variance estimates")
print(p1)

dev.off()
rm(df, p1)


# Deriving differences to SM2 estimates -----------------------------------

est_diff_targeted <- est_df %>%
  filter(nobs != 1000) %>% 
  mutate(est_logHR = -est_logHR) %>% #Note: to get 56 vs 7 days.
  tidyr::pivot_wider(id_cols=c(repetition, DGM, dist, nobs), names_from=model, values_from=est_logHR) %>% 
  mutate(across(c(SM3, SM4, SM5), ~ .x - SM2, .names="{.col}_SM2")) %>% 
  as.data.frame()

save(est_diff_targeted, file="Simulation/simulation_results/estimate_differences_add_targeted_scen.Rdata")


