# ******************************************************************************
# Description: Define simulation settings
# ****************************************************************************** 

# Clean environment (for interactive use)
rm(list = ls())

# Source functions, get packages etc
source("Simulation/settings.R")

# Simulation parameters
reps <- 2000
FU <- 3*365
sample_size <- c(500, 1000)
dist <- c("exponential", "weibull")

# Distributional parameters
baseline_hazard <- -log(0.8) #According to this EuroSpA publication for axSpA, https://ard.bmj.com/content/78/11/1536.long, about 80% of patients retained their TNFi therapy beyond one year. Another EuroSpA publication for PsA, https://pubmed.ncbi.nlm.nih.gov/31665497/, reported a bit lower retention probabilities for one year. ALTERNATIVELY, just say that 80% retention at one year and a median retention of 8 years was inspired by the data used for the application (KM plot for either time to last dose received or first dose missed without differentiating bwt diagnoses or drugs). 

# Weibull parameters (shape and scale) are set such that 1) 80% of patients retain treatment beyond one year and 2) the median retention is eight years.
weibull_k <- log(-log(2)/log(0.8))/log(8) #shape parameter of Weibull distribution, a shape parameter < 1 leads to a decreasing hazard function with time. Kept the same for both treatments.
weibull_lambda <- 1/(-log(0.8))^(1/weibull_k) #scale parameter of Weibull parameter. The scale parameter is changed such that a time-constant HR is obtained. 

# Intervals data frame
intervals <- c(1, 7, 14, 56)
intervals <- expand.grid(intervals, intervals)
setDT(intervals)

# Intervals for hr == 2
intervals_not1 <- replicate(1, intervals, simplify = FALSE)
intervals_not1[[1]]$HR <- 2

# Intervals for hr = 1
intervals1 <- copy(intervals)
intervals1[, row:= .I]
intervals1[, rowtot := sum(Var1, Var2), by = row]
intervals1 <- intervals1[!intervals1[, duplicated(rowtot)]]
intervals1[, c("row", "rowtot") := NULL]
intervals1[, HR := 1]
intervals_not1[[2]]  <- intervals1

# All interval combinations with HR
intervals <- rbindlist(intervals_not1)
intervals[, flag := 1]

# Final set of paramters
params <- copy(intervals)
setDT(params)
params[, flag := NULL]
setnames(params, c("Var1", "Var2"), c("DrugA", "DrugB"))
setorder(params, HR)

# Create simulation table
Sim_table <- as.data.table(do.call("rbind", replicate(reps, params, simplify = FALSE)))

Sim_table[HR == 1, HA := baseline_hazard]
Sim_table[HR == 1, HB := baseline_hazard]
Sim_table[HR == 2, HA := baseline_hazard]
Sim_table[HR == 2, HB := 2*baseline_hazard]

Sim_table[, ':='(k_A = weibull_k, lambda_A = weibull_lambda)]
Sim_table[, ':='(k_B = weibull_k, lambda_B = weibull_lambda/HR^(1/weibull_k))]

Sim_table[, string_para := paste0(DrugA, "_", DrugB, "_", HR)]
setorder(Sim_table, string_para)
Sim_table[, iter := rep(1:reps, length(unique(Sim_table$string_para)))]

# Setting RNG streams
RNGkind("L'Ecuyer-CMRG")
set.seed(23964)

## Four data-generation streams
dg_state <- vector("list", 4)
names(dg_state) <- as.vector(outer(dist, sample_size, paste, sep="_"))
dg_state[[1]] <- .Random.seed

for (i in 2:4) {
  dg_state[[i]] <- parallel::nextRNGStream(dg_state[[i - 1]])
}

## Two main streams for fitting with one fitting substream for each replication in each main stream
fit_state <- vector("list", length=2)

for (i in 1:2) {
  fit_state[[i]] <- vector("list", nrow(Sim_table))
}

fit_state[[1]][[1]] <- parallel::nextRNGStream(dg_state[[4]])
fit_state[[2]][[1]] <- parallel::nextRNGStream(fit_state[[1]][[1]])

for (i in 1:2) {
  for (j in 2:nrow(Sim_table)) {
    fit_state[[i]][[j]] <-
      parallel::nextRNGSubStream(fit_state[[i]][[j - 1]])
  }
}

## Another two data-generation streams
dg_state[["exponential_100"]] <- parallel::nextRNGStream(fit_state[[2]][[1]])
dg_state[["weibull_100"]] <- parallel::nextRNGStream(dg_state[["exponential_100"]])

## One targeted simulation stream
tdg_state <- parallel::nextRNGStream(dg_state[[6]])

save(FU, Sim_table, dg_state, tdg_state, fit_state, file="Simulation/simulation_results/Sim_param.Rdata")
