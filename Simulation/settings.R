# ******************************************************************************
# Description: Package management, settings/options and paths
# ****************************************************************************** 

source("Simulation/functions.R")

# Packages to be loaded
list_of_packages <- c("data.table",
                      "dplyr",
                      "combinat",
                      "tram",
                      "survival",
                      "tidyr",
                      "GGally",
                      "S7",
                      "ggplot2",
                      "parallel")

cat(paste0("-------- Packages to be loaded selected by the user: \n", 
           paste0(list_of_packages, collapse = "\n"), "\n"))

# Load packages, if a package doesnt exist install it, returns a message
package_management(list_of_packages)

data_path <- ""
meta_path <- ""

# Character variables always loaded as strings
options(stringsAsFactors = FALSE)
options(scipen = 999)

cat("\n-------- Global options loaded --------\n")
