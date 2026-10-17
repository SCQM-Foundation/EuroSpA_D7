# ******************************************************************************
# Description: Convenience and scientific functions
# ****************************************************************************** 

# Input a list of packages and downloads the ones that are missing
package_management <- function(list_of_packages){
  
  # Install missing ones
  new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  
  # Load all required packages
  pkg_management <- unlist(lapply(list_of_packages, require, character.only = TRUE))
  
  if (all(pkg_management)){
    cat("-------- All packages loaded --------")
  } else {
    cat(paste0("Packages missing: ",
               paste0(list_of_packages[!pkg_management], collapse = ", "),
               ". Check dependencies!"))
  }
}

