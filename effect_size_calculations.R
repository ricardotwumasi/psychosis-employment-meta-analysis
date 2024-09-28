# ===============================================
# Meta-Analysis Effect Size Calculation Script
# ===============================================
# Description:
# This script calculates effect sizes (Hedges' g) for a meta-analysis 
# on employment attitudes towards people living with psychosis 
# across six studies. It handles different types of effect measures 
# (SMD, RR, PLO) and ensures data validity.
# 
# Date: 28-9-24
# ===============================================

# -------------------------------
# 1. Install and Load Packages
# -------------------------------

# Define required packages
required_packages <- c("metafor", "dplyr", "tidyr", "purrr")

# Install any packages that are not already installed
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed_packages) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load necessary libraries
library(metafor)
library(dplyr)
library(tidyr)
library(purrr)

# -------------------------------
# 2. Define Conversion Function
# -------------------------------

# Function to convert various effect size measures to Hedges' g
# Parameters:
#   yi: Effect size estimate
#   vi: Variance of the effect size
#   measure: Type of measure ("SMD", "RR", "PLO")
#   n1: Sample size of group 1
#   n2: Sample size of group 2 (if applicable)
# Returns:
#   A list containing the converted Hedges' g and its variance
convert_to_hedges_g <- function(yi, vi, measure, n1, n2 = NULL) {
  if (measure == "SMD") {
    # Assuming escalc with measure="SMD" already provides Hedges' g
    return(list(yi = yi, vi = vi))
  } else if (measure == "RR") {
    # Convert log risk ratio to Cohen's d using Chinn (2000) method
    d <- yi * (sqrt(3) / pi)
    var_d <- vi * (3 / pi^2)
  } else if (measure == "PLO") {
    # Convert log odds to Cohen's d using Chinn (2000) method
    d <- yi * (sqrt(3) / pi)
    var_d <- vi * (3 / pi^2)
  } else {
    stop("Unsupported measure type")
  }
  
  # Apply correction factor to get Hedges' g
  n_total <- ifelse(is.null(n2), n1, n1 + n2)
  if (n_total <= 2) {
    stop("Total sample size must be greater than 2 for correction factor.")
  }
  j <- 1 - (3 / (4 * (n_total - 2)))  # Corrected formula
  g <- j * d
  var_g <- j^2 * var_d
  
  return(list(yi = g, vi = var_g))
}

# --------------------------------
# 3. Create Structured Data Frame
# --------------------------------

# Create a data frame with all study data
studies <- tribble(
  ~author,              ~year, ~measure, ~m1,    ~m2,     ~sd1,   ~sd2,   ~n1, ~n2, ~a,  ~b,  ~c,  ~d,  ~xi, ~ni,
  "Andersson et al.",   2015, "SMD",    10.47,  8.23,    3.70,   4.26,   210, 71,  NA,  NA,  NA,  NA,  NA,  NA,
  "Bricout & Bentley", 2000, "SMD",    139.59, 126.15,  17.04,  15.04,  74, 86,  NA,  NA,  NA,  NA,  NA,  NA,
  "Fyhn et al.",        2021, "RR",     NA,     NA,      NA,     NA,     NA, NA, 692, 525, 882, 1074, NA,  NA,
  "Manning & White",    1995, "PLO",    NA,     NA,      NA,     NA,     NA, NA,  NA,  NA,  NA,  NA,  72, 109,
  "Tsang et al.",       2012, "PLO",    NA,     NA,      NA,     NA,     NA, NA,  NA,  NA,  NA,  NA, 161, 183,
  "Zissi et al.",       2007, "PLO",    NA,     NA,      NA,     NA,     NA, NA,  NA,  NA,  NA,  NA, 74,  102
)

# ------------------------------------
# 4. Validate Study Data
# ------------------------------------

# Function to validate study data
# Checks for missing required fields and sufficient sample sizes
# Parameters:
#   study: A single row from the studies data frame
# Returns:
#   Stops execution with an error message if validation fails
validate_study_data <- function(study) {
  if (study$measure == "SMD") {
    required_fields <- c("m1", "m2", "sd1", "sd2", "n1", "n2")
    if (any(is.na(study[required_fields]))) {
      stop(paste("Missing data in study:", study$author))
    }
    if ((study$n1 + study$n2) < 3) {
      stop(paste("Insufficient total sample size in study:", study$author))
    }
  } else if (study$measure == "RR") {
    required_fields <- c("a", "b", "c", "d")
    if (any(is.na(study[required_fields]))) {
      stop(paste("Missing data in study:", study$author))
    }
    if ((study$a + study$b + study$c + study$d) < 3) {
      stop(paste("Insufficient total sample size in study:", study$author))
    }
  } else if (study$measure == "PLO") {
    required_fields <- c("xi", "ni")
    if (any(is.na(study[required_fields]))) {
      stop(paste("Missing data in study:", study$author))
    }
    if (study$ni < 3) {
      stop(paste("Insufficient sample size in study:", study$author))
    }
  } else {
    stop(paste("Unsupported measure type in study:", study$author))
  }
}

# Apply validation to all studies
walk(split(studies, seq(nrow(studies))), validate_study_data)

# -----------------------------------
# 5. Calculate Effect Sizes
# -----------------------------------

# Function to calculate effect size for a single study
# Parameters:
#   study: A single row from the studies data frame (as a list)
# Returns:
#   A data frame with author, year, yi (Hedges' g), vi (variance), and n (sample size)
calculate_effect_size <- function(study) {
  if (study$measure == "SMD") {
    # Calculate standardized mean difference using escalc
    res <- escalc(measure = "SMD", 
                 m1i = study$m1, m2i = study$m2, 
                 sd1i = study$sd1, sd2i = study$sd2, 
                 n1i = study$n1, n2i = study$n2)
    # Convert to Hedges' g (already corrected in escalc with measure="SMD")
    g <- list(yi = res$yi, vi = res$vi)
  } else if (study$measure == "RR") {
    # Calculate risk ratio using escalc
    res <- escalc(measure = "RR", 
                 ai = study$a, bi = study$b, 
                 ci = study$c, di = study$d)
    # Convert log risk ratio to Hedges' g
    g <- convert_to_hedges_g(yi = res$yi, vi = res$vi, 
                            measure = "RR", 
                            n1 = study$a + study$b, 
                            n2 = study$c + study$d)
  } else if (study$measure == "PLO") {
    # Calculate log odds using escalc
    res <- escalc(measure = "PLO", 
                 xi = study$xi, ni = study$ni)
    # Convert log odds to Hedges' g
    g <- convert_to_hedges_g(yi = res$yi, vi = res$vi, 
                            measure = "PLO", 
                            n1 = study$ni)
  } else {
    stop(paste("Unsupported measure type in study:", study$author))
  }
  
  # Determine total sample size for each study
  if (study$measure == "SMD") {
    total_n <- study$n1 + study$n2
  } else if (study$measure == "RR") {
    total_n <- study$a + study$b + study$c + study$d
  } else if (study$measure == "PLO") {
    total_n <- study$ni
  }
  
  # Return a data frame with the results
  return(data.frame(
    author = study$author,
    year = study$year,
    yi = g$yi,
    vi = g$vi,
    n = total_n
  ))
}

# Apply the effect size calculation to all studies using purrr::map
effect_sizes_list <- studies %>%
  split(1:nrow(studies)) %>%
  map(calculate_effect_size)

# Combine all effect sizes into a single data frame
effect_sizes <- bind_rows(effect_sizes_list)

# -----------------------------------
# 6. Display and Save Results
# -----------------------------------

# Print the consolidated effect sizes
print(effect_sizes)

# Define the data directory
data_dir <- "data"

# Create the 'data' directory if it doesn't exist
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

# Define the file path using file.path for cross-platform compatibility
file_path <- file.path(data_dir, "effect_sizes.csv")

# Save the effect sizes to a CSV file
write.csv(effect_sizes, file_path, row.names = FALSE)

# -----------------------------------
# End of Script
# -----------------------------------
