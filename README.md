# Meta-Analysis of Employer Attitudes Towards Individuals with Psychosis

This repository contains R code for replicating our meta-analysis on employer discrimination towards individuals with psychosis (Crestois et al., in preparation). The analysis includes an overall meta-analysis, a subgroup analysis using proportion data, and effect size calculations for individual studies.

This code may be repurposed for other meta-analyses or updated to add studies as they are published.

## Contributors

- Ninon Crestois
- Ricardo Twumasi

## Requirements 

- R (version 4.0.0 or higher recommended)
- R packages: `metafor`

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/ricardotwumasi/psychosis-employment-meta-analysis.git
   ```
2. Install the required R package:
   ```R
   install.packages("metafor")
   ```

## Usage

The analysis is split into three main scripts:

1. `effect_size_calculations.R`: Calculates effect sizes for individual studies
2. `meta_analysis.R`: Conducts the overall meta-analysis
3. `subgroup_analysis.R`: Performs a subgroup analysis using proportion data


Run each script separately in R

## Code Explanation

### 1. Effect Size Calculations (effect_size_calculations.R)

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
# End of Effect Size Script
# -----------------------------------

# Load required packages
library(metafor)
library(dplyr)

# Read the effect size data
data <- read.csv("data/effect_sizes.csv")

# Ensure effect size (yi) and variance (vi) are numeric
data$yi <- as.numeric(data$yi)
data$vi <- as.numeric(data$vi)

# Run the meta-analysis
meta_analysis <- rma(yi = yi, vi = vi, data = data, method = "REML", test = "knha")

# View the results
summary(meta_analysis)

# Create forest plot
forest(meta_analysis, slab = data$author, 
       xlab = "Hedges' g", 
       mlab = "Random-Effects Model",
       addpred = TRUE,
       xlim = c(-2, 3),
       alim = c(-2, 2),
       refline = 0,
       header = "Study")

# Create a Baujat plot
baujat(meta_analysis)

# Funnel plot with trim-and-fill
funnel(tf_meta_analysis)

# Calculate Fail-Safe N (Rosenthal's method)
fsn <- fsn(yi, vi, data = data, type = "Rosenthal")
print(fsn)

# Calculate Orwin's Fail-Safe N
# Define the target effect size (e.g., 0.2 for a small effect)
target_effect <- 0.2

# Calculate Orwin's Fail-Safe N
orwin_fsn <- fsn(yi, vi, data = data, type = "Orwin", target = target_effect)
print(orwin_fsn)

### 3. Subgroup Analysis (subgroup_analysis.R)
# Load necessary library
library(metafor)

# Create data frame with study data
study_data <- data.frame(
  Author = c("Zissi et al.", "Manning et al.", "Tsang et al."),
  Successes = c(74, 72, 161),
  Total = c(102, 109, 183),
  stringsAsFactors = FALSE
)

# Calculate proportions
study_data$Proportion <- study_data$Successes / study_data$Total

# Check for proportions exactly 0 or 1 to avoid issues with logit transformation
if(any(study_data$Proportion == 0 | study_data$Proportion == 1)) {
  stop("Proportions of 0 or 1 detected. Consider adding a continuity correction.")
}

# Calculate variance for raw proportions
study_data$Variance <- (study_data$Proportion * (1 - study_data$Proportion)) / study_data$Total

# Truncate study titles for consistency (if necessary)
study_data$TruncatedTitle <- c(
  "Zissi et al.",
  "Manning & White",
  "Tsang et al."
)

# Meta-analysis using raw proportions (Proportion Method)
meta_raw <- rma(yi = Proportion, vi = Variance, data = study_data, method = "REML")

# Summary of the raw proportions meta-analysis
summary(meta_raw)

# Enhanced Forest Plot for raw proportions without transformation
par(mar = c(4, 8, 4, 2))  # Adjust margins: bottom, left, top, right

forest(meta_raw, 
       slab = study_data$TruncatedTitle, 
       xlab = "Proportion",
       ilab = study_data$Total, # Add Total as an additional column
       ilab.xpos = max(meta_raw$ci.ub, na.rm = TRUE) * 0.6, # Positioning for ilab
       cex = 0.8,               # Text size
       psize = 1,               # Point size
       header = "Study",
       mlab = "Random-Effects Model",
       addfit = TRUE,           # Add summary effect
       showweights = TRUE       # Show weights
)

# Add a legend for the additional information
op <- par(cex=0.8, font=2)
text(x = max(meta_raw$ci.ub, na.rm = TRUE) * 0.6, 
     y = length(study_data$Author) + 2, 
     labels = "n", pos = 3)
par(op)

# Reset plotting parameters to default
par(mar = c(5, 4, 4, 2) + 0.1)

# Baujat plot to identify influential studies
baujat(meta_raw)

# Logit transformation of proportions
study_data$LogitProportion <- log(study_data$Proportion / (1 - study_data$Proportion))

# Calculate variance of the logit-transformed proportions
# Var(logit(p)) = Var(p) / [p(1 - p)]^2
study_data$LogitVariance <- study_data$Variance / (study_data$Proportion * (1 - study_data$Proportion))^2

# Meta-analysis using logit-transformed proportions
meta_logit <- rma(yi = LogitProportion, vi = LogitVariance, data = study_data, method = "REML")

# Summary of the logit-transformed proportions meta-analysis
summary(meta_logit)

# Transform the summary effect back to proportions
transform_proportion <- function(logit) {
  exp(logit) / (1 + exp(logit))
}

# Extract the transformed mean proportion and its confidence interval
logit_summary <- summary(meta_logit)
transformed_mean <- transform_proportion(logit_summary$beta)
transformed_ci_lower <- transform_proportion(logit_summary$ci.lb)
transformed_ci_upper <- transform_proportion(logit_summary$ci.ub)

# Display the transformed results
cat("Transformed Mean Proportion:", round(transformed_mean, 3), "\n")
cat("95% CI:", round(transformed_ci_lower, 3), "-", round(transformed_ci_upper, 3), "\n")

# Optional: Enhanced Forest Plot for logit-transformed proportions
# This requires back-transforming the estimates and variances
# Note: metafor does not natively support back-transformed forest plots,
# so it's often better to interpret on the logit scale or use other visualization tools.

# Reset plotting parameters to default
par(mar = c(5, 4, 4, 2) + 0.1)
This script performs the following steps:
1. Loads proportion data for the subgroup analysis
2. Conducts a meta-analysis using raw proportions
3. Creates a forest plot and Baujat plot for the raw proportion analysis
4. Performs a logit transformation of the proportions
5. Conducts a meta-analysis using the logit-transformed proportions
6. Transforms the results back to the proportion scale


## Interpreting the Results

- The forest plots show the effect sizes and confidence intervals for each study and the overall effect.
- The Baujat plots help identify studies that contribute most to heterogeneity and overall results.
- Egger's test and the funnel plot assess potential publication bias.
- The trim-and-fill method provides an adjusted effect size estimate accounting for potential publication bias.
- For the subgroup analysis, both raw and logit-transformed results are provided. The logit transformation can help handle proportions close to 0 or 1.

## References for data extraction in meta-analysis

Andersson J, Luthra R, Hurtig P, Tideman M. Employer attitudes toward hiring persons with disabilities: A vignette study in Sweden. Journal of Vocational Rehabilitation 2015;43(1):41-50

Bricout JC, Bentley KJ. Disability status and perceptions of employability by employers. Social Work Research 2000;24(2):87-95.

Fyhn T, Sveinsdottir V, Reme SE, Sandal GM. A mixed methods study of employers’ and employees’ evaluations of job seekers with a mental illness, disability, or of a cultural minority. Work 2021;70(1):235-245

Manning C, White PD. Attitudes of employers to the mentally ill. Psychiatric Bulletin 1995;19(9):541-543.

Tsang HW, Corrigan PW, Fung KM, et al. Sino-American employer perspective about behavioral-driven health conditions: Predictive analyses. International Journal of Psychiatry in Clinical Practice 2012;16(4):284-292. 

Zissi A, Rontos C, Papageorgiou D, Pierrakou C, Chtouris S. Greek employers’ attitudes to employing people with disabilities: Effects of the type of disability. Scandinavian Journal of Disability Research 2007;9(1):14-25.

## Contributing

Please feel free to submit issues or pull requests if you have suggestions for improvements or find any bugs.

## License

MIT License

Copyright (c) 2024 Ninon Crestois, Ricardo Twumasi

## References

1. Viechtbauer, W. (2024). metafor: Meta-Analysis Package for R (Version 4.2-0) [Software Manual]. https://cran.r-project.org/web/packages/metafor/metafor.pdf
 
2. Wang, N. (2023). Conducting Meta-Analyses of Proportions in R. Journal of Behavioral Data Science, 3(2), 64-126. https://dx.doi.org/10.35566/jbds/v3n2/wang

3. Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48. https://dx.doi.org/10.18637/jss.v036.i03

4. Cochran, W. G. (1954). The combination of estimates from different experiments. Biometrics, 10(1), 101-129.

5. Chinn S. A simple method for converting an odds ratio to effect size for use in meta-analysis. Statistics in Medicine 2000;19(22):3127-3131.
   
6. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis: Wiley; 2011.

7. Duval, S., & Tweedie, R. (2000). Trim and fill: a simple funnel-plot–based method of testing and adjusting for publication bias in meta-analysis. Biometrics, 56(2), 455-463.

8. Egger, M., Smith, G. D., Schneider, M., & Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. Bmj, 315(7109), 629-634.

9. Higgins, J. P., & Thompson, S. G. (2002). Quantifying heterogeneity in a meta‐analysis. Statistics in medicine, 21(11), 1539-1558.

10. Higgins, J. P., Thompson, S. G., Deeks, J. J., & Altman, D. G. (2003). Measuring inconsistency in meta-analyses. Bmj, 327(7414), 557-560.
