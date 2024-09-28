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
