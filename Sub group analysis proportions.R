# Load necessary library
library(metafor)

# Create data frame with study data
data <- data.frame(
  author = c("Zissi et al.", "Manning et al.", "Tsang"),
  k = c(74, 72, 161),
  n = c(102, 109, 183)
)

# Calculate proportions
data$p <- data$k / data$n

# Calculate variance for raw proportions
data$variance <- (data$p * (1 - data$p)) / data$n

# Truncate study titles
data$TruncatedTitle <- c(
  "Zissi et al.",
  "Manning & White",
  "Tsang et al.")
  
# Meta-analysis using raw proportions
rma_raw <- rma(yi = p, vi = variance, data = data)

# Summary of the meta-analysis
summary(rma_raw)
forest(rma_raw, slab = data$TruncatedTitle, xlab = "Effect Size", 
       mlab = "Random-Effects Model", addpred = TRUE)
baujat(rma_raw)

# Add logit-transformed proportions
data$logit_p <- log(data$p / (1 - data$p))

# Calculate variance of the logit-transformed proportions
data$logit_variance <- (data$variance) / (data$p * (1 - data$p))^2

# Meta-analysis using logit-transformed proportions
rma_logit <- rma(yi = logit_p, vi = logit_variance, data = data)

# Summary of the meta-analysis
summary(rma_logit)
# Calculate the transformed effect size (back to proportions)
transformed_proportion <- function(logit) {
  exp(logit) / (1 + exp(logit))
}

# Get transformed results
logit_summary <- summary(rma_logit)
transformed_mean <- transformed_proportion(logit_summary$yi)
transformed_mean
