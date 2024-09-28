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
