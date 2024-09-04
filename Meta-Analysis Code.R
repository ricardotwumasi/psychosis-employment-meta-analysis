# Load the metafor package
library(metafor)

# Input your data
data <- data.frame(
  author = c("Andersson et al.", "Bricout & Bentley", "Fyhn et al.", "Manning & White", "Tsang et al.", "Zissi et al."),
  year = c(2015, 2000, 2021, 1995, 2012, 2007),
  effect_size = c(0.533, 0.84, 0.473, 0.666, 1.99, 0.972),
  variance = c(0.0161, 0.0273, 0.00128, 0.0409, 0.0517, 0.0492),
  n = c(71, 86, 882, 109, 183, 102)
)
# Truncate study titles
data$TruncatedTitle <- c(
  "Andersson et al",
  "Bricout & Bentley",
  "Fyhn et al.",
  "Manning & White",
  "Tsang et al.",
  "Zissi et al."
)

# Run the meta-analysis
meta_analysis <- rma(yi = effect_size, vi = variance, data = data, method = "REML", test = "knha")

# View the results
summary(meta_analysis)

# Create forest plot with truncated titles
forest(meta_analysis, slab = data$TruncatedTitle, xlab = "Effect Size", 
       mlab = "Random-Effects Model", addpred = TRUE)

# Create a Baujat plot
baujat(meta_analysis)

# Egger's regression test
eggers_test <- regtest(meta_analysis, model="rma", predictor="sei")
print(eggers_test)

# Duval and Tweedie's trim-and-fill method
tf_meta_analysis <- trimfill(meta_analysis)
summary(tf_meta_analysis)

# Funnel plot with trim-and-fill
funnel(tf_meta_analysis)