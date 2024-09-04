# Meta-Analysis of Employer Attitudes Towards Individuals with Psychosis

This repository contains R code for replicating our meta-analysis on employer discrimination towards individuals with psychosis (Crestois et al., in press). The analysis includes an overall meta-analysis, a subgroup analysis using proportion data, and effect size calculations for individual studies.

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

```R
#Calculation for Anderson et al. (2015)
mean_psychosis <- 8.23  # Mean of the psychosis group
SD_psychosis <- 4.26    # Standard deviation of the psychosis group
reference_point <- 10.5 # Midpoint reference, midpoint assumed to represent neutral interest
N <- 71
SMD <- (reference_point - mean_psychosis) / SD_psychosis
Var_SMD <- (1 / N) + (SMD^2 / (2 * (N - 1)))

#Calculating effect size of Bricout et al. (2000)
mean1 <-139.59 #extracted from table 2
mean2 <- 126.15 
sd1 <- 17.04
sd2 <-15.04
n1 <-74
n2 <-86
Sp <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
d <- (mean1 - mean2) / Sp
Var_d <- (n1 + n2) / (n1 * n2) + (d^2 / (2 * (n1 + n2)))

#Calculating effect size for Fyhn et al. (2015)
a <- 692   #negative/neutral assessment in schizophrenic symptom group
b <- 525  #negative/neutral assessment in control group
n1 <- 882  # sample size schizophrenic group
n2 <- 1074  # sample size control group
p_treat <- a/n1
p_contr <- b/n2
RR <- p_treat/p_contr
log_RR <- log(RR)
Var_log_RR <- (1 / a) - (1 / n1) + (1 / b) - (1 / n2)

#Calculating effect size of Manning (1995) 
k<-72 #from table 1, 66% of 109 people 
n<-109 #from table 1, n=109
p<-k/n
SE<-sqrt((p*(1-p))/n)
logit_p <- log(p / (1 - p))
SE_logit_p <- SE / (p * (1 - p))
Var_logit_p <- SE_logit_p^2

#Calculating effect size of Tsang et al.(2012)
k<-161  # from table 1 mental illness, addition of all not offer
n<-183 #from table 1 mental illness, addition of all conditions
p<-k/n
SE<-sqrt((p*(1-p))/n)
logit_p <- log(p / (1 - p))
SE_logit_p <- SE / (p * (1 - p))
Var_logit_p <- SE_logit_p^2

#Calculating effect size of Zissi et al.(2007)
k<- 74 #from attitudes to employing people with disability, 27% willing to employ people with a hospitalisation record for schizophrenia
n<-102 #from sample section
p<-k/n
SE<-sqrt((p*(1-p))/n)
logit_p <- log(p / (1 - p))
SE_logit_p <- SE / (p * (1 - p))
Var_logit_p <- SE_logit_p^2
```

This script calculates effect sizes and variances for individual studies using different methods depending on the data available:
1. Standardized Mean Difference (SMD) for Anderson et al. (2015)
2. Cohen's d for Bricout et al. (2000)
3. Risk Ratio (RR) for Fyhn et al. (2015)
4. Logit-transformed proportions for Manning (1995), Tsang et al. (2012), and Zissi et al. (2007)

### 2. Meta-Analysis (meta_analysis.R)

```R
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
```

This script performs the following steps:
1. Loads the data for the meta-analysis
2. Conducts a random-effects meta-analysis using REML
3. Creates a forest plot of the results
4. Generates a Baujat plot to assess study influence
5. Performs Egger's test for publication bias
6. Applies the trim-and-fill method to adjust for potential publication bias
7. Creates a funnel plot with the trim-and-fill results

### 3. Subgroup Analysis (subgroup_analysis.R)

```R
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
data$logit_variance <- (data$variance) / (data$p * (1 - data$p))^2 ####Ninon please check this ####

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
```

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

## Contributing

Please feel free to submit issues or pull requests if you have suggestions for improvements or find any bugs.

## License

MIT License

Copyright (c) 2024 Ninon Crestois, Ricardo Twumasi

## References

1. Wang, N. (2023). Conducting Meta-Analyses of Proportions in R. Journal of Behavioral Data Science, 3(2), 64-126. https://dx.doi.org/10.35566/jbds/v3n2/wang

2. Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48. https://dx.doi.org/10.18637/jss.v036.i03

3. Cochran, W. G. (1954). The combination of estimates from different experiments. Biometrics, 10(1), 101-129.

4. Duval, S., & Tweedie, R. (2000). Trim and fill: a simple funnel-plot–based method of testing and adjusting for publication bias in meta-analysis. Biometrics, 56(2), 455-463.

5. Egger, M., Smith, G. D., Schneider, M., & Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. Bmj, 315(7109), 629-634.

6. Higgins, J. P., & Thompson, S. G. (2002). Quantifying heterogeneity in a meta‐analysis. Statistics in medicine, 21(11), 1539-1558.

7. Higgins, J. P., Thompson, S. G., Deeks, J. J., & Altman, D. G. (2003). Measuring inconsistency in meta-analyses. Bmj, 327(7414), 557-560.
