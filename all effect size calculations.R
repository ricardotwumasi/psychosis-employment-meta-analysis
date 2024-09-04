#Calculation for Anderson et al. (2015)
#Data extracted from table 1
mean_psychosis <- 8.23  # Mean of the psychosis group
SD_psychosis <- 4.26    # Standard deviation of the psychosis group
reference_point <- 10.5 # Midpoint reference, midpoint assumed to represent neutral interest
N <- 71
# Calculate the Standardized Mean Difference (SMD)
SMD <- (reference_point - mean_psychosis) / SD_psychosis
SMD
# Calculate the variance of the SMD
Var_SMD <- (1 / N) + (SMD^2 / (2 * (N - 1)))
Var_SMD


#Calculating effect size of Bricout et al. (2000)
# Calculate the standardised mean difference
mean1 <-139.59 #extracted from table 2
mean2 <- 126.15 
sd1 <- 17.04
sd2 <-15.04
n1 <-74
n2 <-86
# Calculate pooled standard deviation (Sp)
Sp <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
# Calculate Cohen's d
d <- (mean1 - mean2) / Sp
d
# Calculate the variance of Cohen's d
Var_d <- (n1 + n2) / (n1 * n2) + (d^2 / (2 * (n1 + n2)))
Var_d

#Calculating effect size for Fyhn et al. (2015)
# Define data #from table 2
a <- 692   #negative/neutral assessment in schizophrenic symptom group
b <- 525  #negative/neutral assessment in control group
n1 <- 882  # sample size schizophrenic group
n2 <- 1074  # sample size control group

# Calculate the risks
p_treat <- a/n1
p_contr <- b/n2

# Calculate the risk ratio
RR <- p_treat/p_contr
RR
# Calculate the log-risk ratio and its standard error
log_RR <- log(RR)
log_RR
# Calculate the variance of the log risk ratio
Var_log_RR <- (1 / a) - (1 / n1) + (1 / b) - (1 / n2)
Var_log_RR

#Calculating effect size of Manning (1995) Attitudes of employers to the mentally ill
k<-72 #from table 1, 66% of 109 people 
n<-109 #from table 1, n=109
p<-k/n
#Calculating standard error
SE<-sqrt((p*(1-p))/n)
# Logit transformation of the proportion
logit_p <- log(p / (1 - p))
logit_p
# Standard Error of the logit-transformed proportion using delta method
SE_logit_p <- SE / (p * (1 - p))
# Variance of the logit-transformed proportion
Var_logit_p <- SE_logit_p^2
Var_logit_p


#Calculating effect size of Tsang et al.(2012)
k<-161  # from table 1 mental illness, addition of all not offer
n<-183 #from table 1 mental illness, addition of all conditions
p<-k/n
#Calculating standard error
SE<-sqrt((p*(1-p))/n)
# Logit transformation of the proportion
logit_p <- log(p / (1 - p))
logit_p
# Standard Error of the logit-transformed proportion using delta method
SE_logit_p <- SE / (p * (1 - p))
# Variance of the logit-transformed proportion
Var_logit_p <- SE_logit_p^2
Var_logit_p


#Calculating effect size of Zissi et al.(2007)
k<- 74 #from attitudes to employing people with disability, 27% willing to employ people with a hospitalisation record for schizophrenia
n<-102 #from sample section
p<-k/n
#Calculating standard error
SE<-sqrt((p*(1-p))/n)
# Logit transformation of the proportion
logit_p <- log(p / (1 - p))
logit_p
# Standard Error of the logit-transformed proportion using delta method
SE_logit_p <- SE / (p * (1 - p))
# Variance of the logit-transformed proportion
Var_logit_p <- SE_logit_p^2
Var_logit_p


