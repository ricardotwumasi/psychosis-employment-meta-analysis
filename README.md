# Meta-Analysis of Employer Attitudes Towards Individuals with Psychosis

This repository contains R code for replicating our meta-analysis on employer discrimination towards individuals with psychosis (Crestois et al., in preparation). The analysis includes an overall meta-analysis, a subgroup analysis using proportion data, and effect size calculations for individual studies.

This code may be repurposed for other meta-analyses or updated to add studies as they are published.

## Contributors

- Ninon Crestois
- Ricardo Twumasi

## Requirements 

- R (version 4.0.0 or higher recommended)
- R packages: `metafor`, `dplyr`, `tidyr`, `purrr`

## Installation

1. Clone this repository: git clone https://github.com/ricardotwumasi/psychosis-employment-meta-analysis.git

2. Install the required R packages:

install.packages(c("metafor", "dplyr", "tidyr", "purrr"))

Usage
The analysis is split into three main scripts:

effect_size_calculations.R: Calculates effect sizes for individual studies
meta_analysis.R: Conducts the overall meta-analysis
subgroup_analysis.R: Performs a subgroup analysis using proportion data

Run each script separately in R.
Code Explanation
1. Effect Size Calculations (effect_size_calculations.R)
This script:

Loads necessary packages
Defines functions for converting various effect size measures to Hedges' g
Creates a structured data frame with study data
Validates study data for completeness and sufficient sample sizes
Calculates effect sizes (Hedges' g) for each study
Saves the calculated effect sizes to a CSV file

2. Meta-Analysis (meta_analysis.R)
This script:

Loads the effect size data calculated in the previous step
Conducts a random-effects meta-analysis using the rma function from metafor
Creates a forest plot to visualize study effect sizes and the overall effect
Generates a Baujat plot to identify influential studies
Produces a funnel plot with trim-and-fill method to assess publication bias
Calculates Rosenthal's and Orwin's Fail-Safe N to evaluate the robustness of results

3. Subgroup Analysis (subgroup_analysis.R)
This script:

Loads proportion data for the subgroup analysis
Conducts a meta-analysis using raw proportions
Creates a forest plot and Baujat plot for the raw proportion analysis
Performs a logit transformation of the proportions
Conducts a meta-analysis using the logit-transformed proportions
Transforms the results back to the proportion scale for interpretation

Interpreting the Results

Forest plots: Show effect sizes and confidence intervals for each study and the overall effect
Baujat plots: Help identify studies that contribute most to heterogeneity and overall results
Egger's test and funnel plot: Assess potential publication bias
Trim-and-fill method: Provides an adjusted effect size estimate accounting for potential publication bias
Subgroup analysis: Both raw and logit-transformed results are provided to handle proportions close to 0 or 1

## References for data extraction in meta-analysis

Andersson J, Luthra R, Hurtig P, Tideman M. Employer attitudes toward hiring persons with disabilities: A vignette study in Sweden. Journal of Vocational Rehabilitation 2015;43(1):41-50

Bricout JC, Bentley KJ. Disability status and perceptions of employability by employers. Social Work Research 2000;24(2):87-95.

Fyhn T, Sveinsdottir V, Reme SE, Sandal GM. A mixed methods study of employers’ and employees’ evaluations of job seekers with a mental illness, disability, or of a cultural minority. Work 2021;70(1):235-245

Manning C, White PD. Attitudes of employers to the mentally ill. Psychiatric Bulletin 1995;19(9):541-543.

Tsang HW, Corrigan PW, Fung KM, et al. Sino-American employer perspective about behavioral-driven health conditions: Predictive analyses. International Journal of Psychiatry in Clinical Practice 2012;16(4):284-292. 

Zissi A, Rontos C, Papageorgiou D, Pierrakou C, Chtouris S. Greek employers’ attitudes to employing people with disabilities: Effects of the type of disability. Scandinavian Journal of Disability Research 2007;9(1):14-25.

## Contributing

Please feel free to submit issues or pull requests if you have suggestions for improvements or find any bugs, or would like to continue this meta-analysis with new studies.

## AI Statement

This code was edited with the assistance of Claude Sonnet 3.5 (Anthropic, San Francisco: CA)

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
