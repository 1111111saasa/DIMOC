DIMOC
================

This study introduces a novel framework for detecting confounders in high-dimensional mediation models, addressing challenges in mediation analysis caused by confounding variables. The data-driven method uses symmetric data aggregation to construct test statistic pairs and applies a double filter for confounder identification. This approach controls the false discovery rate (FDR) while enhancing confounder detection. Theoretical and empirical analyses demonstrate that the method controls FDR asymptotically and ensures high detection power. Simulations show robust performance, emphasizing the importance of confounder detection in improving the accuracy of mediation effect estimates.

OVERVIEW
================
This file contains the exact code necessary to reproduce the numerical results presented in the paper: Large-Scale Mediator-Outcome Confounder Detection for Correcting Mediation Effects in High-Dimensional Mediation Models.

The code is organized into three main folders:

1.DIMOC: This folder generates the toy example in Section 2 and the simulation results in Section 4.1.

2.Mediation Correction: This folder produces the simulation results discussed in Section 4.2.

3.Extension on GLM Model: This folder contains the results for the GLM model presented in the Supplementary Material.

4.Additional Exprements：This folder contains scripts for some additional experimental results.

Guides for the codes in folder \code\DIMOC
=========================================
## Enviroment 
R-4.3.2

## Codes for reproducing results in Section 2
\code\DIMOC\total_simu_code_DIMOC.R

Please note that the code contains several parameters, which we will explain one by one:

  `n:sample size`
  
  `B:refinement times`
  
  `p:dimension`
  
  `ratio:splitting ratio`
  
  `pi10:portion of alpha!=0,beta=0 in the mediaton model`
  
  `pi01:portion of alpha=0,beta!=0 in the mediaton model`
  
  `nonzero:Locations of mediators`
  
  `q:targrt FDR level`
  
  `tau:signal strength`
  
  `corr:dependence`
  
  To reproduce the results of the toy example, set 𝑛=400 and 𝑝=4000.



## Codes for reproducing results in Section 4.1

\code\DIMOC\total_simu_code_DIMOC.R

\code\DIMOC\oracle_DIMOC.R
To reproduce the results in figure 3, we need to fix n=400 and alter the p from 1000 to 5000 in the interval of 500.

For the results in Table 1, we need to alter the `tau` and `pi01`.   

Guides for the codes in folder \code\Mediation effect correction
===============================================

## Codes for reproducing results in Section 4.2

\code\Mediation effect correction\mediation effect correction.R

To reproduce the results in section 4.2, we need to alter the `corr` and `iota`.
Note that this part applys the method of DACT [Liu, Z., Shen, J., Barfield, R., Schwartz, J., Baccarelli, A. A., & Lin, X. (2022). *Large-scale hypothesis testing for causal mediation effects with applications in genome-wide epigenetic studies*. *Journal of the American Statistical Association, 117*(537), 67–81.](https://doi.org/10.1080/01621459.2021.1944104).

Additionally, this part utilizes the method HIMA2. [Perera, C., Zhang, H., Zheng, Y., Hou, L., Qu, A., Zheng, C., Xie, K., & Liu, L. (2022). *Hima2: High-dimensional mediation analysis and its application in epigenome-wide DNA methylation data*. *BMC Bioinformatics, 23*(1), 296.](https://doi.org/10.1186/s12859-022-04884-w)



Guides for the codes in folder \code\Extension on GLM
===============================================

## Codes for reproducing results in SUPPLEMENTARY MATERIAL

\code\Extension on GLM\glm_total_BB.R

\code\Extension on GLM\glm_total_BC.R

\code\Extension on GLM\oracle_glm_total_BB.R

\code\Extension on GLM\oracle_glm_total_BC.R



