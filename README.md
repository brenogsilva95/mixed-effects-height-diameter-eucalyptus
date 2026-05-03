# mixed-effects-height-diameter-eucalyptus
Mixed-effects modeling framework for height–diameter relationships in Eucalyptus urograndis using hierarchical data structures and advanced residual diagnostics.


# A Mixed-Effects Model Approach to Height–Diameter Relationships

## Introduction

Height–diameter models are widely used in forest inventories to estimate tree height from diameter at breast height (DBH), significantly reducing fieldwork costs and improving operational efficiency. However, modeling this relationship involves several statistical challenges, including nonlinearity, heteroscedasticity, non-normality, and hierarchical data structures.

In forest datasets, observations are naturally grouped across multiple levels such as forest regions, stands, and trees. Ignoring this structure can lead to biased estimates and reduced predictive accuracy.

This study proposes a mixed-effects modeling framework that extends Scolforo’s model by incorporating random effects and variance structures to better represent the complexity of forestry data :contentReference[oaicite:0]{index=0}.

---

## Materials and Methods

### Dataset

The data consist of observations from *Eucalyptus urograndis* plantations in Brazil, structured hierarchically across:

- Forest regions  
- Stands  
- Trees  

The dataset includes:

- Diameter at breast height ($D$)  
- Total height ($H$)  
- Dominant height ($H_d$)  
- Mean quadratic diameter ($D_g$)  
- Age ($I$)  

This hierarchical structure is fundamental for modeling variability across different spatial and biological levels :contentReference[oaicite:1]{index=1}.

---

### Modeling Approach

The initial fixed-effects model is based on Scolforo (1998), defined as:

$$
\log(H_{rsi}) = \beta_0 + \beta_1 \log(H_d) + \beta_2 \log\left(\frac{D_g}{D}\right) + \beta_3 \frac{1}{D \cdot I} + \beta_4 \frac{1}{D} + \varepsilon_{rsi}
$$

To account for hierarchical dependence, mixed-effects models were introduced:

$$
\log(H_{rsi}) = (\beta_0 + u_{2(rs)}) + \beta_1 \log(H_d) + \beta_2 \log\left(\frac{D_g}{D}\right) + \beta_3 \frac{1}{D \cdot I} + \beta_4 \frac{1}{D} + \varepsilon_{rsi}
$$

where:

- $u_{2(rs)} \sim N(0, \sigma^2_s)$ represents random effects at the stand level  
- $\varepsilon_{rsi} \sim N(0, \sigma^2)$ is the residual error  

---

### Estimation and Model Selection

Model selection followed a structured approach:

- Likelihood Ratio Tests (LRT)  
- Akaike Information Criterion (AIC)  
- Conditional AIC (cAIC)  
- Wald tests for fixed effects  

The estimation was performed using **REML (Restricted Maximum Likelihood)**.

The final selected model (**M4**) includes random effects at the stand level, providing the best fit according to AIC and likelihood criteria :contentReference[oaicite:2]{index=2}.

---

### Diagnostics and Residual Analysis

Model adequacy was assessed using:

- Least-confounded residuals  
- Randomized quantile residuals  
- Deviance residuals  
- Half-normal plots with simulated envelopes  

Local influence analysis was also applied to identify influential observations, revealing the impact of specific trees on parameter estimates :contentReference[oaicite:3]{index=3}.

---

### Simulation Study

A simulation study was conducted to evaluate residual performance under different scenarios involving:

- Number of forest regions  
- Number of stands  
- Number of trees  

Results showed that:

- Least-confounded residuals consistently outperformed other residual types  
- Randomized quantile residuals showed stable performance  
- Deviance residuals presented higher variability  

This reinforces the importance of appropriate residual diagnostics in mixed models :contentReference[oaicite:4]{index=4}.

---

## Results

The mixed-effects model with random effects at the stand level (M4) demonstrated:

- Superior model fit (lowest AIC and cAIC)  
- Improved residual behavior  
- Better representation of variability across hierarchical levels  

The removal of influential observations further improved model stability and parameter estimation.

The results highlight that incorporating hierarchical structure and variance heterogeneity significantly enhances predictive performance and model interpretability :contentReference[oaicite:5]{index=5}.

---

## Conclusion

This study demonstrates that mixed-effects models provide a robust and flexible framework for modeling height–diameter relationships in forestry.

The proposed approach:

- Accounts for hierarchical data structure  
- Improves prediction accuracy  
- Enhances interpretability  
- Provides robust diagnostic tools  

Additionally, the simulation study confirms the superiority of least-confounded residuals for model evaluation.

This work contributes a practical and reproducible framework for forest modeling, bridging theory and application in forest biometrics.

---

## Reference

- Akaike, H. (1974). A new look at the statistical model identification. IEEE transactions on automatic control, 19(6), 716-723. DOI: https://doi.org/10.1109/TAC.1974.1100705
- Atkinson, A. C. (1985). Plots, transformations, and regression an introduction to graphical methods of diagnostic regression analysis
- Baey, C., & Kuhn, E. (2023). varTestnlme: An R package for variance components testing in linear and nonlinear mixed-effects models. Journal of Statistical Software, 107, 1-32. DOI: https://doi.org/10.18637/jss.v107.i06
- Baey, C., Cournède, P. H., & Kuhn, E. (2019). Asymptotic distribution of likelihood ratio test statistics for variance components in nonlinear mixed effects models. Computational Statistics & Data Analysis, 135, 107-122. 
DOI: https://doi.org/10.1016/j.csda.2019.01.014
- Beckman, R. J., Nachtsheim, C. J., & Cook, R. D. (1987). Diagnostics for mixed–model analysis of variance. Technometrics, 29(4), 413-426. DOI: https://doi.org/10.1080/00401706.1987.10488269
- Butler, D. G., Cullis, B. R., Gilmour, A. R., & Gogel, B. J. (2009). ASReml-R reference manual. The State of Queensland, Department of primary industries and fisheries, Brisbane.
- Comets, E., Lavenu, A., & Lavielle, M. (2017). Parameter estimation in nonlinear mixed effect models using saemix, an R implementation of the SAEM algorithm. Journal of Statistical Software, 80, 1-41. DOI: https://doi.org/10.18637/jss.v080.i03
- Cook, R. D. (1986). Assessment of local influence. Journal of the Royal Statistical Society Series B: Statistical Methodology, 48(2), 133-155. DOI: https://doi.org/10.1111/j.2517-6161.1986.tb01398.x
- Cordeiro, G. M., Demétrio, C. G. B., & Moral, R. D. A. (2024). Modelos lineares generalizados e aplicações. Editora Blucher.
- Verbeke, G., & Lesaffre, E. (1996). A linear mixed-effects model with heterogeneity in the random-effects population. Journal of the American Statistical Association, 91(433), 217-221. DOI: https://doi.org/10.1080/01621459.1996.10476679
- Verbeke, G., & Molenberghs, G. (2000). Estimation of the marginal model. Linear mixed models for longitudinal data, 41-54
- Verbeke, G., Molenberghs, G., & Verbeke, G. (2009). Linear mixed models for longitudinal data. Springer New York
- West, B. T., Welch, K. B., & Galecki, A. T. (2022). Linear mixed models: a practical guide using statistical software. Chapman and Hall/CRC
- Webster, R., Welham, S. J., Potts, J. M., & Oliver, M. A. (2006). Estimating the spatial scales of regionalized variables by nested sampling, hierarchical analysis of variance and residual maximum likelihood. Computers & Geosciences, 32(9), 1320-1333. DOI: https://doi.org/10.1016/j.cageo.2005.12.002
- Waternaux, C., Laird, N. M., & Ware, J. H. (1989). Methods for analysis of longitudinal data: blood-lead concentrations and cognitive development. Journal of the American Statistical Association, 84(405), 33-41. DOI: https://doi.org/10.2307/2289844
- Wolfinger, R. (1993). Covariance structure selection in general mixed models. Communications in statistics-Simulation and computation, 22(4), 1079-1106. DOI: https://doi.org/10.1080/03610919308813143
