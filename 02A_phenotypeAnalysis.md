# Phenotype Analysis Workflow 

## Table of Contents

1. [Data](#data)
2. [Analysis 1 - Extra-pallial fluid all timepoints](#one)
3. [Analysis 2 - Calcification Response](#two)

---

## Data <a name="data"></a>

[**Link to data**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/)  

## Analysis 1 - Extra-pallial fluid all timepoints <a name="one"></a>

Analysis looking at extra-pallial fluid pH across time (24H,48H,9D,22D,50D,80D) considering time and treatment as fixed factors and several random effects. Include full and final models and final figure.

**Description**  
A linear mixed model was used to determine the effect of the explanatory variables of treatment and time on the response variable of pHEPF. The full model included the explanatory variables of treatment and time and their interaction as categorical fixed effects, and tank (nested within shelf) and oyster collection site as random effects. The model was performed in R using the lme4 package (Bates et al., 2015). A step-down strategy with likelihood ratio tests (based on the degrees of freedom estimated using the Satterthwaite method) implemented in lmerTest, following the principle of marginality (Kuznetsova et al., 2017), was used to select the most parsimonious linear mixed effects model. This model selection approach evaluates two different interpretations of pHEPF: (1) the measured value of pHEPF and (2) pHEPF relative to treatment seawater pH (ΔpH = pHEPF - pHseawater), which may be a better indicator of active pHEPF regulation. To evaluate significant differences between pHEPF and ΔpH under the moderate and high OA treatments for each time point, a series of post hoc comparisons with a Tukey correction were performed using the multcomp package (Hothorn et al., 2008). In addition, the change in ΔpH over time was examined by performing a series of single sample t-tests to determine whether each treatment at each time point was significantly different from seawater pH. We corrected for testing multiple hypotheses using a Benjamini-Hochberg correction (Benjamini and Hochberg, 1995).

* [Script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/AE17_fig1_EPFtimeseries.R)

![](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/results/manuscript/figures/Figure1/figure1.png)

## Analysis 2 - Long term EPF pH and Calcification Response <a name="two"></a>

Analysis examining long term EPF pH and calcification rate in response to treatment or EPF pH using individuals collect at the final two timepoints of the exposure (50D and 80D) and for bouyant weights collected from those individuals at day 33-34 of the exposure.

**Description**  
Linear mixed effect models were used to examine the effect of pCO2 on long term EPF response, both measured pHEPF  and ΔpH. The full models included EPF response (pHEPF or ΔpH ) as the response variable, treatment and time and their interaction as categorical fixed effects and tank (nested within shelf) and oyster collection site as random effects. We included samples from our two long term time points (day 50 and 80) in the model (n = 35). Model selection and post hoc comparisons were performed using the approach described above.  
  
Another set of linear mixed effect models were used to examine the effect of either treatment or pHEPF on calcification rates. In the full models, tank (nested in shelf) and collection site were included as random effects. Calcification rate was based on buoyant weights measured at a single time (day 33) for all oysters (n = 35), so time could not be included as a fixed effect. Treatment was handled as a continuous variable and was calculated as the average pCO2 for each tank across the duration of the exposure prior to sampling. Model selection was performed using the approach described above. Regression was used to examine the effect of either treatment or pHEPF and calcification rate.  

* [Script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/AE17_fig2_EPFfinal_calcification.R)

![](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/results/manuscript/figures/Figure2/Figure2.png)

