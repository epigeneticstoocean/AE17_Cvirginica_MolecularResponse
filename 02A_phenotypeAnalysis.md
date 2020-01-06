# Phenotype Analysis Workflow 

## Table of Contents

1. [Data](#data)
2. [Analysis 1 - Extra-pallial fluid all timepoints](#one)
3. [Analysis 2 - Calcification Response](#two)

---

## Data <a name="data"></a>

* [**Link to scripts**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/src/)  
* [**Link to data**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/)  

## Analysis 1 - Extra-pallial fluid all timepoints <a name="one"></a>

Analysis looking at extra-pallial fluid pH across time (24H,48H,9D,22D,50D,80D) considering time and treatment as fixed factors and several random effects. Include full and final models and final figure.

**Description**  
A linear mixed model was used to determine the effect of the explanatory variables of treatment and time on the response variable of pHEPF. The full model included the explanatory variables of treatment and time and their interaction as categorical fixed effects, and tank (nested within shelf) and site as random effects. The model was performed in R using the lme4 package (v.1.1-21; Bates et al., 2015). A model selection approach was used to determine the most parsimonious linear mixed model. We used a step-down strategy with likelihood ratio tests (based on the degrees of freedom estimated using Satterthwaite method) implemented in the lmerTest, following the principle of marginality (v.3.1-0; Kuznetsova et al., 2017). This model selection approach was used to evaluate two different interpretations of EPF pH: (1) the measured value (pHEPF) and (2) as the EPF pH relative to external seawater (pHEPF, relative = pHEPF - pHseawater), which may be a better indicator of active pHEPF regulation. To evaluate significant differences between ambient pHEPF and OA pHEPF (OA 1000 or OA 2800) for each timepoint a series of post hoc comparisons with a Tukey correction were performed using the multcomp package (v.1.4-10; Hothorn et al., 2008). Next we examined the change in pHEPF, relative over time relative to the environment by performing a series of one sided t-tests (mu = 0) to determine whether each treatment at each timepoint was significantly different than the seawater pH. We corrected for testing multiple hypotheses using a Benjamini-Hochberg correction.  

* [Script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/Phenotype/EPFpH_timeseries.R)

![](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/results/figures/Figure1.jpg)

## Analysis 2 - Calcification Response <a name="two"></a>

Analysis examining calcification rate in response to treatment or EPF pH using individuals collect at the final two timepoints of the exposure (50D and 80D) and for bouyant weights collected from those individuals at day 33-34 of the exposure.

**Description**  
Calcification rate was determined by estimating calcification over an approximately 33 day window using the buoyant weight method described above. Measurements were performed on all oysters remaining in the experiment which included oysters sampled on day 50 and 80 of the exposure (n = 35, 5-6 samples per treatment per time point). To examine the effect of OA treatment on calcification rate, a linear mixed effect model was used including treatment as a fixed effect and tank nested in shelf and site as a random effects. Importantly, calcification was based on bouyant weights collected at a single time (day 33) for all oysters, so time was not included as a fixed effect. Here, treatment was handled as a continuous variable and was calculated as the average pCO2 for each tank across the duration of the exposure prior to sampling. The same model testing approach described above was used to determine the best model. Regression was used to examine the relationship between treatment and calcification. 

Linear regression was also used to examine the relationship between pHEPF and calcification for oyster samples where both pieces of information overlapped (n = 35). The linear mixed model for this regression included pHEPF as the primary explanatory variable and tank nested in shelf and site as a random effects. Model selection was performed using the approach described above to identify the best model.

* [Script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/Phenotype/Calcification.R)

![](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/results/figures/Figure2.pdf)

