```DESCRIPTION```
The goal is to fit 4 parameters to a linear model: a, the intrisic growth rate of the prey population, b, the rate of predation, c, the mortality rate of the predator population in the absence of prey, and d, the growth rate of the predator population per prey, to produce a posterior density function for each of the parameters.
Prior distributions are assumed to be uniform, and we use a Gaussian likelihood function. 

```Reference```
The code here is adapted from M.J.Laine code found in the folder MJLaine_Algae_example.
The mhsample function is implemented in Reuel_Smith_mhsample.m and based on code from: https://crr.umd.edu/bayesian-parameter-estimation-single-data-set-example-problem-52-matlab


