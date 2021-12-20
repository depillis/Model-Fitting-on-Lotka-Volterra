# Model-Fitting-on-Lotka-Volterra

We implement all parametrization techniques (Joint UKF, PSO, DRAM) on 2 datasets: [Mahaffy](https://jmahaffy.sdsu.edu/courses/f09/math636/lectures/lotka/qualde2.html) and [Hudson Bay](https://gist.github.com/michaelosthege/27315631c1aedbe55f5affbccabef1ca)

**Instruction**: Run either ```Run_MahaffyData``` or ```Run_HudsonBayData``` to generate all the figures in the manuscript. All the figures inlcude: the predicted prey and predator population by each method, error magntidue, time-variant parameters, parameter distributions, and correlation on both datasets.

Repository overview: 
1. ```Chow Lorenz```: the folder contains the original code from [[Chow et al.]](#1). in application of UKF on Lorenz system. They combine MLE and UKF to obtain the time-invariant parameters.
2. ```Dual```: The folder contains the modified REU student Dual UKF on prey-predator model adapted from  [[Chow et al.]](#1). It turned out that Dual UKF is not as efficient as Joint UKF.
3. ```Joint```: The folder contains the modified REU student on prey-predator model. The main function is ```Run_JointUKF(dataset)``` to run UKF for any particular dataset. 
4. ```MCMC```: The folder contains the DRAM implementation. The main function to run is ``Run_DRAM(dataset)``, which runs DRAM for any particular dataset. 
5. ```PSO```: The folder contains the PSO implementation. The main function to run is ``Run_PSO(dataset)`` to run DRAM for any particular dataset. 
7. ```Run_HudsonBayData```: the main function to run all methods and their analysis on Hudson_bay dataset 
8. ```Run_MahaffyData```: the main function to run all methods and their analysis on Mahaffy dataset 
9. ```Notes```: meeting presentations and results from what An learned from previous work. 
10. ``Figures``: The folder contains all the figures (fig files and eps files) that appear in the manuscript. 



## References
<a id="1">[1]</a> 
Chow, Sy-Miin, Emilio Ferrer, and John R. Nesselroade. "An unscented Kalman filter approach to the estimation of nonlinear dynamical systems models." Multivariate Behavioral Research 42.2 (2007): 283-321.
