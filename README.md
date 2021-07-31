# Model-Fitting-on-Lotka-Volterra
An modified all methods.
First time: run either combine_mahafy or combine_Hudson to see what this folder can generate. Comparison of all methods UKF, DRAM, PSO on Mahafy and Hudson_Bay dataset 

1. Chow Lorenz: comment the original code from Chow et al. in application of UKF on Lorenz system. They use MLE to get the time-invariant parameters 
2. Dual: modified REU student Dual UKF on prey-predator model. addressed the issues of intial parameters, covariance matrices and dataset (Mahafy). dual is not as efficient as joint UKF
3. Joint: modified REU student on prey-predator model. Added Mahafy and Hudson_bay dataset, used Mahafy's intial guesses, covariance matrices are simpler (scaled identity matrix). Created Joint_AD(dataset) function to run UKF for any particular dataset. All data are saved in UKFpar, UKF_data, UKF_time 
4. MCMC: only kept the folder with DRAM in it. Added mahafy and hudson_bay dataset, created DRAM_AD function to run DRAM for any particular dataset. All data are saved in csv file 
5. PSO: Added mahafy and hudson_bay dataset, created PSO_AD function to run DRAM for any particular dataset. All data are saved in csv file 
7. combine_Hudson: run all methods and compare their error for mahafy and Hudson_bay dataset 
8. Notes: Powerpoint presentations of what An learned from the codes and their resuts. 
