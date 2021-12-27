``Description``:        Script to estimate parameters of the Lotka-Volterra system.
               This script can be used to perform the Metropolis-Hastings 
               and DRAM MCMC methods. 
	     
 ``Direction``:   This script is meant to be run in full, ~15 min. User 
               specifications exist for DRAM vs Metropolis-Hastings and 
               naming the workspace. 

Execute function: Run_DRAM(dataset) where dataset = 'Mahaffy' or 'HudsonBay'. The datasets are in Excel format and contain year, hares and lynx population. 

E.g.: Run_DRAM('Mahaffy')

```Reference```
We use a script adapted from an algae example by M.J. Laine and rely on the 
functions from the mcmcstat library.


	lotkaVolterrafun.m
		Function to solve the ODE
	
	lotkaVolterrass.m
		Calculates a sum of squares for use in the acceptance criteria
	
	lotkaVolterrasys.m
		Function that describes the system of equations.
	
	LVmse.m
		Calculates mean squared error of the model prediction
	
	[error_est.m
		Calculates standard error of the estimate (another measure of algorithm 
		performance)] * not used

	
	fitInitialParams.m 
		Using fmincon finds minimum value of parameter values for which system is 
		viable. 
		These values are used as initial guesses for parameter values. Also 
		calculates an estimate for the variance of the prior function (model.sigma2).
	
	logprior.m
		Function for a user-specified prior function. Default in this implementation is a
		log-uniform prior.
		
	DRAM_pars: estimated parameters by DRAM 
	DRAM_time: reported runtime 
		

