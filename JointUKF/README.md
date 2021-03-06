``Description``: This folder is meant for using the Joint UKF to perform estimation of states and parameters on the Lotka-Volterra system. We adapted and modified Chow et al. into the Lotka Volterra model, which has not been done before.

1. Run_JointUKF: the main file to run the estimation
2. Lotka_Volterra_Model - contains set up of ODEs for the model that is used in the ODE solver throughout
3. MeasurementFcn - set up of the measurement function for the UKF
4. UKF - performs ALL iterations of the Joint UKF algorithm
5. Joint_UKF_Real_Data_Figures - produces figures found in the paper (as well as some extras)
6. Figures folder contains pre-generated figures reported in the manuscript.

``Reference`` Chow, Sy-Miin, Emilio Ferrer, and John R. Nesselroade. "An unscented Kalman filter approach to the estimation of nonlinear dynamical systems models." Multivariate Behavioral Research 42.2 (2007): 283-321.


