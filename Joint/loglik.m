% AUTHOR: An Do Dela 
% Negative Log likelihood fuction for UKF estimation for Lotera model
% INPUT: 
%       initPar initial parameter guess 
%       y: observable data
%       x0: initial observables data
%       Px: convariance matrix for state variables
%       Pparam: covariance matrix for parameters
%       InfDS: information structure for state

%       pNoise: process noise for state
%       oNoise: measurement noise for state
%OUTPUT: 
% negative log likelihood function

function negloglik= loglik(x0,y,Px,InfDS,pNoise,oNoise)
oNoise.cov = diag(x0(1:2)); % why update the measurement noise cov matrix this way?  
initPar=x0(4:6); 
InfDS.par = initPar;
    % Full Joint Unscented Kalman Filter step
  % =================================
  [xhat, Px, xh_, Px_, yh_, inov, Py, KG, negloglik,Pxall, Pyall, Px_all]=...
      UKF(x0,Px,pNoise, oNoise, y,[],[],InfDS); 
  
  
