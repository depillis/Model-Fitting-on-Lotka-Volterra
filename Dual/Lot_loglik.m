% AUTHOR: An Do Dela 
% Negative Log likelihood fuction for UKF estimation for Lotera model
% INPUT: 
%       initParL initial parameter guess 
%       y: observable data
%       x0: initial observables data
%       Px: convariance matrix for state variables
%       Pparam: covariance matrix for parameters
%       InfDS: information structure for state
%       InfDSparams: information structure for params

%       pNoise: process noise for state
%       pNoise_params: : process noise for parameters
%       oNoise: measurement noise for state
%       oNoise_params: measurement noise for parameters
%OUTPUT: 
% negative log likelihood function


function negloglik = Lot_loglik(initPar,y,x0,Px,Pparams,InfDS,InfDSparams, pNoise,pNoise_params, oNoise, oNoise_params)
initPar
oNoise_params.cov = diag([initPar]);
InfDS.par = initPar;

T = length(y); 
    
    for t = 2:T
        [param_next, P_param_next] =...
           UKF_Params(initPar, Pparams, pNoise_params, oNoise_params,...
           y(:,t), InfDSparams, x0, InfDSparams.Xdim); %estimate parameters
        
       %set new values based on output
        params0 = param_next;
        P0_params_real = P_param_next;
        InfDS_X_real.par = params0;
        InfDS_Param_real.par = params0;
        [x_next, Px_next] = UKF(x0_real,P0_x_real,Q_x_real, R_x_real,...
           y(:,t),InfDS_X_real); %estimate states
        
       %set new values based on output
        x0_real = x_next;
        P0_x_real = Px_next;
        xhat_real(:, t) = x_next;
        paramshat_real(:,t) = param_next;
    end
    
    
end 





    % Full Unscented Kalman Filter step
  % =================================
  [xhat, Px, xh_, Px_, yh_, inov, Py, KG, negloglik,Pxall, Pyall, Px_all]=...
      ukfaddmultiP(x0,Px,pNoise, oNoise, y,[],[],InfDS); 
end

