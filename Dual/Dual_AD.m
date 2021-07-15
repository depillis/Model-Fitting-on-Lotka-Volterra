% AUTHOR: An Do Dela 
% Date: June, 3, 2021
% Purpose: Fitting the Prey-predator model to noisy time series using UKF 
% Adapted from Sy-Miin Chow (schow@nd.edu)
% The appropriate references and descriptions of the algorithm can be found in
% Chow, Ferrer & Nesselroade (2005). 
% An unscented Kalman filter approach to the estimation 
% of nonlinear dynamical systems models.


clear all;
clc;
echo off;
ftsz= 20; % fontsize
nloop = 1;

Ny = 2;  %Number of observables
Nx = 2;  %Number of states
Nstates=2; %number states
Npar=4;   %number parameters
Nsubj=1;

%===== load data =====% 
HLData = load('HaresLynxData_AD.mat');  %Load dataset
rawData = HLData.Mahafy;
time_vector = rawData(:,1);

T = length(rawData);  %Number of data points
x(1:2,:) = rawData(:, 2:3)'; %Get just the predator and prey populations
y = zeros(Ny, T); %will hold observables

rand('state',sum(100*clock));   
randn('state',sum(100*clock));   

%===== True parameters from Mahafy =======%
% alpha, gamma, beta, delta 
% a1, b1, a2, b2

truepar = [0.4; 0.8; 0.018; 0.023];

%========SET UKF PARAMETERS==============%
alpha = 10e-3; % UKF : point scaling parameter
beta = 2; %Prior knowledge of distribution (for gaussian use 2)
kappa = 0;    %second scaling parameter, set to 0 or 3 - L 


% GENERATE THE DATA
% no need to 
% This should be the data extraction function
for t = 1:T
    y(:,t) = feval('MeasurementFcn', x(:,t), [], [], [], []); %Call measurement function
end

%============ PLOT THE DATA ========================
figure(1)
 plot(time_vector, x(1,:), '-*', 'Linewidth', 2); hold on;
ylabel('Prey Populations');
xlabel('Year');
set(gca, 'fontsize', ftsz)

figure(2)
 plot(time_vector, x(2,:), '-*', 'Linewidth', 2); hold on;
ylabel('Predator Populations');
xlabel('Year');
set(gca, 'fontsize', ftsz)

%===================================================

fprintf('Now estimating the model parameters...')
fprintf('\n')

x0_real = y(:,1);
params0 = [0.4; 0.8; 0.018; 0.023]; % alpha, gamma, beta, delta 
                                    % a1, b1, a2, b2
P0_x_real= diag([repmat(1,Nx,1)]); 
P0_params_real = diag([repmat(1,Npar,1)]);
Q_x_real.cov = diag(diag([repmat(10e-4,Nx)])); %[]; %initialize structs for process noise
Q_param_real.cov =  diag(diag([repmat(10e-4,Npar)])); %[];
R_x_real = []; %Initialize structs for measurement noise
R_param_real = [];
R_param_real.cov = diag(params0)


%==ALLOCATION  state and parameter vectors========%
xhat_real = zeros(Nx,T);
paramshat_real = zeros(Npar,T);
xhat_real(:,1) = x(:,1);
paramshat_real(:,1) = params0';

%SET UP STATE DATA STRUCTURE AND INFO FOR STATES
InfDS_X_real = [];
InfDS_X_real.spkfParams = [alpha beta kappa];
InfDS_X_real.hfun = @MeasurementFcn;
InfDS_X_real.par = params0;
InfDS_X_real.obsdim = Ny;
InfDS_X_real.Nsubj = Nsubj;
InfDS_X_real.pNoiseAdaptpar = [];    %Parameters to shape process noise to speed convergence. Set to empty matrix if don't want to shape process noise.
InfDS_X_real.NxNoPar = 2;	          %Number of state variables
InfDS_X_real.partQflag = 0; % A flag--if only estimating part of the state covariance matrix, set to 1; otherwise set to 0 
InfDS_X_real.Xdim = 2;
InfDS_X_real.ffun = @Lotka_Volterra_Model;
InfDS_X_real.Xbool = 1;

%SET UP PARAMETER DATA STRUCTURE AND INFO FOR PARAMETERS
InfDS_Param_real = [];
InfDS_Param_real.spkfParams = [alpha beta kappa];
InfDS_Param_real.hfun = @MeasurementFcn;
InfDS_Param_real.par = params0;
InfDS_Param_real.obsdim = Ny;
InfDS_Param_real.Nsubj = Nsubj;
InfDS_Param_real.pNoiseAdaptpar = [];    %Parameters to shape process noise to speed convergence. Set to empty matrix if don't want to shape process noise.
InfDS_Param_real.NxNoPar = 2;	          %Number of state variables
InfDS_Param_real.partQflag = 0; % A flag--if only estimating part of the state covariance matrix, set to 1; otherwise set to 0 
InfDS_Param_real.Xdim = Npar;
InfDS_Param_real.ffun = @ParamTransitionFunction;
InfDS_Param_real.Xbool = 0;


starttime = cputime;
options=optimset('Display','final','TolX',1e-7,'MaxIter',...
    10000,'MaxFunEvals',10000);


%% need to work on this 
% Maximize likelihood function Lor_loglik
% input: 
% par0: initial parameters 
% y: full ODE system dataset 
% x0: initial state variables and parameters guess 
% P0: initial covariance matrix 
% InfDs: saved information structure 
% Q: process noise matrix 
% R: measurement noise matrix 

[finalpar,fval,exitflag,output] = ...
    fminsearch('Lor_loglik',par0,options,y,x0,P0,InfDS,Q, R);
time_ML = cputime-starttime

%% 
%======== old ==========%


% June 3, 2021, P0_x and P0_params are chosen to be some identity matrix  

% P0_x_real = [20 0;0 15]; %initial covariance for states
% P0_params_real = [0.5 .001 .001 .001; .001 0.7 .001 .001;...
%     .001 .001 0.02 .001; .001 .001 .001 0.02]; %initial covariance for parameters




R_param_real.cov = R_x_real.cov; %set measurement noise (should be same as R_x_real.cov)
Q_param_real.cov = .01 * P0_params_real; %set using formula from ____

xhat_real = zeros(Nx, T); %will hold state estimates
paramshat_real = zeros(Npar, T); %will hold parameters estimates
xhat_real(:,1) = x0_real; %set first time entry
paramshat_real(:,1) = params0; %set first time entry

%===== initial guess of param 
params0 = [0.4; 0.8; 0.018; 0.023]; % alpha, gamma, beta, delta 
                                            % a1, b1, a2, b2
                                            
est_param(:,1)= params0;
                                            
for iter = 1:nloop
    iter

    
    for t = 2:T
        [param_next, P_param_next] =...
           UKF_Params(params0, P0_params_real, Q_param_real, R_param_real,...
           y(:,t), InfDS_Param_real, x0_real, Nx); %estimate parameters
        
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
    
    % alpha, gamma, beta, delta 
    est_param(:, iter+1) = mean(paramshat_real,2);
    mean_param = mean(paramshat_real,2);
    param.alpha(:,iter) = paramshat_real(1,:); 
    param.gamma(:,iter) = paramshat_real(2,:); 
            param.beta(:,iter) = paramshat_real(3,:); 
                param.delta(:,iter) = paramshat_real(4,:); 
    
    
    
    params0 = mean_param % question?
    paramshat_real(:,1) = params0;
    
    
    %pause
end 

    mean_param = mean(paramshat_real,2); 
    std_param = std(paramshat_real'); 




%PLOT REAL DATA AND PREDICTIONS ON TOP OF EACH OTHER

figure(4)
 plot(time_vector, x(1,:), '*', 'Linewidth', 2); hold on;
 plot(time_vector, xhat_real(1,:), 'blue', 'Linewidth', 2); hold off;
ylabel('Populations');
xlabel('Year');
set(gca, 'fontsize', ftsz)
legend('Prey Real', 'Prey Predicted');

figure(5)
 plot(time_vector, x(2,:), '*', 'Linewidth', 2); hold on;
 plot(time_vector, xhat_real(2,:), 'r', 'Linewidth', 2); hold off;
ylabel('Populations');
xlabel('Year');
set(gca, 'fontsize', ftsz)
legend('Predator Real', 'Predator Predicted');

% %PLOT PARAMETER ESTIMATES
% alpha_vec = truepar(1) * ones(1, T);
% gamma_vec = truepar(2) * ones(1, T);
% beta_vec = truepar(3) * ones(1, T);
% delta_vec = truepar(4) * ones(1, T);




figure

subplot(2,2,1);
plot(est_param(1,:), 'b', 'Linewidth', 2); hold on;
yline(truepar(1),':b', 'Linewidth', 2)
ylim([0,0.6])
title('\alpha');
xlabel('Iteration');
legend('Initial', 'Predicted');
set(gca, 'fontsize', ftsz)



subplot(2,2,2);
plot(est_param(2,:), 'b', 'Linewidth', 2); hold on;
yline(truepar(2),':b', 'Linewidth', 2)
ylim([0,1.2])
title('\gamma');
xlabel('Iteration');
legend('Initial', 'Predicted');
set(gca, 'fontsize', ftsz)


subplot(2,2,3);

plot( est_param(3,:), 'b', 'Linewidth', 2); hold on;
yline(truepar(3),':b', 'Linewidth', 2)
ylim([0,0.05])

title('\beta');
xlabel('Iteration');
legend('Initial', 'Predicted');
set(gca, 'fontsize', ftsz)

subplot(2,2,4);
plot(est_param(4,:), 'b', 'Linewidth', 2); hold on;
yline(truepar(3),':b', 'Linewidth', 2)
ylim([0,0.05])

title('\delta');
xlabel('Iteration');
legend('Initial', 'Predicted');
set(gca, 'fontsize', ftsz)

%======PLOT ODE WITH parameters drawn from a normal distribution 
tspan = [0 T-1];
x0 = y(:,1);


for j = 1:30
    r = normrnd(mean_param,std_param')
    %[t,sol] = ode45(@(t, y)Lotka_Volterra_Model(t, y, r),tspan, x0);
    [t,sol] = ode45(@(t, y)Lotka_Volterra_Model(t, y, mean_param'),...
        tspan, x0);
    
    figure(1)
    ylim([0,200])
    plot(0:T-1,x(1,:), '.-','LineWidth',2,'Color',[1.00,0.00,1.00],...
    'MarkerSize', 20); hold on; 
    %plot(t, sol(:,1), '-.', 'Linewidth', 1, 'Color', [0.00,0.45,0.74]);
    plot(t, sol(:,1), '-.blue', 'Linewidth', 1);
    set(gca, 'XTick', 1:5:T ,'Xticklabel', 1900: 5: 1920,'fontsize', ftsz)
    ylabel('population');
    xlabel('year');
    legend('ODE Fit', 'Raw');
    title('Prey');
    
    figure(2)
        ylim([0,200])

    plot(0:T-1,x(2,:), '.-','LineWidth',2,'Color',[1.00,0.00,1.00],...
    'MarkerSize', 20); hold on; 
    %plot(t, sol(:,2), '-.', 'Linewidth', 1, 'Color', [0.00,0.45,0.74]);
    plot(t, sol(:,2),'-.blue', 'Linewidth', 1);
    set(gca, 'XTick', 1:5:T ,'Xticklabel', 1900: 5: 1920,'fontsize', ftsz)
    ylabel('population');
    xlabel('year');
    legend('ODE Fit', 'Raw');
    title('Predator');
end




