function Joint_AD(dataset)
% July 10 2021 An modify this routine into a function
% input: dataset e.g 'Mahafy' or 'Hudson_Bay'
% output: none. 
% All data are saved in UKFpar.csv and UKF_data.csv and UKF_time.csv


% We started with the REU2020 student code   
% we did not know where they find the initial parameters from 
% the dataset they worked on from Shtylla but we also lost track of their
% source 
% An modified 
% We obtained 2 dataset: 
% 1. mahafy https://jmahaffy.sdsu.edu/courses/f09/math636/lectures/lotka/qualde2.html
% 2. Hudson_bay: https://gist.github.com/michaelosthege/27315631c1aedbe55f5affbccabef1ca
% Intial guesses are chosen to be the same as Mahafy's 
% Covariance matrix is WAY MORE SIMPLER
% measurement covariance matrix
% UKF is extremely sensitive to initial guesses btw
% 

tic 
ftsz= 20; % fontsize
close all 
 
% INITIALISATION AND PARAMETERS:
% ==============================

Ny = 2;  %Number of observables
Nx = 6;  %Number of values in x vector (# states + # parameters)
Nstates = 2; %number states
Npar = 4;   %number parameters
Nsubj = 1; %Number of subjects to run the UKF on, in our case unused

%************** START HERE *********************************

% Load dataset of interest
% extract time data, prey and predator data (in units of thousands)

HLData = load('HaresLynxData_AD.mat');  
rawData = eval(['HLData.',dataset]);
T = length(rawData);  %Number of time data points
x(1:2,:) = rawData(:, 2:3)'; % predator and prey populations

truepar = [0.4807; 0.9272; 0.02482; 0.02756];% alpha, gamma, beta, delta 


%================UKF PARAMETERS, adapted from Chow et al. 
alpha = 10e-3; % UKF : point scaling parameter
beta = 2;  % UKF : scaling parameter for higher order terms of Taylor series expansion 
kappa = 0;   % UKF : sigma point selection scaling parameter (best to leave this = 0)


% allocation for observable data
% MeasurementFcn should add some noise to the data 
% if more infomration is known 
% Currently, it does nothing. 
% but we still keep in here for future use

y = zeros(Ny, T); 
Q = []; %initialize struct


 %Generate measurements
  for t = 1:T
      y(:,t) = feval('MeasurementFcn', x(:,t), [], [], [], []); %Call measurement function
  end

%==== Initiate parameters and observables 
x0_real = y(:,1);
par0 = [0.4 ;0.8; 0.02; 0.02]; %truepar; % alpha, gamma, beta, delta 


% June 7: An Do 
% REU student original code reads: 
% P0 = [20 .01 .01 .01 .01 .01;
%       .01 15 .01 .01 .01 .01;
%       .01 .01 0.5 .01 .01 .01;
%       .01 .01 .01 0.3 .01 .01;
%       .01 .01 .01 .01 0.02 .01;
%       .01 .01 .01 .01 .01 0.02];
% An modified it to just identity matrix as in Chow et al.
  
P0 = .01*eye(Nx); % diagonal matrix of entries =1 

no_of_runs = 1; % the algorithm is deterministic, 
                % multiple runs will give the same result  
time_vector = rawData(:,1);
ftsz=20;

for j=1:no_of_runs
    j 

% Chow's code
% rand('state',sum(100*clock));   % Shuffle the pack!
% randn('state',sum(100*clock));   % Shuffle the pack!


%intial state + parameters
x0 = [x0_real; par0]; 

fprintf('Now estimating the model parameters...')
fprintf('\n')

% Chow et al. 
Q.cov = 10e-5*eye(Nx); % Add a little bit of process
R = []; %Initialize covariance of measurement noise as struct
%Set R covariance matrix using known parameter values for noise (determined ad hoc)
R.cov = diag([30,15]); %diag(x0(1:2)); %Set R covariance matrix using known parameter values for noise (determined ad hoc)
% 



%Set up struct of information needed to run the UKF
InfDS.par = par0; % initial values for model parameters
InfDS.spkfParams = [alpha beta kappa]; %UKf scaling parameters
InfDS.hfun = @MeasurementFcn; %The measurement function H
InfDS.obsdim = Ny; %Number of observables
InfDS.Nsubj = Nsubj;
InfDS.pNoiseAdaptpar = [];    %Parameters to shape process noise to speed convergence. Set to empty matrix if don't want to shape process noise.
InfDS.NxNoPar = 2;	          %Number of state variables
InfDS.partQflag = 0;	% A flag--if only estimating part of the state covariance matrix, set to 1; otherwise set to 0 
InfDS.Xdim = Nx; %x dimension 3 for 3 states

%covariance of measurement noise
starttime = cputime; %Initialize starttime (to calculate how long UKF takes to run)

%== test UKF routine to see if it works before we start fmin routine ===%
% the UKF routine alone takes about 18 seconds
[xhat, Px, xh_, Px_, yh_, inov, Py, KG,negloglik, Pxall, Pyall, Px_all]= ...
    UKF(x0,P0,Q, R, y,[],[],InfDS); 
time_ukf = cputime-starttime %Time that UKF took

%======== Minimize negative log likelihood routine =============% 
% GOAL: minimize the "loglik.m" function with respect to the intial guess
% 
% input: 
%       par0: initial parameters 
%       y: full dataset 
%       x0: initial state variables and parameters guess 
%       P0: initial covariance matrix 
%       InfDS: saved information structure 
%       Q: process noise matrix 
%       R: measurement noise matrix 
% output:   
%       finalpar: optimal x0 to minimize the cost function 
% take about 18 mins 
% 
% options=optimset('Display','final','TolX',1e-7,'MaxIter',...
%     10000,'MaxFunEvals',10000);
% 
% [finalpar,fval,exitflag,output] = ...
%     fminsearch('loglik',x0,options,y,P0,InfDS,Q, R);
% time_ML = cputime-starttime
% 
% % update measurement error matrix with optimal initial observables 
% % update initial parameters with optimal initial parameters 
% % x0 = finalpar;
% R.cov = diag(finalpar(1:2)); % diag(diag([repmat(10e-4,Ny)])) %
% InfDS.par = finalpar(3:6); 
% 
% %Re-run UKF routine with the updated information 
% starttime = cputime;
% [xhat, Px, xh_, Px_, yh_, inov, Py, KG, negloglik,Pxall, Pyall, Px_all]=...
%     UKF(x0,P0,Q, R, y,[],[],InfDS); %Use UKF to do state estimates; 
% xhat(:,1)=x0;

%===================== PLOTS============================%

figure(1)
 plot(time_vector, x(1,:), '*', 'Linewidth', 2); hold on;
 plot(time_vector, xhat(1,:), '--', 'Linewidth', 2); 
ylabel('Populations');
xlabel('Year');
set(gca, 'fontsize', ftsz)
legend('Prey Real', 'Prey Predicted');

figure(2)
 plot(time_vector, x(2,:), '*', 'Linewidth', 2); hold on;
 plot(time_vector, xhat(2,:), '--', 'Linewidth', 2); 
ylabel('Populations');
xlabel('Year');
set(gca, 'fontsize', ftsz)
legend('Predator Real', 'Predator Predicted');


% estimated parameters 
par.alpha(j,:) = xhat(3,:); 
par.gamma(j,:) = xhat(4,:); 
par.beta(j,:) = xhat(5,:); 
par.delta(j,:) = xhat(6,:); 

 figure(3)
    subplot(2,2,1);
    plot(par.alpha(j,:),'.-b', 'Linewidth', 2, 'MarkerSize',20);
    hold on
    yline(truepar(1), 'Linewidth', 2)
    title('\alpha');
    xlabel('years');
    set(gca, 'fontsize', ftsz)

    subplot(2,2,2);
    plot(par.gamma(j,:), '.-b', 'Linewidth', 2, 'MarkerSize',20); 
    hold on
    yline(truepar(2), 'Linewidth', 2)
    title('\gamma');
    xlabel('years');
    set(gca, 'fontsize', ftsz)

    subplot(2,2,3);
    plot(par.beta(j,:), '.-b', 'Linewidth', 2, 'MarkerSize',20); 
    hold on
    yline(truepar(3), 'Linewidth', 2)
    title('\beta');
    xlabel('years');
    set(gca, 'fontsize', ftsz)

    subplot(2,2,4);
    plot(par.delta(j,:),'.-b', 'Linewidth', 2, 'MarkerSize',20); 
    hold on
    yline(truepar(4), 'Linewidth', 2)
    title('\delta');
    xlabel('years');
    set(gca, 'fontsize', ftsz)
% 
% %Calsulate mean squared error for prey
% prey_error = xhat(1,:) - x(1,:);
% prey_error = prey_error.^2;
% %prey_error = sum(prey_error(:,45:91));
% prey_error = sum(prey_error);
% MSE_prey(:,j) = prey_error /T;
% 
% %Calculate mean squared error for predator
% predator_error = xhat(2,:) - x(2,:);
% predator_error = predator_error.^2;
% %predator_error = sum(predator_error(:, 45:91));
% 
% predator_error = sum(predator_error);
% 
% MSE_predator(:,j) = predator_error / T;
end

x_error = x(1:2,:) - xhat(1:2,:); %use matrix subtraction to get error
error_norm = vecnorm(x_error); %take columnwise norm

save('UKF_data.mat','xhat','error_norm')
csvwrite('UKFpar.csv',xhat(3:end,end))

elapse = toc; 
csvwrite('UKF_time.csv', elapse)

% param_final %Print out final parameters
% MSE_prey %Print out prey mean-squared error
% MSE_predator %Print out predator mean-squared error

% ======= Generate all figures =========% 
%figures_gen

hold off; 

end

