%% Find initial parameter values and ranges using fmincon
% Authors:      C. Catlett
%
% Date:         June 2020
%
% Descr:        Script to find initial parameter values for the
%               Lotka-Volterra system. Uses MATLAB fmincon


function [mse, params] = fitInitialParams(params, data)
% Extract the initial guess for your parameters (must be in vector form)
for i = 1:length(params)
    initGuess(i) = params{i}{2}; % value of parameter
    lb(i) = params{i}{3}; % lower bound of uniform prior
    ub(i) = params{i}{4}; % upper bound of uniform prior
end

%Run a quick fminsearch to get initial parameter and s2 guesses for DRAM
if length(data(:,1)) > length(params)
    fun = @(params) lotkaVolterrass(params, data);
    
%     options=optimset('Display','final','TolX',1e-7,'MaxIter',...
%     10000,'MaxFunEvals',10000);
% 
%     [theta,ss0] = ...
%         fminsearch('lotkaVolterrass',initGuess,options,data);

    %finalpar: parameters to optimize sum of square error 
    %ss0: sum of square error 
    
    opt = optimset('Display','iter','Algorithm','sqp',...
        'TolX',1e-7,'MaxIter',10000,'MaxFunEvals',10000);

    [theta, ss0] = fmincon(fun, initGuess, [], [], [], [], lb, ub, [], opt);


    mse = ss0/(length(data(:,1))-length(params));
    
    %Now replace params starting values with fminsearch output
    for i = 1:length(params)
        params{i}{2} = theta(i);
    end
else
    mse = 1; %THIS MAY NEED TO BE CHANGED! Normal guess for mse breaks down! Need more data pts than parameters
end
