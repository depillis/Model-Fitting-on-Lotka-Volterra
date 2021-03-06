%AUTHOR: An Do 
% Date: 6/16/2021
% Compare estimated parameters across all the methods 

% Data for HudsonBay comes from https://gist.github.com/michaelosthege/27315631c1aedbe55f5affbccabef1ca

clear ; 
close all; 
clc; 

dataset = 'Hudson_Bay'; 
%% 

%very slow around 40-45 mins
cd ./DRAM
Run_DRAM(dataset) 

%============ Load data =================% 
cd .. 
cd ./JointUKF % Change directory
Run_JointUKF(dataset)


%change directory to PSO
cd ..
cd ./PSO
Run_PSO(dataset)


%% call dataset here 
cd ..
rawData = xlsread(['./Data/',dataset,'.xlsx']);
T = length(rawData);  %Number of time data points
x(1:2,:) = rawData(:, 2:3)'; % predator and prey populations
time_vector = rawData(:,1);


%==========Parameters ================%
UKFpar = readmatrix('./JointUKF/UKFpar.csv');
PSOpar = readmatrix('./PSO/PSOpar.csv');
DRAMpar = readmatrix('./DRAM/DRAMpar.csv');

% 2nd and 3rd entries of PSOpar and DRAM needs to swap 
PSOpar([2,3]) = PSOpar([3 2]); 
DRAMpar([2 3]) = DRAMpar([3 2]); 



%============PLOT ODE WITH FINAL PARAMETERS OVER DATA
tspan = time_vector;  %[0,T-1];
x0 = x(:,1);

sol = ode45(@(t, y)Lotka_Volterra_Model(t, y,DRAMpar), tspan, x0);
sol_DRAM= deval(sol,tspan);

sol = ode45(@(t, y)Lotka_Volterra_Model(t, y,PSOpar), tspan, x0);
sol_PSO= deval(sol,tspan);

load('./JointUKF/UKF_data.mat'); 
sol_UKF = xhat(1:2,:);

ftsz=20; 

figure(1)
    ylim([0,200])
    plot(time_vector,x(1,:), '.','MarkerSize', 20); hold on; 
    plot(time_vector,sol_UKF(1,:), '.-','LineWidth',2,'Color',[1.00,0.00,1.00],...
    'MarkerSize', 20); hold on;
    plot(time_vector, sol_DRAM(1,:), '-.', 'Linewidth', 1, 'Color', [0.00,0.45,0.74]);
    plot(time_vector, sol_PSO(1,:), '-.', 'Linewidth', 1, 'Color', [0.00,0.45,0.74]);
    
    set(gca,'Xtick', time_vector,'Xticklabel',time_vector, 'fontsize', ftsz)
    ylabel('population');
    xlabel('year');
    legend('Raw data', 'UKF','DRAM','PSO');
    title('Prey');
    
figure(2)
    ylim([0,200])
    plot(time_vector,x(2,:), '.','MarkerSize', 20); hold on; 
    plot(time_vector,sol_UKF(2,:), '.-','LineWidth',2,'Color',[1.00,0.00,1.00],...
    'MarkerSize', 20); hold on;
    plot(time_vector, sol_DRAM(2,:), '-.', 'Linewidth', 1, 'Color', [0.00,0.45,0.74]);
    plot(time_vector, sol_PSO(2,:), '-.', 'Linewidth', 1, 'Color', [0.00,0.45,0.74]);
    
    set(gca,'Xtick', time_vector,'Xticklabel',time_vector, 'fontsize', ftsz)
    ylabel('population');
    xlabel('year');
    legend('Raw data','UKF', 'DRAM','PSO');
    title('Predator');


%============PLOT THE ERROR

DRAM_error = abs(x(1:2,:) - sol_DRAM); %use matrix subtraction to get error
DRAM_error_norm =vecnorm(vecnorm(DRAM_error));

PSO_error = abs(x(1:2,:) - sol_PSO); %use matrix subtraction to get error
PSO_error_norm = vecnorm(vecnorm(PSO_error));

UKF_error = abs(x(1:2,:) - xhat(1:2,:));
UKF_error_norm= vecnorm(UKF_error);



figure(20) %Plot norm of the error over time
plot(time_vector, UKF_error(1,:)', 'Linewidth', 2); hold on;
plot(time_vector, DRAM_error(1,:)', 'Linewidth', 2); hold on; 
plot(time_vector, PSO_error(1,:)', 'Linewidth', 2); hold on; 

ylabel('Prey point-wise error ');
xlabel('Years');
legend('UKF','DRAM','PSO');
set(gca, 'fontsize', ftsz)

figure(21) %Plot norm of the error over time
plot(time_vector, UKF_error(2,:)', 'Linewidth', 2); hold on;
plot(time_vector, DRAM_error(2,:)', 'Linewidth', 2); hold on; 
plot(time_vector, PSO_error(2,:)', 'Linewidth', 2); hold on; 

ylabel('Predator point-wise error ');
xlabel('Years');
legend('UKF','DRAM','PSO');
set(gca, 'fontsize', ftsz)

