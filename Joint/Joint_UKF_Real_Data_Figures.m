function [] = Joint_UKF_Real_Data_Figures(xhat, x, T, truepar, y, param_final)


ftsz= 25; % fontsize
time_vector = rawData(:,1);



%PLOT REAL DATA AND PREDICTIONS ON TOP OF EACH OTHER
figure
 plot(time_vector, x(1,:), '*', 'Linewidth', 2); hold on;
 plot(time_vector, xhat_real(1,:), 'r', 'Linewidth', 2); hold off;
ylabel('Populations');
xlabel('Year');
set(gca, 'fontsize', ftsz)
legend('Prey Real', 'Prey Predicted');

figure
 plot(time_vector, x(2,:), '*', 'Linewidth', 2); hold on;
 plot(time_vector, xhat_real(2,:), 'r', 'Linewidth', 2); hold off;
ylabel('Populations');
xlabel('Year');
set(gca, 'fontsize', ftsz)
legend('Predator Real', 'Predator Predicted');



%PLOT PARAMETER ESTIMATES
alpha_vec = truepar(1) * ones(1, T);
gamma_vec = truepar(2) * ones(1, T);
beta_vec = truepar(3) * ones(1, T);
delta_vec = truepar(4) * ones(1, T);


figure (11)
subplot(2,2,1);
plot(time_vector, alpha_vec, 'k+', 'Linewidth', 2); hold on;
plot(time_vector, xhat(3,:), 'b', 'Linewidth', 2); hold off;
ylabel('\alpha');
xlabel('year');
legend('Initial', 'Predicted');
title('\alpha Parameter Prediction');
set(gca, 'fontsize', ftsz)


subplot(2,2,2);

plot(time_vector, gamma_vec, 'k+', 'Linewidth', 2); hold on;
plot(time_vector, xhat(4,:), 'b', 'Linewidth', 2); hold off;
ylabel('\gamma');
xlabel('year');
legend('Initial', 'Predicted');
title('\gamma Parameter Prediction');
set(gca, 'fontsize', ftsz)




subplot(2,2,3);
p1 = plot(time_vector, beta_vec, 'k+', 'Linewidth', 2); hold on;
p2 = plot(time_vector, xhat(5,:), 'b', 'Linewidth', 2); hold off;
ylabel('\beta');
xlabel('year');
legend('Initial', 'Predicted');
title('\beta Parameter Prediction');
set(gca, 'fontsize', ftsz)


subplot(2,2,4);
p1 = plot(time_vector, delta_vec, 'k+', 'Linewidth', 2); hold on;
p2 = plot(time_vector, xhat(6,:), 'b', 'Linewidth', 2); hold off;
ylabel('\delta');
xlabel('year');
legend('Initial', 'Predicted');
title('\delta Parameter Prediction');
set(gca, 'fontsize', ftsz)



%PLOT THE ERROR
x_error = x(1:2,:) - xhat_real(1:2,:); %use matrix subtraction to get error
error_norm = vecnorm(x_error); %take columnwise norm

figure(20) %Plot norm of the error over time
plot(time_vector, error_norm(1,:), 'b', 'Linewidth', 2);
ylabel('Error Norm');
xlabel('Time');
title('Norm of Error Over Time');
set(gca, 'fontsize', ftsz)

%PLOT ODE WITH FINAL PARAMETERS OVER DATA
tspan = [0 T-1];
x0 = y(:,1);
[t,sol] = ode45(@(t, y) ...
    Lotka_Volterra_Model(t, y, param_final), tspan, x0);


figure(100)
subplot(2,1,1);
plot(0:T-1,x(1,:), '*black'); hold on; 
plot(t, sol(:,1), 'blue', 'Linewidth', 2);
set(gca, 'XTick', 1:5:T ,'Xticklabel', 1900: 5: 1920,'fontsize', ftsz)
hold off; 
ylabel('population');
xlabel('year');
legend('ODE Fit', 'Raw');
title('Prey');

subplot(2,1,2);
plot(0:T-1,  x(2,:), '*black'); hold on; 
set(gca, 'XTick', 1:5:T ,'Xticklabel', 1900: 5: 1920,'fontsize', ftsz)
plot(t, sol(:,2), 'red', 'Linewidth', 2);
hold off; 
ylabel('population');
xlabel('year');
legend('ODE Fit', 'Raw');
title('Predator');
set(gca, 'fontsize', ftsz)



% 
% subplot(2,1,1);
% p1 = plot(1845:1935, x_with_final(1,:), 'b', 'Linewidth', 2); hold on;
% p2 = plot(1845:1935, x(1,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold off;
% p2.MarkerSize = 3;
% ylabel('population (thousands)', 'fontsize',15);
% xlabel('year', 'fontsize', 15);
% legend('ODE Fit', 'Raw');
% title('Prey');
% 
% subplot(2,1,2);
% p1 = plot(1845:1935, x_with_final(2,:), 'r', 'Linewidth', 2); hold on;
% p2 = plot(1845:1935, x(2,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold off;
% p2.MarkerSize = 3;
% ylabel('population (thousands)', 'fontsize',15);
% xlabel('year', 'fontsize', 15);
% legend('ODE Fit', 'Raw');
% title('Predator');


% %CALCULATE MODEL ERROR WITH FINAL PARAMETERS
% prey_final_error = x(1,:) - x_with_final(1,:);
% prey_final_error = prey_final_error.^2;
% prey_final_error = sum(prey_final_error);
% MSE_prey_final = prey_final_error /91;
% MSE_prey_final
% 
% predator_final_error = x(2,:) - x_with_final(2,:);
% predator_final_error = predator_final_error.^2;
% predator_final_error = sum(predator_final_error);
% MSE_predator_final = predator_final_error /91;
% MSE_predator_final
% 
% %CALCULATE MODEL ERROR WITH STARTING PARAMETERS
% prey_start_error = x(1,:) - x_start(1,:);
% prey_start_error = prey_start_error.^2;
% prey_start_error = sum(prey_start_error);
% MSE_prey_start = prey_start_error /91;
% MSE_prey_start
% 
% predator_start_error = x(2,:) - x_start(2,:);
% predator_start_error = predator_start_error.^2;
% predator_start_error = sum(predator_start_error);
% MSE_predator_start = predator_start_error / 91;
% MSE_predator_start


%=======================

params0 = [0.4; 0.8; 0.018; 0.023]; % alpha, gamma, beta, delta 

[t,sol] = ode45(@(t, y) ...
    Lotka_Volterra_Model(t, y, params0), tspan, x0);
% 
% 
% 
% sol = ode45(@(t, y) Lotka_Volterra_Model(t, y, params0), tspan, x0); %Use ODE solver
% x_start = zeros(6, 91);
% for j = 1:91
%     x_start(:,j) = deval(sol, j);
% end

% figure(102)
% p1 = plot(1:T, x_start(1,:), 'b', 'Linewidth', 2); hold on;
% p2 = plot(1:T, x(1,:), 'b+', 'Linewidth', 2); hold off;
% ylabel('Prey Population', 'fontsize', 15);
% xlabel('Year', 'fontsize', 15);
% title('Prey Simulated Using Initial Paramaters');
% 
% figure(103)
% p1 = plot(1:T, x_start(2,:), 'b', 'Linewidth', 2); hold on;
% p2 = plot(1:T, x(2,:), 'b+', 'Linewidth', 2); hold off;
% ylabel('Predator Population', 'fontsize', 15);
% xlabel('Year', 'fontsize', 15);
% title('Predator Simulated Using Initial Paramaters');


figure(101)
subplot(2,1,1);
plot(0:T-1,x(1,:), '*black'); hold on; 
plot(t, sol(:,1), 'blue', 'Linewidth', 2);
set(gca, 'XTick', 1:5:T ,'Xticklabel', 1900: 5: 1920,'fontsize', ftsz)
hold off; 
ylabel('population');
xlabel('year');
legend('ODE Fit with Mahafy estimated value ', 'Raw');
title('Prey');

subplot(2,1,2);
plot(0:T-1,  x(2,:), '*black'); hold on; 
set(gca, 'XTick', 1:5:T ,'Xticklabel', 1900: 5: 1920,'fontsize', ftsz)
plot(t, sol(:,2), 'red', 'Linewidth', 2);
hold off; 
ylabel('population');
xlabel('year');
legend('ODE Fit with Mahafy estimated value ', 'Raw');
title('Predator');
set(gca, 'fontsize', ftsz)





