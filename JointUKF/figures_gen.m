
%========== PLOTS ================%

finalpar = xhat(3:6,end); 

time_vector = rawData(:,1);
ftsz=20;
%============ PLOT REAL DATA AND PREDICTIONS ON TOP OF EACH OTHER
figure(1)
 plot(time_vector, x(1,:), '*', 'Linewidth', 2); hold on;
 plot(time_vector, xhat(1,:), 'r', 'Linewidth', 2); hold off;
ylabel('Populations');
xlabel('Year');
set(gca, 'fontsize', ftsz)
legend('Prey Real', 'Prey Predicted');

figure(2)
 plot(time_vector, x(2,:), '*', 'Linewidth', 2); hold on;
 plot(time_vector, xhat(2,:), 'r', 'Linewidth', 2); hold off;
ylabel('Populations');
xlabel('Year');
set(gca, 'fontsize', ftsz)
legend('Predator Real', 'Predator Predicted');



%============PLOT THE ERROR
x_error = x(1:2,:) - xhat(1:2,:); %use matrix subtraction to get error
error_norm = vecnorm(x_error); %take columnwise norm

figure(20) %Plot norm of the error over time
plot(time_vector, error_norm(1,:), 'b', 'Linewidth', 2);
ylabel('Error Norm');
xlabel('Time');
title('Norm of Error Over Time');
set(gca, 'fontsize', ftsz)
% 
%============PLOT ODE WITH FINAL PARAMETERS OVER DATA
finalpar=xhat(3:6,end);
tspan = [0 T-1];
x0 = y(:,1);
[t,sol] = ode45(@(t, y)Lotka_Volterra_Model(t, y, finalpar), tspan, x0);





figure(30)
    ylim([0,200])
    plot(0:T-1,x(1,:), '.-','LineWidth',2,'Color',[1.00,0.00,1.00],...
    'MarkerSize', 20); hold on; 
    plot(t, sol(:,1), '-.', 'Linewidth', 1, 'Color', [0.00,0.45,0.74]);
    %plot(t, sol(:,1), '-blue', 'Linewidth', 1);
    set(gca, 'XTick', 1:5:T ,'Xticklabel', 1900: 5: 1920,'fontsize', ftsz)
    ylabel('population');
    xlabel('year');
    legend('ODE Fit', 'Raw');
    title('Prey');
    
figure(31)
    ylim([0,200])
    plot(0:T-1,x(2,:), '.-','LineWidth',2,'Color',[1.00,0.00,1.00],...
    'MarkerSize', 20); hold on; 
    plot(t, sol(:,2), '-.', 'Linewidth', 1, 'Color', [0.00,0.45,0.74]);
    %plot(t, sol(:,2),'-.red', 'Linewidth', 1);
    set(gca, 'XTick', 1:5:T ,'Xticklabel', 1900: 5: 1920,'fontsize', ftsz)
    ylabel('population');
    xlabel('year');
    legend('ODE Fit', 'Raw');
    title('Predator');
% 
% %====== Estimated ODE plot using Mahafy's parameters
% 
% params0 = [0.4; 0.8; 0.018; 0.023]; % alpha, gamma, beta, delta 
% 
% [t,mahafy] = ode45(@(t, y) ...
%     Lotka_Volterra_Model(t, y, params0), tspan, x0);
% 
% figure(101)
% subplot(2,1,1);
% plot(0:T-1,x(1,:), '*black'); hold on; 
% plot(t, mahafy(:,1), 'blue', 'Linewidth', 2);
% set(gca, 'XTick', 1:5:T ,'Xticklabel', 1900: 5: 1920,'fontsize', ftsz)
% hold off; 
% ylabel('population');
% xlabel('year');
% legend('Mahafy data','ODE fit using Mahafy parameters');
% title('Prey');
% 
% subplot(2,1,2);
% plot(0:T-1,  x(2,:), '*black'); hold on; 
% set(gca, 'XTick', 1:5:T ,'Xticklabel', 1900: 5: 1920,'fontsize', ftsz)
% plot(t, mahafy(:,2), 'red', 'Linewidth', 2);
% hold off; 
% ylabel('population');
% xlabel('year');
% legend('Mahafy data','ODE fit using Mahafy parameters');
% title('Predator');
% set(gca, 'fontsize', ftsz)
