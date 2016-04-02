
% Figures of the simulated series (50000 indep. simulations of 200 periods)

% Simulate series
simul_RSP;



%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 16/18: simulated series, avg. allocations/1   %
%%%%%%%%%%%%%%%%%%%%%%%

figure(16+avg);
subplot(2,2,1);
plot(1:periods_simulations,c1_simu,'b',1:periods_simulations,c2_simu,'r:','LineWidth',1);
title('consumption'); xlabel('t'); legend('c_1','c_2'); %,'avg c1','avg c2'
subplot(2,2,2);
plot(1:periods_simulations,a1_simu,'b',1:periods_simulations,a2_simu,'r:','LineWidth',2);
title('effort');xlabel('t'); legend('a_1','a_2');
subplot(2,2,3);
plot(1:periods_simulations,exp_disc_util1_simu,'b',1:periods_simulations,exp_disc_util2_simu,'r:','LineWidth',2);
title('U(y,\theta)');xlabel('t'); legend('U_1','U_2');
subplot(2,2,4);
plot(1:periods_simulations,lambda1_simu,'b',1:periods_simulations,lambda2_simu,'r:','LineWidth',2);
title('lambda');xlabel('t'); legend('\lambda_1','\lambda_2');


%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 17/19: simulated series, avg. allocations/2   %
%%%%%%%%%%%%%%%%%%%%%%%

figure(17+avg);
subplot(2,2,1);
plot(1:periods_simulations,phi1_simu(1:end-1),'b','LineWidth',2);
title('\theta');xlabel('t');
subplot(2,2,2);
plot(1:periods_simulations,k1next_simu(1:end-1),'b',1:periods_simulations,k2next_simu(1:end-1),'r:','LineWidth',2);
title('capital');xlabel('t'); legend('k_1','k_2');
subplot(2,2,3);
plot(1:periods_simulations,exp_planner_disc_util_simu,'b','LineWidth',2);
title('J(y,\theta)');xlabel('t')
subplot(2,2,4);
plot(1:periods_simulations,invest1_simu,'b',1:periods_simulations,invest2_simu,'r:','LineWidth',2);
title('investment');xlabel('t'); legend('i_1','i_2');



