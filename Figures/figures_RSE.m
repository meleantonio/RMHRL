%%%%%%%%%%%%%%
% Policy functions
%%%%%%%%%%%%%%%%

plot_grid_results_RSE;


%%%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 11: Plot of the policy functions on the grid/1  %
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(11);
subplot(2,2,1);
plot(grid,c1_grid_g1g1,'b--',grid,c1_grid_g1g2,'r-.',grid,c1_grid_g2g1,'k-',grid,c1_grid_g2g2,'m:','LineWidth',2);
title('consumption'); xlabel(' \theta');
legend('c1_H_H','c1_H_L','c1_L_H','c1_L_L');
subplot(2,2,2);
plot(grid,a1_grid,'b--',grid,a2_grid,'r','LineWidth',2);
title('effort'); xlabel('  \theta');
legend('a1','a2');
subplot(2,2,3);
plot(grid,phiNext1_grid_g1_g1 - grid ,'b--',grid,phiNext1_grid_g1_g2 - grid ,'r-.',grid,phiNext1_grid_g2_g1 - grid ,'k-',...
    grid,phiNext1_grid_g2_g2 - grid ,'m:','LineWidth',2);
hold on;
title('Future \theta'); xlabel('  \theta');
hold off;
legend('\theta^H^H - \theta  ',' \theta^H^L - \theta ','\theta^L^H - \theta ' , '\theta^L^L - \theta ');%,'\theta'
subplot(2,2,4);
plot(grid,lambda1_grid,'b--',grid,lambda2_grid,'r','LineWidth',2);
title('lambda'); xlabel(' \theta');
legend('lambda1','lambda2');


%%%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 12: Plot of the policy functions on the grid/2  %
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12)
subplot(2,2,1);
plot(grid,c2_grid_g1g1,'b--',grid,c2_grid_g1g2,'r-.',grid,c2_grid_g2g1,'k-',...
    grid,c2_grid_g2g2,'m:','LineWidth',2);
title('consumption'); xlabel(' \theta');
legend('c2_H_H','c2_H_L','c2_L_H','c2_L_L');
subplot(2,2,2);
plot(grid,exp_planner_disc_utility_g1,'b--',...
    grid,(c1_grid_g1g2.^(1-sig))./(1-sig) - (c1_grid_g1g1.^(1-sig))./(1-sig) + ...
    grid.*((c2_grid_g1g2.^(1-sig))./(1-sig) - (c2_grid_g1g1.^(1-sig))./(1-sig)) + ...
    exp_planner_disc_utility_g1,'r-.',...
    grid,(c1_grid_g2g1.^(1-sig))./(1-sig) - (c1_grid_g1g1.^(1-sig))./(1-sig) + ...
    grid.*((c2_grid_g2g1.^(1-sig))./(1-sig) - (c2_grid_g1g1.^(1-sig))./(1-sig)) + ...
    exp_planner_disc_utility_g1,'k-',...
    grid,(c1_grid_g2g2.^(1-sig))./(1-sig) - (c1_grid_g1g1.^(1-sig))./(1-sig) + ...
    grid.*((c2_grid_g2g2.^(1-sig))./(1-sig) - (c2_grid_g1g1.^(1-sig))./(1-sig)) + ...
    exp_planner_disc_utility_g1,'m:','LineWidth',2);%,grid,exp_planner_disc_utility_g2,'r',
title('J(y,\theta)'); xlabel(' \theta');
legend('J_H_H','J_H_L','J_L_H','J_L_L');

subplot(2,2,3);
plot(grid,exp_disc_utility1,'b--',...
    grid,(c1_grid_g1g2.^(1-sig))./(1-sig) - (c1_grid_g1g1.^(1-sig))./(1-sig) +...
    exp_disc_utility1,'r-.',...
    grid,(c1_grid_g2g1.^(1-sig))./(1-sig) - (c1_grid_g1g1.^(1-sig))./(1-sig) +...
    exp_disc_utility1,'k-',...
    grid,(c1_grid_g2g2.^(1-sig))./(1-sig) - (c1_grid_g1g1.^(1-sig))./(1-sig) +...
    exp_disc_utility1,'m:','LineWidth',2);
title('U_1(y,\theta)  '); xlabel(' \theta')

subplot(2,2,4);
plot(grid,exp_disc_utility2,'b--',...
    grid,(c2_grid_g1g2.^(1-sig))./(1-sig) - (c2_grid_g1g1.^(1-sig))./(1-sig) +...
    exp_disc_utility2,'r-.',...
    grid,(c2_grid_g2g1.^(1-sig))./(1-sig) - (c2_grid_g1g1.^(1-sig))./(1-sig) +...
    exp_disc_utility2,'k-',...
    grid,(c2_grid_g2g2.^(1-sig))./(1-sig) - (c2_grid_g1g1.^(1-sig))./(1-sig) +...
    exp_disc_utility2,'m:','LineWidth',2);
title('U_2(y,\theta)  '); xlabel(' \theta')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures of the simulated series (50000 indep. simulations of 200 periods)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate series
simul_RSE;

%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 13: simulated series, avg. allocations/1   %
%%%%%%%%%%%%%%%%%%%%%%%

figure(13);
subplot(2,2,1);
plot(1:periods_simulations,c1_simu,'b',1:periods_simulations,c2_simu,'r--','LineWidth',2);
title('consumption'); xlabel('t'); legend('c_1','c_2'); %,'avg c1','avg c2'
subplot(2,2,2);
plot(1:periods_simulations,a1_simu,'b',1:periods_simulations,a2_simu,'r--','LineWidth',2);
title('effort');xlabel('t'); legend('a_1','a_2');
subplot(2,2,3);
plot(1:periods_simulations,exp_disc_util1_simu,'b',1:periods_simulations,exp_disc_util2_simu,'r--','LineWidth',2);
title('U(y,\theta): lifetime utility (agent)');xlabel('t'); legend('U_1','U_2');
subplot(2,2,4);
plot(1:periods_simulations,lambda1_simu,'b',1:periods_simulations,lambda2_simu,'r--','LineWidth',2);
title('\lambda_1 and \lambda_2');xlabel('t')


%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 14: simulated series, avg. allocations/2   %
%%%%%%%%%%%%%%%%%%%%%%%

figure(14);
subplot(2,2,1);
plot(1:periods_simulations,phi1_simu(1:end-1),'b','LineWidth',2);
title('\theta');xlabel('t')
subplot(2,2,3);
plot(1:periods_simulations,exp_planner_disc_util_simu,'b','LineWidth',2);
title('J(y,\theta): planner value');xlabel('t')


%%%%%%%%%%%%%%%%
%    FIGURE 15: the Pareto frontier %
%%%%%%%%%%%%%%%%

gam1_vec = linspace(phi_min,phi_max,1000);

figure(15);
plot(    funeval(parexp1,fspace,gam1_vec'),funeval(parexp2,fspace,gam1_vec'),'b-',...
    'LineWidth',2);
title('Pareto Frontier');
xlabel('Agent 1 expect. utility');
ylabel('Agent 2 expect. utility');
