% plot grid results
gridpoints = 1000;


% % grid over phi
% grid1 = linspace(phi_min,phi_max,gridpoints)';
% grid = [grid1 .5*zetabar*ones(1000,1)];

% grid over zeta
grid2 = linspace(zeta_min,zeta_max,gridpoints)';
grid = [phibar*ones(1000,1) grid2];
grid1 = grid2;

g1_grid = s(1)*ones(gridpoints,1);
g2_grid = s(2)*ones(gridpoints,1);

a1_grid  = funeval(para1,fspace,   [grid   ]);
lambda1_grid  = funeval(parlambda1,fspace,   [grid   ]);


prob1_grid_g1 =a1_grid.^nu;
prob1_grid_g2 =1 - prob1_grid_g1;


dprob1_grid_g1 =nu.*(a1_grid.^(nu-1));
dprob1_grid_g2 = - dprob1_grid_g1;


c1_grid  = funeval(parc1,fspace,grid);


phiNext1_grid_g1_g1  =  grid(:,1) + lambda1_grid.*(dprob1_grid_g1./prob1_grid_g1 );
phiNext1_grid_g1_g2  =  grid(:,1) + lambda1_grid.*(dprob1_grid_g2./prob1_grid_g2  );

 eta1_grid = grid(:,2)./betta + (  1 - grid(:,1).*(c1_grid.^(-sig)))./((-sig).*(c1_grid.^(-sig-1)))   ;


exp_planner_disc_utility_g1 = funeval(parexp_planner1,fspace ,  [grid  ]);
exp_planner_disc_utility_g2 = s(1) - s(2) + exp_planner_disc_utility_g1;%funeval(parexp_planner2,fspace ,  [grid   ]);


exp_disc_utility_g1 = funeval(parexp1,fspace ,  [grid   ]);


transfer_grid_g1 = c1_grid + g1_grid -1;
transfer_grid_g2 = c1_grid + g2_grid -1;


figure(1);
subplot(2,2,1);
plot(grid1,c1_grid,'b');
title('c1');
legend('c1 ');
subplot(2,2,2);
plot(grid1,a1_grid,'b');
title('a1');
legend('a1');
subplot(2,2,3);
plot(grid1,phiNext1_grid_g1_g1,'b',grid1,phiNext1_grid_g1_g2,'r',grid1,grid1,'k');
title('phi_t_+_1');
legend('phi_t_+_1 when low','phi_t_+_1 when high','phi_t');
subplot(2,2,4);
plot(grid1,lambda1_grid,'b');
title('lambda1');
legend('lambda1');

figure(2)
subplot(2,2,1);
plot(grid1,transfer_grid_g1,'b',grid1,transfer_grid_g2,'r');
title('transfers');
legend('transfer low','transfer high');
subplot(2,2,2);
plot(grid1,exp_planner_disc_utility_g1,'b',grid1,exp_planner_disc_utility_g2,'r');
title('planner expected utility');
legend('planner expected utility low','planner expected utility high');
subplot(2,2,3);
plot(grid1,exp_disc_utility_g1,'b');
title('expected utility agent');
legend('expected utility agent');
% subplot(2,2,4);
% plot(grid,exp_planner_disc_utility_g1,'r',grid,exp_planner_disc_utility_g2,'b',grid,-epsil*alpha*(a1_grid_g1.^(epsil-1)),'g');
% title('planner expected utility');
% legend('planner expected utility low','planner expected utility high','-v_a(a)');
subplot(2,2,4);
plot(grid1,eta1_grid,'b');
title('eta');
legend('eta');

figure(3)
subplot(2,1,1);
plot(grid1,a1_grid.^nu,'b');
title('prob1_s1');

subplot(2,1,2);
plot(grid1,1-a1_grid.^nu,'b');
title('prob1_s2');


% figure(1);
% subplot(2,1,1);
% plot(grid,c1_grid_g1,'b',grid,c1_grid_g2,'r');
% title('c1');
% legend('c1 low','c1 high');
% % subplot(2,2,2);
% % plot(1:100,c2_simu);
% % title('c2');
% % subplot(2,1,2);
% % plot(grid,leisure1_grid_g1,'b',grid,leisure1_grid_g2,'r');
% % title('leisure1');
% % legend('leisure1 low','leisure1 high');
% % % subplot(2,2,4);
% % plot(1:100,1-leisure2_simu);
% % title('labor2');
% 
% 
% figure(2);
% subplot(2,2,1);
% plot(grid,a1_grid_g1,'b',grid,a1_grid_g2,'r');
% title('a1');
% legend('a1 low','a1 high');
% 
% % subplot(2,2,2);
% % plot(1:100,a2_simu);
% % title('a2');
% subplot(2,2,2);
% plot(grid,lambda1_grid_g1,'b',grid,lambda1_grid_g2,'r');
% title('lambda1');
% legend('lambda1 low','lambda1 high');
% 
% % subplot(2,2,3);
% % plot(grid,value_grid_g1,'b',grid,value_grid_g2,'r');
% % title('value');
% % legend('value low','value high');
% % 
% 
% figure(3);
% subplot(2,2,1);
% plot(grid,phiNext1_grid_g1_g1,'b',grid,phiNext1_grid_g1_g2,'r',grid,phiNext1_grid_g2_g1,'g',grid,phiNext1_grid_g2_g2,'y');
% title('phi_t_+_1');
% legend('phi_t_+_1 low when low','phi_t_+_1 high when low','phi_t_+_1 low when high','phi_t_+_1 high when high');
% 
% % subplot(2,2,2);
% % plot(grid,eta_grid_g1,'b',grid,eta_grid_g2,'r');
% % title('eta');
% % legend('eta low','eta high');
% 
% % subplot(2,2,3);
% % plot(grid,eta_grid_g1,'b',grid,eta_grid_g2,'r');
% % title('eta');
% % legend('eta low','eta high');
% % subplot(2,2,4);
% % plot(1:periods_simulations,eta_simu.*(c1_simu.^(sig))+ gam1_scal*(1:periods_simulations),'r', 1:periods_simulations, eta_simu.*(c1_simu.^(sig)),'b');
% % title('eta_t/u_c_t');
% 
% 
% 
% 
% 
% 
% % figure(4);
% % subplot(2,1,1);
% % plot(1:50,c1_simu(1:50));
% % title('c1');
% % % subplot(2,2,2);
% % % plot(1:50,c2_simu(1:50));
% % % title('c2');
% % subplot(2,1,2);
% % plot(1:50,1-leisure1_simu(1:50));
% % title('labor1');
% % % subplot(2,2,4);
% % % plot(1:50,1-leisure2_simu(1:50));
% % % title('labor2');
% 
% 
% % figure(5);
% % subplot(2,1,1);
% % plot(1:50,a1_simu(1:50));
% % title('a1');
% % % subplot(2,2,2);
% % % plot(1:50,a2_simu(1:50));
% % % title('a2');
% % subplot(2,1,2);
% % plot(1:50,lambda1_simu(1:50));
% % title('lambda1');
% % % subplot(2,2,4);
% % % plot(1:50,lambda2_simu(1:50));
% % % title('lambda2');
% % 
% % figure(6);
% % % subplot(2,2,1);
% % plot(1:50,phi1_simu(1:50));
% % title('phi1');
% % % subplot(2,2,2);
% % % plot(1:50,phi2_simu(1:50));
% % % title('phi2');
% % % subplot(2,2,3);
% % % plot(1:100,g1_simu(1:50));
% % % title('g1');
% % % subplot(2,2,4);
% % % plot(1:100,g2_simu(1:50));
% % % title('g2');
% 
% 
% % figure(7);
% % plot(1:periods_simulations-1,eta_simu(2:end).*((c1_simu(2:end)).^(sig)) - ...
% %     eta_simu(1:end-1).*((c1_simu(1:end-1)).^(sig)) ,'r');
% % title('eta_t_+_1/u_c_,_t_+_1 -eta_t/u_c_t  ');
% 
% 
% 
% 
% 
% % % show first realization
% % 
% % figure(11);
% % subplot(2,1,1);
% % plot(1:periods_simulations,c1_simu_real(1,:));
% % title('c1');
% % % subplot(2,2,2);
% % % plot(1:100,c2_simu_real(1,:));
% % % title('c2');
% % subplot(2,1,2);
% % plot(1:periods_simulations,1-leisure1_simu_real(1,:));
% % title('labor1');
% % % subplot(2,2,4);
% % % plot(1:100,1-leisure2_simu_real(1,:));
% % % title('labor2');
% % 
% % 
% % figure(12);
% % subplot(2,1,1);
% % plot(1:periods_simulations,a1_simu_real(1,:));
% % title('a1');
% % % subplot(2,2,2);
% % % plot(1:100,a2_simu_real(1,:));
% % % title('a2');
% % subplot(2,1,2);
% % plot(1:periods_simulations,lambda1_simu_real(1,:));
% % title('lambda1');
% % subplot(2,2,3);
% % plot(1:periods_simulations,value_simu_real(1,:));
% % title('planner value function');
% % 
% % figure(13);
% % subplot(2,1,1);
% % plot(1:periods_simulations,phi1_simu_real(1,1:end-1));
% % title('phi1');
% % % subplot(2,2,2);
% % % plot(1:100,phi2_simu_real(1,1:end-1));
% % % title('phi2');
% % subplot(2,1,2);
% % plot(1:periods_simulations,g1_simu(1,:));
% % title('g1');
% % % subplot(2,2,4);
% % % plot(1:100,g2_simu(1,:));
% % % title('g2');
% % 
% % figure(14);
% % subplot(2,2,1);
% % plot(1:periods_simulations,gam1_simu_real(1,:),'b');
% % title('gam1');
% % subplot(2,2,2);
% % plot(1:periods_simulations,eta_simu_real(1,:),'r');
% % title('eta');
% % subplot(2,2,3);
% % plot(1:periods_simulations-1,eta_simu_real(1,2:end)./eta_simu_real(1,1:end-1),'r');
% % title('eta_t_+_1/eta_t');
% % subplot(2,2,4);
% % plot(1:periods_simulations,eta_simu_real(1,:).*((c1_simu_real(1,:)).^(sig)),'r');
% % title('eta_t/u_c_t');
% % 
% % 
% % 
% % 
% % figure(15);
% % subplot(2,1,1);
% % plot(1:periods_simulations,g1_simu(1,:),'b',1:periods_simulations,c1_simu_real(1,:),'g');
% % title('c1 and g1');
% % legend('g1','c1');
% % subplot(2,1,2);
% % plot(1:periods_simulations-1,eta_simu_real(1,2:end).*((c1_simu_real(1,2:end)).^(sig)) - ...
% %     eta_simu_real(1,1:end-1).*((c1_simu_real(1,1:end-1)).^(sig)),'r');
% % title('eta_t_+_1/u_c_,_t_+_1 -eta_t/u_c_t  ');
% 
% 
% 
% 
% % figure(14);
% % subplot(2,2,1);
% % plot(1:50,c1_simu(1:50));
% % title('c1');
% % subplot(2,2,2);
% % plot(1:50,c2_simu(1:50));
% % title('c2');
% % subplot(2,2,3);
% % plot(1:50,1-leisure1_simu(1:50));
% % title('labor1');
% % subplot(2,2,4);
% % plot(1:50,1-leisure2_simu(1:50));
% % title('labor2');
% % 
% % 
% % figure(15);
% % subplot(2,2,1);
% % plot(1:50,a1_simu(1:50));
% % title('a1');
% % subplot(2,2,2);
% % plot(1:50,a2_simu(1:50));
% % title('a2');
% % subplot(2,2,3);
% % plot(1:50,lambda1_simu(1:50));
% % title('lambda1');
% % subplot(2,2,4);
% % plot(1:50,lambda2_simu(1:50));
% % title('lambda2');
% % 
% % figure(16);
% % subplot(2,1,1);
% % plot(1:50,phi1_simu(1:50));
% % title('phi1');
% % subplot(2,1,2);
% % plot(1:50,phi2_simu(1:50));
% % title('phi2');
% % 
% 
% 
% 
% 
% 
% 
% 
% 
% 
