%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Policy functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% the grid
phigrid_nodes_test = linspace(phi_min,phi_max,33)';%
zetagrid_nodes_test = linspace(zeta_min,zeta_max,33)';%
k1_test = zetagrid_nodes_test;
phi_test = phigrid_nodes_test;

% generate policy functions on the grid
for uu =1:length(phi_test)
    for hh=1: length(k1_test)
        zed1(uu,hh) = funeval(parc1,fspace,[phi_test(uu) k1_test(hh)]);
        zed2(uu,hh) = funeval(parexp_planner1,fspace,[phi_test(uu) k1_test(hh)]);
        zed3(uu,hh) = funeval(para1,fspace,[phi_test(uu) k1_test(hh)]);
        zed4(uu,hh) = funeval(parexp1,fspace,[phi_test(uu) k1_test(hh)]);
        zed5(uu,hh) = funeval(parexp_planner1,fspace,[phi_test(uu) k1_test(hh)]) ...
            - phi_test(uu).*zed4(uu,hh) ;
        zed55(uu,hh) = funeval(parexp_planner1,fspace,[phi_test(uu) k1_test(hh)])  ;
        zed6(uu,hh) = funeval(parlambda1,fspace,[phi_test(uu) k1_test(hh)]) ;
        prob1(uu,hh) = zed3(uu,hh).^nu;
        prob2(uu,hh) = 1 - zed3(uu,hh).^nu;
        dprob1(uu,hh) = nu.*(zed3(uu,hh).^(nu-1));
        dprob2(uu,hh) =  - dprob1(uu,hh);
        zed8(uu,hh) = phi_test(uu)+...
            + zed6(uu,hh).*dprob1(uu,hh)./prob1(uu,hh)  ;
        zed9(uu,hh) = phi_test(uu)+...
            + zed6(uu,hh).*dprob2(uu,hh)./prob2(uu,hh)  ;
        zed10(uu,hh) = phi_test(uu);
        zed11(uu,hh)  =   k1_test(hh)./betta + (  1 - ...
            phi_test(uu).*((zed1(uu,hh)).^(-sig)))./((-sig).*((zed1(uu,hh).^(-sig-1))))   ;
    end
end
zed12 = zed11;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6: Policy functions/1
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6);
subplot(2,2,1);
surf(k1_test,phi_test,zed1);
title('consumption'); ylabel(' \phi'); xlabel(' \zeta');
axis([zeta_min zeta_max phi_min phi_max min(min(zed1)) max(max(zed1))]);
subplot(2,2,2);
surf(k1_test,phi_test,zed11);
title('\eta'); ylabel(' \phi'); xlabel('\zeta');
axis([zeta_min zeta_max phi_min phi_max min(min(zed11)) max(max(zed11))]);
subplot(2,2,3);
surf(k1_test,phi_test,zed3);
title('effort'); ylabel(' \phi'); xlabel(' \zeta');
axis([zeta_min zeta_max phi_min phi_max min(min(zed3)) max(max(zed3))]);
subplot(2,2,4);
surf(k1_test,phi_test,zed4);
title('U(y,\phi,\zeta):agent continuation value ');ylabel(' \phi'); xlabel(' \zeta');
axis([zeta_min zeta_max phi_min phi_max min(min(zed4)) max(max(zed4))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7: Policy functions/2
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7);
subplot(2,2,1);
surf(k1_test,phi_test,zed9);
title('Pareto weights: Low'); ylabel(' \phi'); xlabel(' \zeta');
axis([zeta_min zeta_max phi_min phi_max min(min(zed9)) max(max(zed8))]);
subplot(2,2,2);
surf(k1_test,phi_test,zed8);
title('Pareto weights: High'); ylabel(' \phi'); xlabel(' \zeta');
axis([zeta_min zeta_max phi_min phi_max min(min(zed9)) max(max(zed8))]);
subplot(2,2,3);
surf(k1_test,phi_test,zed6);
title('\lambda'); ylabel(' \phi'); xlabel(' \zeta');
axis([zeta_min zeta_max phi_min phi_max min(min(zed6)) max(max(zed6))]);
subplot(2,2,4);
surf(k1_test,phi_test,zed55);
title('J(y,\phi,\zeta): planner value ');ylabel(' \phi'); xlabel(' \zeta');
axis([zeta_min zeta_max phi_min phi_max min(min(zed55)) max(max(zed55))]);


%%%%%%%%%%%%%%%%
%    FIGURE 10: the Pareto frontier %
%%%%%%%%%%%%%%%%
figure(10);
plot(zed4(:,1),zed5(:,1),'b',zed4(:,14),zed5(:,14),'r',zed4(:,30),zed5(:,30),'g','LineWidth',2);
title('Pareto frontier '); legend('\zeta_0 =0','\zeta_0 =0.1774','\zeta_0 =0.3825');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figures of the simulated series (50000 indep. simulations of 200 periods)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate series
simul_HA;

%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 8: simulated series, avg. allocations/1   %
%%%%%%%%%%%%%%%%%%%%%%%

figure(8);
subplot(2,2,1);
plot(1:periods_simulations,phi1_simu(1:end-1),'b','LineWidth',2);
title('average \phi_t');xlabel('t');
subplot(2,2,2);
plot(1:periods_simulations,zeta_simu(1:end-1),'b','LineWidth',2);
title('average \zeta_t');xlabel('t')
subplot(2,2,3);
plot(1:periods_simulations,exp_planner_disc_util_simu,'b','LineWidth',2);
title('J(y,\phi,\zeta): average planner value');xlabel('t')
subplot(2,2,4);
plot(1:periods_simulations-1,bonds_simu(1:end-1),'b','LineWidth',2);%
title('average asset holdings');xlabel('t')


%%%%%%%%%%%%%%%%%%%%%%%
%    Figure 9: simulated series, avg. allocations/2   %
%%%%%%%%%%%%%%%%%%%%%%%

figure(9);
subplot(4,2,[1 3]);
plot(1:periods_simulations,c1_simu,'b','LineWidth',2);
title('average consumption'); xlabel('t')
subplot(4,2,[2 4]);
plot(1:periods_simulations,a1_simu,'b','LineWidth',2);
title('average effort');xlabel('t')
subplot(4,2,[5 7]);
plot(1:periods_simulations,exp_disc_util_simu,'b','LineWidth',2);
title('U(y,\phi,\zeta): average lifetime utility (agent)');xlabel('t')
subplot(4,2,6);
plot(1:periods_simulations,lambda1_simu,'b','LineWidth',2);
title('average \lambda ');xlabel('t');
subplot(4,2,8);
plot(1:periods_simulations,eta_simu,'r','LineWidth',2);
title('average  \eta');xlabel('t');
