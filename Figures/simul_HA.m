% simulations
disp(' ');disp(' ');
disp(sprintf(' simulate series          '));
disp(' ');
rand('state',0);


% allocate memory
phi1_simuNext_s1 = zeros(1, periods_simulations);
phi1_simuNext_s2 = zeros(1, periods_simulations);
phi1_simuNext =   zeros(1, periods_simulations);
phi2_simuNext_s1 = zeros(1, periods_simulations);
phi2_simuNext_s2 = zeros(1, periods_simulations);
phi2_simuNext =   zeros(1, periods_simulations);
phi1_simu = zeros(1, periods_simulations+1);
phi1_simu(:,1) = gam1;
a1_simu = zeros(1, periods_simulations);
eta1_simu = zeros(1, periods_simulations);
c_simu =zeros(1, periods_simulations);
a_simu = zeros(1, periods_simulations);
lambda_simu = zeros(1, periods_simulations);
bonds_simu = zeros(1, periods_simulations);
eta_simu  = zeros(1, periods_simulations);
zeta_simu  = zeros(1, periods_simulations+1);
exp_disc_util_simu =  zeros(1, periods_simulations);
exp_planner_disc_util_simu =  zeros(1, periods_simulations);
lambda1_simu= zeros(1, periods_simulations);
prob1_s1_simu= zeros(1, periods_simulations);
prob1_s2_simu= zeros(1, periods_simulations);
prob2_s1_simu= zeros(1, periods_simulations);
prob2_s2_simu= zeros(1, periods_simulations);
transfer_simu= zeros(1, periods_simulations);
g1_simu = zeros(1, periods_simulations);

phi1_simuNext_s1_real = zeros(random_generations, periods_simulations);
phi1_simuNext_s2_real = zeros(random_generations, periods_simulations);
phi1_simuNext_real =   zeros(random_generations, periods_simulations);
phi2_simuNext_s1_real = zeros(random_generations, periods_simulations);
phi2_simuNext_s2_real = zeros(random_generations, periods_simulations);
phi2_simuNext_real =   zeros(random_generations, periods_simulations);
phi1_simu_real = zeros(random_generations, periods_simulations+1);
phi1_simu_real(:,1) = phizero;
a1_simu_real = zeros(random_generations, periods_simulations);
a2_simu_real = zeros(random_generations, periods_simulations);
a_simu_real = zeros(random_generations, periods_simulations);
eta1_simu_real = zeros(random_generations, periods_simulations);
eta2_simu_real = zeros(random_generations, periods_simulations);
eta_simu_real = zeros(random_generations, periods_simulations);
zeta_simu_real = zeros(random_generations, periods_simulations+1);
zeta_simu_real(:,1) = 0;
lambda1_simu_real= zeros(random_generations, periods_simulations);
lambda2_simu_real= zeros(random_generations, periods_simulations);
lambda_simu_real= zeros(random_generations, periods_simulations);
prob1_s1_simu_real= zeros(random_generations, periods_simulations);
prob1_s2_simu_real= zeros(random_generations, periods_simulations);
prob2_s1_simu_real= zeros(random_generations, periods_simulations);
prob2_s2_simu_real= zeros(random_generations, periods_simulations);
transfer_simu_real= zeros(random_generations, periods_simulations);
g1_simu_real = zeros(random_generations, periods_simulations);
bonds_simu_real = zeros(random_generations, periods_simulations);
bonds_simu_real(:,1) =  - betta.*( funeval(parexp_planner1,fspace,[phizero 0])  ...
    - phizero*funeval(parexp1,fspace,[phizero 0]));%
exp_disc_util_simu_real =  zeros(random_generations, periods_simulations);
exp_planner_disc_util_simu_real =  zeros(random_generations, periods_simulations);



for ipet=1:metit
    g1_simu_old = g1_simu(1,:);
    phi1_simu_old = phi1_simu(1,:);
    a1_simu_old = a1_simu(1,:);
    lambda1_simu_old= lambda1_simu(1,:);
    c_simu_old = c_simu;
    a_simu_old= a_simu;
    lambda_simu_old = lambda_simu;
    eta_simu_old = eta_simu;
    zeta_simu_old= zeta_simu;
    bonds_simu_old = bonds_simu;
    transfer_simu_old = transfer_simu;
    exp_disc_util_simu_old =      exp_disc_util_simu;
    exp_planner_disc_util_simu_old =      exp_planner_disc_util_simu;

    state_simu1 = zeros(random_generations, periods_simulations);%sparse(zeros(random_generations, periods_simulations));
    state_ss=zeros(random_generations,1);

    g1_simu_real = zeros(random_generations, periods_simulations);
    g1_simu_real(:,1)=s(1);

    % generate series
    for i =1:periods_simulations;

        c1_simu_real(:,i)  = funeval(parc1,fspace, [phi1_simu_real(:,i)   zeta_simu_real(:,i) ] );
        eta_simu_real(:,i)  =  zeta_simu_real(:,i)./betta + (  1 - ...
            phi1_simu_real(:,i).*((c1_simu_real(:,i)).^(-sig)))./((-sig).*((c1_simu_real(:,i).^(-sig-1))))   ;
        a1_simu_real(:,i)  = funeval(para1,fspace,  [phi1_simu_real(:,i)   zeta_simu_real(:,i) ]  );
        prob1_s1_simu_real(:,i) = a1_simu_real(:,i).^nu;%1 - exp(-rho.*a1_simu(:,i));
        prob1_s2_simu_real(:,i) = 1 - prob1_s1_simu_real(:,i);
        lambda1_simu_real(:,i)  = funeval(parlambda1,fspace, [phi1_simu_real(:,i)   zeta_simu_real(:,i) ] );

        exp_disc_util_simu_real(:,i) = funeval(parexp1,fspace,  [phi1_simu_real(:,i)   zeta_simu_real(:,i) ] ) ;
        exp_planner_disc_util_simu_real(:,i) = state_simu1(:,i).*(s(2) - s(1) ) + funeval(parexp_planner1,fspace,  [phi1_simu_real(:,i)   zeta_simu_real(:,i) ] ) ;

        if i<=periods_simulations-1

            state_ss = rand(random_generations,1);
            [strow1] = find(state_ss>prob1_s1_simu_real(:,i));
            state_simu1(strow1,i+1) =1;
            g1_simu_real(:,i+1) = (1- state_simu1(:,i+1))*s(1) +state_simu1(:,i+1)*s(2);


            phi1_simuNext_s1_real(:,i) =  phi1_simu_real(:,i) ...
                + lambda1_simu_real(:,i).*(nu./a1_simu_real(:,i));
            phi1_simuNext_s2_real(:,i) = phi1_simu_real(:,i) ...
                + lambda1_simu_real(:,i).*(-nu.*(a1_simu_real(:,i).^(nu-1)))./prob1_s2_simu_real(:,i);
            phi1_simuNext_real(:,i) = (1-state_simu1(:,i+1)).*phi1_simuNext_s1_real(:,i) + ...
                state_simu1(:,i+1).* phi1_simuNext_s2_real(:,i) ;
            phi1_simu_real(:,i+1) = phi1_simuNext_real(:,i);

            zeta_simu_real(:,i+1) = eta_simu_real(:,i);
            
            bonds_simu_real(:,i+1) =-betta*(prob1_s1_simu_real(:,i).*(...
                (zeta_simu_real(:,i+1)./phi1_simuNext_s1_real(:,i) )./betta + ...
                funeval(parexp_planner1,fspace,  [phi1_simuNext_s1_real(:,i)   zeta_simu_real(:,i+1)  ] ) - ...
                phi1_simuNext_s1_real(:,i).* funeval(parexp1,fspace,   [phi1_simuNext_s1_real(:,i)   zeta_simu_real(:,i+1) ] )) +...
                prob1_s2_simu_real(:,i).*(...
                (zeta_simu_real(:,i+1)./phi1_simuNext_s2_real(:,i))./betta  + ...
                (s(1) -s(2) +funeval(parexp_planner1,fspace,  [phi1_simuNext_s2_real(:,i)   zeta_simu_real(:,i+1)  ] )) - ...
                phi1_simuNext_s2_real(:,i).* funeval(parexp1,fspace,  [phi1_simuNext_s2_real(:,i)   zeta_simu_real(:,i+1)  ] ) ));
        
        end;
    end  

    c_simu_real = c1_simu_real;
    transfer_simu_real  =   - 1 + c1_simu_real  + g1_simu_real;
    exp_disc_util_simu =     ((ipet-1).*exp_disc_util_simu_old + mean(exp_disc_util_simu_real ,1))./ipet;
    exp_planner_disc_util_simu =     ((ipet-1).*exp_planner_disc_util_simu_old + ...
        mean(exp_planner_disc_util_simu_real ,1))./ipet;
    c_simu =     ((ipet-1).*c_simu_old + mean(c1_simu_real ,1))./ipet;
    a1_simu=    ((ipet-1).*a1_simu_old + mean(a1_simu_real ,1))./ipet;
    prob_s1_simu = a1_simu.^nu;
    prob_s2_simu = 1 - prob_s1_simu;
    lambda1_simu = ((ipet-1).*lambda1_simu_old + mean(lambda1_simu_real ,1))./ipet;
    transfer_simu = ((ipet-1).*transfer_simu_old + mean(transfer_simu_real ,1))./ipet;

    eta_simu =   ((ipet-1).*eta_simu_old + mean(eta_simu_real ,1))./ipet;
    zeta_simu =   ((ipet-1).*zeta_simu_old + mean(zeta_simu_real ,1))./ipet;

    bonds_simu =       ((ipet-1).*bonds_simu_old + mean(bonds_simu_real ,1))./ipet;

    phi1_simu =((ipet-1).*phi1_simu_old + mean(phi1_simu_real ,1))./ipet;

    c1_simu = c_simu;
    g1_simu =((ipet-1).*g1_simu_old + mean(g1_simu_real ,1))./ipet;

end
