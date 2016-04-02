% simulations
disp(' ');disp(' ');
disp(sprintf(' simulate series          '));
disp(' ');
rand('state',0);


phi1_simuNext_s1s1 = zeros(1, periods_simulations);
phi1_simuNext_s1s2 = zeros(1, periods_simulations);
phi1_simuNext_s2s2 = zeros(1, periods_simulations);
phi1_simuNext_s2s1 = zeros(1, periods_simulations);
phi1_simuNext =   zeros(1, periods_simulations);

phi1_simu = zeros(1, periods_simulations+1);
phi1_simu(:,1) = omega2/omega1;
a1_simu = zeros(1, periods_simulations);
a2_simu = zeros(1, periods_simulations);

c1_simu =zeros(1, periods_simulations);
c2_simu =zeros(1, periods_simulations);


exp_disc_util1_simu =  zeros(1, periods_simulations);
exp_disc_util2_simu =  zeros(1, periods_simulations);
exp_planner_disc_util_simu =  zeros(1, periods_simulations);

lambda1_simu= zeros(1, periods_simulations);
lambda2_simu= zeros(1, periods_simulations);

prob1_s1_simu= zeros(1, periods_simulations);
prob1_s2_simu= zeros(1, periods_simulations);
prob2_s1_simu= zeros(1, periods_simulations);
prob2_s2_simu= zeros(1, periods_simulations);

g1_simu = zeros(1, periods_simulations);
g2_simu = zeros(1, periods_simulations);

phi1_simuNext_s1s1_real = zeros(random_generations, periods_simulations);
phi1_simuNext_s1s2_real = zeros(random_generations, periods_simulations);
phi1_simuNext_s2s1_real = zeros(random_generations, periods_simulations);
phi1_simuNext_s2s2_real = zeros(random_generations, periods_simulations);


phi1_simu_real = zeros(random_generations, periods_simulations+1);
phi1_simu_real(:,1) = omega2/omega1;
a1_simu_real = zeros(random_generations, periods_simulations);
a2_simu_real = zeros(random_generations, periods_simulations);
lambda1_simu_real= zeros(random_generations, periods_simulations);
lambda2_simu_real= zeros(random_generations, periods_simulations);

prob1_s1_simu_real= zeros(random_generations, periods_simulations);
prob1_s2_simu_real= zeros(random_generations, periods_simulations);

prob2_s1_simu_real= zeros(random_generations, periods_simulations);
prob2_s2_simu_real= zeros(random_generations, periods_simulations);


g1_simu_real = zeros(random_generations, periods_simulations);
g2_simu_real = zeros(random_generations, periods_simulations);


exp_disc_util1_simu_real =  zeros(random_generations, periods_simulations);
exp_disc_util2_simu_real =  zeros(random_generations, periods_simulations);
exp_planner_disc_util_simu_real =  zeros(random_generations, periods_simulations);


for ipet=1:metit

    g1_simu_old = g1_simu(1,:);
    g2_simu_old = g2_simu(1,:);

    phi1_simu_old = phi1_simu(1,:);
    a1_simu_old = a1_simu(1,:);
    lambda1_simu_old= lambda1_simu(1,:);
    a2_simu_old = a2_simu(1,:);
    lambda2_simu_old= lambda2_simu(1,:);
    c1_simu_old = c1_simu;
    c2_simu_old = c2_simu;

    exp_disc_util1_simu_old =      exp_disc_util1_simu;
    exp_disc_util2_simu_old =      exp_disc_util2_simu;
    exp_planner_disc_util_simu_old =      exp_planner_disc_util_simu;


    state_simu1 = zeros(random_generations, periods_simulations);%sparse(zeros(random_generations, periods_simulations));
    state_simu2 = zeros(random_generations, periods_simulations);%sparse(zeros(random_generations, periods_simulations));
    state_ss1=zeros(random_generations,1);
    state_ss2=zeros(random_generations,1);

    g1_simu_real = zeros(random_generations, periods_simulations);
    g1_simu_real(:,1)=s1(1);
    g2_simu_real = zeros(random_generations, periods_simulations);
    g2_simu_real(:,1)=s2(1);


    for i =1:periods_simulations;

        a1_simu_real(:,i)  = funeval(para1,fspace, phi1_simu_real(:,i)   );
        prob1_s1_simu_real(:,i) = a1_simu_real(:,i).^nu;%1 - exp(-rho.*a1_simu(:,i));
        prob1_s2_simu_real(:,i) = 1 - prob1_s1_simu_real(:,i);

        a2_simu_real(:,i)  = funeval(para2,fspace, phi1_simu_real(:,i)   );
        prob2_s1_simu_real(:,i) = a2_simu_real(:,i).^nu2;
        prob2_s2_simu_real(:,i) = 1 - prob2_s1_simu_real(:,i);

        dprob1_s1_simu_real(:,i) = nu.*(a1_simu_real(:,i).^(nu-1));
        dprob1_s2_simu_real(:,i) =  - dprob1_s1_simu_real(:,i);

        dprob2_s1_simu_real(:,i) = nu2.*(a2_simu_real(:,i).^(nu2-1));
        dprob2_s2_simu_real(:,i) =  - dprob2_s1_simu_real(:,i);


        lambda1_simu_real(:,i)  = funeval(parlambda1,fspace,phi1_simu_real(:,i) );
        lambda2_simu_real(:,i)  = funeval(parlambda2,fspace,phi1_simu_real(:,i) );

        c1_simu_real(:,i)  = (2*output_y - g1_simu_real(:,i) - g2_simu_real(:,i))./(1 + ...
            phi1_simu_real(:,i).^(1/sig)); %
        c2_simu_real(:,i)  = (2*output_y - g1_simu_real(:,i) - g2_simu_real(:,i)) - c1_simu_real(:,i); %

        exp_disc_util1_simu_real(:,i) = (c1_simu_real(:,i).^(1-sig))./(1-sig) -  ...
            (((2*output_y - s1(1) - s2(1))./(1 + ...
            phi1_simu_real(:,i).^(1/sig))).^(1-sig))./(1-sig)  + funeval(parexp1,fspace,phi1_simu_real(:,i) ) ;
        exp_disc_util2_simu_real(:,i) = (c2_simu_real(:,i).^(1-sig))./(1-sig) ...
            -  ( (2*output_y - s1(1) - s2(1) - ...
            (2*output_y - s1(1) - s2(1))./(1 + ...
            phi1_simu_real(:,i).^(1/sig))).^(1-sig))./(1-sig)   +funeval(parexp2,fspace,phi1_simu_real(:,i) );

        exp_planner_disc_util_simu_real(:,i) = (c1_simu_real(:,i).^(1-sig))./(1-sig) -  ...
            (((2*output_y - s1(1) - s2(1))./(1 + ...
            phi1_simu_real(:,i).^(1/sig))).^(1-sig))./(1-sig)  + ...
            phi1_simu_real(:,i).*((c2_simu_real(:,i).^(1-sig))./(1-sig) ...
            -  ( (2*output_y - s1(1) - s2(1) - ...
            (2*output_y - s1(1) - s2(1))./(1 + ...
            phi1_simu_real(:,i).^(1/sig))).^(1-sig))./(1-sig)  )   +...
            funeval(parexp_planner1,fspace,phi1_simu_real(:,i));

        if i<=periods_simulations-1

            state_ss1 = rand(random_generations,1);
            state_ss2 = rand(random_generations,1);

            [strow1] = find(state_ss1>prob1_s1_simu_real(:,i));
            state_simu1(strow1,i+1) =1;
            g1_simu_real(:,i+1) = (1- state_simu1(:,i+1))*s1(1) +state_simu1(:,i+1)*s1(2);

            [strow2] = find(state_ss2>prob2_s1_simu_real(:,i));
            state_simu2(strow2,i+1) =1;
            g2_simu_real(:,i+1) = (1- state_simu2(:,i+1))*s2(1) +state_simu2(:,i+1)*s2(2);


            phi1_simuNext_s1s1_real(:,i) =  phi1_simu_real(:,i).*(1+...
                lambda2_simu_real(:,i).*dprob2_s1_simu_real(:,i)./prob2_s1_simu_real(:,i) )./(1+...
                lambda1_simu_real(:,i).*dprob1_s1_simu_real(:,i)./prob1_s1_simu_real(:,i) ) ;

            phi1_simuNext_s1s2_real(:,i) =  phi1_simu_real(:,i).*(1+...
                lambda2_simu_real(:,i).*dprob2_s2_simu_real(:,i)./prob2_s2_simu_real(:,i) )./(1+...
                lambda1_simu_real(:,i).*dprob1_s1_simu_real(:,i)./prob1_s1_simu_real(:,i) ) ;

            phi1_simuNext_s2s1_real(:,i) =  phi1_simu_real(:,i).*(1+...
                lambda2_simu_real(:,i).*dprob2_s1_simu_real(:,i)./prob2_s1_simu_real(:,i) )./(1+...
                lambda1_simu_real(:,i).*dprob1_s2_simu_real(:,i)./prob1_s2_simu_real(:,i) ) ;

            phi1_simuNext_s2s2_real(:,i) =  phi1_simu_real(:,i).*(1+...
                lambda2_simu_real(:,i).*dprob2_s2_simu_real(:,i)./prob2_s2_simu_real(:,i) )./(1+...
                lambda1_simu_real(:,i).*dprob1_s2_simu_real(:,i)./prob1_s2_simu_real(:,i) ) ;


            phi1_simu_real(:,i+1) = (1- state_simu1(:,i)).*(1- state_simu2(:,i)).*phi1_simuNext_s1s1_real(:,i)...
                + (1- state_simu1(:,i)).*( state_simu2(:,i)).*phi1_simuNext_s1s2_real(:,i)...
                + ( state_simu1(:,i)).*(1- state_simu2(:,i)).*phi1_simuNext_s2s1_real(:,i)...
                +  ( state_simu1(:,i)).*(state_simu2(:,i)).*phi1_simuNext_s2s2_real(:,i);

        end;
    end

    exp_disc_util1_simu =     ((ipet-1).*exp_disc_util1_simu_old + mean(exp_disc_util1_simu_real ,1))./ipet;
    exp_disc_util2_simu =     ((ipet-1).*exp_disc_util2_simu_old + mean(exp_disc_util2_simu_real ,1))./ipet;
    exp_planner_disc_util_simu =     ((ipet-1).*exp_planner_disc_util_simu_old + ...
        mean(exp_planner_disc_util_simu_real ,1))./ipet;

    c1_simu =     ((ipet-1).*c1_simu_old + mean(c1_simu_real ,1))./ipet;
    c2_simu =     ((ipet-1).*c2_simu_old + mean(c2_simu_real ,1))./ipet;
    a1_simu =     ((ipet-1).*a1_simu_old + mean(a1_simu_real ,1))./ipet;
    a2_simu =     ((ipet-1).*a2_simu_old + mean(a2_simu_real ,1))./ipet;
    lambda1_simu =     ((ipet-1).*lambda1_simu_old + mean(lambda1_simu_real ,1))./ipet;
    lambda2_simu =     ((ipet-1).*lambda2_simu_old + mean(lambda2_simu_real ,1))./ipet;

    phi1_simu =((ipet-1).*phi1_simu_old + mean(phi1_simu_real ,1))./ipet;

    g1_simu =((ipet-1).*g1_simu_old + mean(g1_simu_real ,1))./ipet;
    g2_simu =((ipet-1).*g2_simu_old + mean(g2_simu_real ,1))./ipet;
end
