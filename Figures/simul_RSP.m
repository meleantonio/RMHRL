% simulations for focs
% tic;


% disp(sprintf(' simulate series          '));
rand('state',0);

cum_change = zeros(random_generations, periods_simulations);

phi1_simuNext_s1s1 = zeros(1, periods_simulations);
phi1_simuNext_s1s2 = zeros(1, periods_simulations);
phi1_simuNext_s2s2 = zeros(1, periods_simulations);
phi1_simuNext_s2s1 = zeros(1, periods_simulations);
phi1_simuNext =   zeros(1, periods_simulations);

phi1_simu = zeros(1, periods_simulations+1);
phi1_simu(:,1) = omega2/omega1;
% phi2_simu = zeros(random_generations, 101);
a1_simu = zeros(1, periods_simulations);
a2_simu = zeros(1, periods_simulations);

c1_simu =zeros(1, periods_simulations);
c2_simu =zeros(1, periods_simulations);

y1_simu =zeros(1, periods_simulations);
y2_simu =zeros(1, periods_simulations);

k1next_simu =zeros(1, periods_simulations+1);
k2next_simu =zeros(1, periods_simulations+1);
% a_simu = zeros(1, periods_simulations);
% lambda_simu = zeros(1, periods_simulations);
% bonds_simu = zeros(1, periods_simulations+1);

invest1_simu =zeros(1, periods_simulations);
invest2_simu =zeros(1, periods_simulations);

nx1_simu =zeros(1, periods_simulations);
nx2_simu =zeros(1, periods_simulations);



varphi1_simu = zeros(1, periods_simulations+1);
% phi2_simu = zeros(random_generations, 101);
vara1_simu = zeros(1, periods_simulations);
vara2_simu = zeros(1, periods_simulations);

varc1_simu =zeros(1, periods_simulations);
varc2_simu =zeros(1, periods_simulations);

vark1next_simu =zeros(1, periods_simulations+1);
vark2next_simu =zeros(1, periods_simulations+1);

varlambda1_simu= zeros(1, periods_simulations);
varlambda2_simu= zeros(1, periods_simulations);


exp_disc_util1_simu =  zeros(1, periods_simulations);
exp_disc_util2_simu =  zeros(1, periods_simulations);
exp_planner_disc_util_simu =  zeros(1, periods_simulations);


% a2_simu= zeros(random_generations, 100);
% eta_simu= zeros(random_generations, periods_simulations);
lambda1_simu= zeros(1, periods_simulations);
lambda2_simu= zeros(1, periods_simulations);

% value_simu = zeros(random_generations,periods_simulations);
% lambda2_simu= zeros(random_generations, 100);
prob1_s1_simu= zeros(1, periods_simulations);
prob1_s2_simu= zeros(1, periods_simulations);
prob2_s1_simu= zeros(1, periods_simulations);
prob2_s2_simu= zeros(1, periods_simulations);

%     transfer_simu= zeros(1, periods_simulations);
g1_simu = zeros(1, periods_simulations);
g2_simu = zeros(1, periods_simulations);



phi1_simuNext_s1s1_real = zeros(random_generations, periods_simulations);
phi1_simuNext_s1s2_real = zeros(random_generations, periods_simulations);
phi1_simuNext_s2s1_real = zeros(random_generations, periods_simulations);
phi1_simuNext_s2s2_real = zeros(random_generations, periods_simulations);


phi1_simu_real = zeros(random_generations, periods_simulations+1);
phi1_simu_real(:,1) = omega2/omega1;
% phi2_simu = zeros(random_generations, 101);
a1_simu_real = zeros(random_generations, periods_simulations);
a2_simu_real = zeros(random_generations, periods_simulations);
%                 a_simu_real = zeros(random_generations, periods_simulations);
k1next_simu_real = zeros(random_generations, periods_simulations+1);
k2next_simu_real = zeros(random_generations, periods_simulations+1);
k1next_simu_real(:,1) = 3.1;%k_min;%.1;%3.1;%.2125;%2;%.5*k_ss;
k2next_simu_real(:,1) = 3.1;%k_min;%.1;%3.1;%.2125;%3;%.5*k_ss;



dprob1_s1_simu_real= zeros(random_generations, periods_simulations);
dprob1_s2_simu_real= zeros(random_generations, periods_simulations);

dprob2_s1_simu_real= zeros(random_generations, periods_simulations);
dprob2_s2_simu_real= zeros(random_generations, periods_simulations);



invest1_simu_real = zeros(random_generations, periods_simulations);
invest2_simu_real = zeros(random_generations, periods_simulations);

c1_simu_real= zeros(random_generations, periods_simulations);
c2_simu_real= zeros(random_generations, periods_simulations);
y1_simu_real= zeros(random_generations, periods_simulations);
y2_simu_real= zeros(random_generations, periods_simulations);

nx1_simu_real= zeros(random_generations, periods_simulations);
nx2_simu_real= zeros(random_generations, periods_simulations);


% a2_simu= zeros(random_generations, 100);
% eta_simu= zeros(random_generations, periods_simulations);
lambda1_simu_real= zeros(random_generations, periods_simulations);
lambda2_simu_real= zeros(random_generations, periods_simulations);
%         lambda_simu_real= zeros(random_generations, periods_simulations);

% value_simu = zeros(random_generations,periods_simulations);
% lambda2_simu= zeros(random_generations, 100);
prob1_s1_simu_real= zeros(random_generations, periods_simulations);
prob1_s2_simu_real= zeros(random_generations, periods_simulations);

prob2_s1_simu_real= zeros(random_generations, periods_simulations);
prob2_s2_simu_real= zeros(random_generations, periods_simulations);

%     transfer_simu_real= zeros(random_generations, periods_simulations);
g1_simu_real = zeros(random_generations, periods_simulations);
g2_simu_real = zeros(random_generations, periods_simulations);

%     bonds_simu_real = zeros(random_generations, periods_simulations+1);
%     bonds_simu_real(:,1) = -betta*funeval(parexp_planner1,fspace,gam1_2(1,1)) + gam1_2(1,1)*funeval(parexp1,fspace,gam1_2(1,1));%


exp_disc_util1_simu_real =  zeros(random_generations, periods_simulations);
exp_disc_util2_simu_real =  zeros(random_generations, periods_simulations);
exp_planner_disc_util_simu_real =  zeros(random_generations, periods_simulations);

exp_transfer= 0;
exp_transfer2 = 0 ;
exp_value_planner_simu = 0 ;
exp_value_planner_simu2 = 0 ;
exp_value_agent_simu = 0 ;
exp_value_agent_simu2 =0;


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
    k1next_simu_old = k1next_simu;
    k2next_simu_old = k2next_simu;

        invest1_simu_old = invest1_simu;
    invest2_simu_old = invest2_simu;
        nx1_simu_old = nx1_simu;
    nx2_simu_old = nx2_simu;
    y1_simu_old = c1_simu;
    y2_simu_old = c2_simu;


varphi1_simu_old = varphi1_simu(1,:);
    vara1_simu_old = vara1_simu(1,:);
    varlambda1_simu_old= varlambda1_simu(1,:);
    vara2_simu_old = vara2_simu(1,:);
    varlambda2_simu_old= varlambda2_simu(1,:);
    varc1_simu_old = varc1_simu;
    varc2_simu_old = varc2_simu;
    vark1next_simu_old = vark1next_simu;
    vark2next_simu_old = vark2next_simu;

    % a_simu_old= a_simu;
    % lambda_simu_old = lambda_simu;
    %
    % bonds_simu_old = bonds_simu;
    % transfer_simu_old = transfer_simu;

    exp_disc_util1_simu_old =      exp_disc_util1_simu;
    exp_disc_util2_simu_old =      exp_disc_util2_simu;
    exp_planner_disc_util_simu_old =      exp_planner_disc_util_simu;

    exp_transfer_old =   exp_transfer;
    exp_transfer2_old =       exp_transfer2 ;
    exp_value_planner_simu_old=   exp_value_planner_simu ;
    exp_value_planner_simu2_old =    exp_value_planner_simu2 ;
    exp_value_agent_simu_old =     exp_value_agent_simu ;
    exp_value_agent_simu2_old =     exp_value_agent_simu2 ;


    state_simu1 = zeros(random_generations, periods_simulations);%sparse(zeros(random_generations, periods_simulations));
    state_simu2 = zeros(random_generations, periods_simulations);%sparse(zeros(random_generations, periods_simulations));
    state_ss1=zeros(random_generations,1);
    state_ss2=zeros(random_generations,1);

    %        statess1=rand(random_generations,periods_simulations);
    %
    %         [sssrow1,ssscol1] = find(statess1>.5);
    %          sssval1 = ones(length(sssrow1),1);
    %         state_simu1 = (sparse(sssrow1,ssscol1,sssval1,random_generations,periods_simulations));
    g1_simu_real = zeros(random_generations, periods_simulations);
    %          g1_simu = (1-state_simu1)*s(1) + state_simu1*s(2);
    g1_simu_real(:,1)=s1(1);
    g2_simu_real = zeros(random_generations, periods_simulations);
    %          g1_simu = (1-state_simu1)*s(1) + state_simu1*s(2);
    g2_simu_real(:,1)=s2(2);%s2(1);%


    % for j=2:1001
    for i =1:periods_simulations;% 99;% 100
        
       gridd = [phi1_simu_real(:,i)  k1next_simu_real(:,i) k2next_simu_real(:,i)];
        
       k1next_simu_real(:,i+1)  = funeval(park1,fspace, gridd);
       k2next_simu_real(:,i+1)  = funeval(park2,fspace, gridd);

        invest1_simu_real(:,i) = k1next_simu_real(:,i+1) - (1-delta1).*k1next_simu_real(:,i);
       invest2_simu_real(:,i) = k2next_simu_real(:,i+1) - (1-delta2).*k2next_simu_real(:,i);
        
        a1_simu_real(:,i)  = funeval(para1,fspace, gridd  );

        prob1_s1_simu_real(:,i) = a1_simu_real(:,i).^nu;%1 - exp(-rho.*a1_simu(:,i));
        prob1_s2_simu_real(:,i) = 1 - prob1_s1_simu_real(:,i);

        a2_simu_real(:,i)  = funeval(para2,fspace,     gridd);

        prob2_s1_simu_real(:,i) = a2_simu_real(:,i).^nu2;
        prob2_s2_simu_real(:,i) = 1 - prob2_s1_simu_real(:,i);

        dprob1_s1_simu_real(:,i) = nu.*(a1_simu_real(:,i).^(nu-1));%1 - exp(-rho.*a1_simu(:,i));
        dprob1_s2_simu_real(:,i) =  - dprob1_s1_simu_real(:,i);

        dprob2_s1_simu_real(:,i) = nu.*(a2_simu_real(:,i).^(nu-1));%1 - exp(-rho.*a1_simu(:,i));
        dprob2_s2_simu_real(:,i) =  - dprob2_s1_simu_real(:,i);


        lambda1_simu_real(:,i)  = funeval(parlambda1,fspace,  gridd );
        lambda2_simu_real(:,i)  = funeval(parlambda2,fspace,  gridd );
        %              [phi1_simu(:,i)  phi2_simu(:,i)  g1_simu(:,i)   g2_simu(:,i)]);


    

        c1_simu_real(:,i)  = (- k1next_simu_real(:,i+1) - k2next_simu_real(:,i+1) ...
            +  (1-delta1).*gridd(:,2) + (1-delta2).*gridd(:,3)  +  g1_simu_real(:,i).*(gridd(:,2).^elast1)...
            + g2_simu_real(:,i).*(gridd(:,3).^elast2) )./(1 + phi1_simu_real(:,i).^(1/sig)); %
         c2_simu_real(:,i)  = (- k1next_simu_real(:,i+1) - k2next_simu_real(:,i+1) ...
            +  (1-delta1).*gridd(:,2) + (1-delta2).*gridd(:,3)  +  g1_simu_real(:,i).*(gridd(:,2).^elast1)...
            + g2_simu_real(:,i).*(gridd(:,3).^elast2) ) - c1_simu_real(:,i); %
        %         bonds_simu_real(:,i+1) =-(1-g1_simu_real(:,i) - c1_simu_real(:,i)) +(betta).*bonds_simu_real(:,i);% funeval(parplanner1,fspace, phi1_simu_real(:,i)) - gam1_2(1,1)*V_outside_option;%

        
        y1_simu_real(:,i) = g1_simu_real(:,i).*(gridd(:,2).^elast1); 
        y2_simu_real(:,i) = g2_simu_real(:,i).*(gridd(:,3).^elast2);

        
        
        exp_disc_util1_simu_real(:,i) = funeval(parexp1,fspace, gridd  ) ;
        exp_disc_util2_simu_real(:,i) = funeval(parexp2,fspace, gridd );

        exp_planner_disc_util_simu_real(:,i) = funeval(parexp_planner1,fspace, gridd );
        
        nx1_simu_real(:,i) = (y1_simu_real(:,i)  - c1_simu_real(:,i) - invest1_simu_real(:,i))./(y1_simu_real(:,i) );
        nx2_simu_real(:,i) = (y2_simu_real(:,i)  - c2_simu_real(:,i) - invest2_simu_real(:,i))./(y2_simu_real(:,i) );


        if i<=periods_simulations-1

            state_ss1 = rand(random_generations,1);
            state_ss2 = rand(random_generations,1);

            %             [strow1] = find(state_ss>((1- state_simu1(:,i)).*prob_s1_simu_real(:,i)+...
            %                 state_simu1(:,i).*prob_s2_simu_real(:,i)));
            [strow1] = find(state_ss1>prob1_s1_simu_real(:,i));
            state_simu1(strow1,i+1) =1;
            g1_simu_real(:,i+1) = (1- state_simu1(:,i+1))*s1(1) +state_simu1(:,i+1)*s1(2); % s1(2);%

            [strow2] = find(state_ss2>prob2_s1_simu_real(:,i));
            state_simu2(strow2,i+1) =1;
            g2_simu_real(:,i+1) =(1- state_simu2(:,i+1))*s2(1) +state_simu2(:,i+1)*s2(2); % s2(2);%


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

    end  % clear sssval;



  
    exp_disc_util1_simu =     ((ipet-1).*exp_disc_util1_simu_old + mean(exp_disc_util1_simu_real ,1))./ipet;
    exp_disc_util2_simu =     ((ipet-1).*exp_disc_util2_simu_old + mean(exp_disc_util2_simu_real ,1))./ipet;
    exp_planner_disc_util_simu =     ((ipet-1).*exp_planner_disc_util_simu_old + ...
        mean(exp_planner_disc_util_simu_real ,1))./ipet;



    c1_simu =     ((ipet-1).*c1_simu_old + mean(c1_simu_real ,1))./ipet;
    c2_simu =     ((ipet-1).*c2_simu_old + mean(c2_simu_real ,1))./ipet;
    a1_simu =     ((ipet-1).*a1_simu_old + mean(a1_simu_real ,1))./ipet;
    a2_simu =     ((ipet-1).*a2_simu_old + mean(a2_simu_real ,1))./ipet;
    k1next_simu =     ((ipet-1).*k1next_simu_old + mean(k1next_simu_real ,1))./ipet;
    k2next_simu =     ((ipet-1).*k2next_simu_old + mean(k2next_simu_real ,1))./ipet;

   lambda1_simu =     ((ipet-1).*lambda1_simu_old + mean(lambda1_simu_real ,1))./ipet;
    lambda2_simu =     ((ipet-1).*lambda2_simu_old + mean(lambda2_simu_real ,1))./ipet;
    invest1_simu =     ((ipet-1).*invest1_simu_old + mean(invest1_simu_real ,1))./ipet;
    invest2_simu =     ((ipet-1).*invest2_simu_old + mean(invest2_simu_real ,1))./ipet;

    nx1_simu =     ((ipet-1).*nx1_simu_old + mean(nx1_simu_real ,1))./ipet;
    nx2_simu =     ((ipet-1).*nx2_simu_old + mean(nx2_simu_real ,1))./ipet;
    y1_simu =     ((ipet-1).*y1_simu_old + mean(y1_simu_real ,1))./ipet;
    y2_simu =     ((ipet-1).*y2_simu_old + mean(y2_simu_real ,1))./ipet;

    phi1_simu =((ipet-1).*phi1_simu_old + mean(phi1_simu_real ,1))./ipet;

    g1_simu =((ipet-1).*g1_simu_old + mean(g1_simu_real ,1))./ipet;
g2_simu =((ipet-1).*g2_simu_old + mean(g2_simu_real ,1))./ipet;


varc1_simu = ((ipet-1).*varc1_simu_old + var(c1_simu_real ,0,1))./ipet;

    varc2_simu =     ((ipet-1).*varc2_simu_old + var(c2_simu_real ,0,1))./ipet;
    vara1_simu =     ((ipet-1).*vara1_simu_old + var(a1_simu_real ,0,1))./ipet;
    vara2_simu =     ((ipet-1).*vara2_simu_old + var(a2_simu_real ,0,1))./ipet;
    vark1next_simu =     ((ipet-1).*vark1next_simu_old + var(k1next_simu_real ,0,1))./ipet;
    vark2next_simu =     ((ipet-1).*vark2next_simu_old + var(k2next_simu_real ,0,1))./ipet;

   varlambda1_simu =     ((ipet-1).*varlambda1_simu_old + var(lambda1_simu_real ,0,1))./ipet;
    varlambda2_simu =     ((ipet-1).*varlambda2_simu_old + var(lambda2_simu_real ,0,1))./ipet;


    varphi1_simu =((ipet-1).*varphi1_simu_old + var(phi1_simu_real ,0,1))./ipet;


end
