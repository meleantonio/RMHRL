function equ = RSE_focs(par00, Gridphi,fspace);


global alpha betta sig  epsil  nu nu2 output_y s1 s2


aggr_output = 2*output_y - [s1(2)+s2(2); s1(1)+s2(2); s1(2)+s2(1); s1(1)+s2(1)];


par0 = reshape(par00,length(par00)/7,7 );
%                     par0 = reshape(par00(1:end),length(par00)/5,5 );

parlambda1 =               par0(:,1);
para1 =                         par0(:,2);
parexp_planner1 =      par0(:,7);
%                     parexp_planner2 =      par0(:,4);
parexp1 =                     par0(:,3);
parexp2 =                     par0(:,4);
parlambda2 =               par0(:,5);
para2 =                     par0(:,6);


a1 =                     funeval(para1,fspace, Gridphi);
a2 =                     funeval(para2,fspace, Gridphi);
lambda1 = funeval(parlambda1,fspace, Gridphi);
lambda2 = funeval(parlambda2,fspace, Gridphi);


phi1 = Gridphi(:,1);




% consumption and leisure when phi's are not zero
aggr_y   = aggr_output(4);
c1 = aggr_y./(  1+ phi1.^(1/sig));
c2 = aggr_y - c1;

c1_s1s2 = aggr_output(2)./(  1+ phi1.^(1/sig)); %;
c1_s2s1 = aggr_output(3)./(  1+ phi1.^(1/sig)); %;
c1_s2s2 = aggr_output(1)./(  1+ phi1.^(1/sig)); %;

c2_s1s2 = aggr_output(2) - c1_s1s2 ;
c2_s2s1 = aggr_output(3) - c1_s2s1 ;
c2_s2s2 = aggr_output(1) - c1_s2s2 ;


% probabilities
prob1_s1 =  a1.^nu;
prob1_s2 = 1 - prob1_s1;

dprob1_s1 =  nu.*(a1.^(nu-1));
dprob1_s2 =  - dprob1_s1;


% probabilities

prob2_s2 = (1 - a2.^nu2);
prob2_s1 = (a2.^nu2);

dprob2_s2 =  (-nu2.*(a2.^(nu2-1)));
dprob2_s1 =  (nu2.*(a2.^(nu2-1)));




% likelihood ratios
like_ratio1_s1 =  dprob1_s1./prob1_s1;
like_ratio1_s2 = dprob1_s2./prob1_s2;

like_ratio2_s1 =  dprob2_s1./prob2_s1;
like_ratio2_s2 = dprob2_s2./prob2_s2;



% derivatives of likelihood ratios
dlike_ratio1_s1 =  -nu./(a1.^2);%
dlike_ratio1_s2 = (-nu*(nu-1)*(a1.^(nu-2)) - nu*(a1.^(2*nu-2)))./(prob1_s2.^2) ;


dlike_ratio2_s2 =  ((-nu2*(nu2-1)*(a2.^(nu2-2)) - nu2*(a2.^(2*nu2-2)))./(prob2_s2.^2)) ;%
dlike_ratio2_s1 = ( -nu2./(a2.^2))  ;



PHINext_s1s1 = phi1.*(1+lambda2.*like_ratio2_s1)./(1+lambda1.*like_ratio1_s1);
PHINext_s1s2 = phi1.*(1+lambda2.*like_ratio2_s2)./(1+lambda1.*like_ratio1_s1);
PHINext_s2s1 = phi1.*(1+lambda2.*like_ratio2_s1)./(1+lambda1.*like_ratio1_s2);
PHINext_s2s2 = phi1.*(1+lambda2.*like_ratio2_s2)./(1+lambda1.*like_ratio1_s2);



c1Next_s1s1 = aggr_output(4)./(1+PHINext_s1s1.^(1/sig)); %
c1Next_s1s2 = aggr_output(2)./(1+PHINext_s1s2.^(1/sig)); %
c1Next_s2s1 = aggr_output(3)./(1+PHINext_s2s1.^(1/sig)); %
c1Next_s2s2 = aggr_output(1)./(1+PHINext_s2s2.^(1/sig)); %

c2Next_s1s1 = aggr_output(4)-c1Next_s1s1; %
c2Next_s1s2 = aggr_output(2)-c1Next_s1s2; %
c2Next_s2s1 = aggr_output(3)-c1Next_s2s1; %
c2Next_s2s2 = aggr_output(1)-c1Next_s2s2; %



a1Next_s1s1 = funeval(para1,fspace,PHINext_s1s1 );
a1Next_s1s2 = funeval(para1,fspace,PHINext_s1s2 );
a1Next_s2s1 = funeval(para1,fspace,PHINext_s2s1 );
a1Next_s2s2 = funeval(para1,fspace,PHINext_s2s2  );
%
a2Next_s1s1 = funeval(para2,fspace,PHINext_s1s1 );
a2Next_s1s2 = funeval(para2,fspace,PHINext_s1s2 );
a2Next_s2s1 = funeval(para2,fspace,PHINext_s2s1);
a2Next_s2s2 = funeval(para2,fspace,PHINext_s2s2 );


if sig==1
    utilityNext1_s1s1 = log(c1Next_s1s1) - alpha*(a1Next_s1s1.^epsil) ;
    utilityNext1_s2s1 = log(c1Next_s2s1) - alpha*(a1Next_s2s1.^epsil) ;
    utilityNext1_s1s2 = log(c1Next_s1s2) - alpha*(a1Next_s1s2.^epsil) ;
    utilityNext1_s2s2 = log(c1Next_s2s2) - alpha*(a1Next_s2s2.^epsil) ;

    utilityNext2_s1s1 = log(c2Next_s1s1) - alpha*(a2Next_s1s1.^epsil);
    utilityNext2_s2s1 = log(c2Next_s2s1) - alpha*(a2Next_s2s1.^epsil);
    utilityNext2_s1s2 = log(c2Next_s1s2) - alpha*(a2Next_s1s2.^epsil);
    utilityNext2_s2s2 = log(c2Next_s2s2) - alpha*(a2Next_s2s2.^epsil);

else

    utilityNext1_s1s1 = (c1Next_s1s1.^(1-sig))./(1-sig)  - alpha*(a1Next_s1s1.^epsil) ;
    utilityNext1_s2s1 = (c1Next_s2s1.^(1-sig))./(1-sig)  - alpha*(a1Next_s2s1.^epsil) ;
    utilityNext1_s1s2 = (c1Next_s1s2.^(1-sig))./(1-sig)  - alpha*(a1Next_s1s2.^epsil) ;
    utilityNext1_s2s2 = (c1Next_s2s2.^(1-sig))./(1-sig)  - alpha*(a1Next_s2s2.^epsil) ;

    utilityNext2_s1s1 = (c2Next_s1s1.^(1-sig))./(1-sig)  - alpha*(a2Next_s1s1.^epsil);
    utilityNext2_s2s1 = (c2Next_s2s1.^(1-sig))./(1-sig)  - alpha*(a2Next_s2s1.^epsil);
    utilityNext2_s1s2 = (c2Next_s1s2.^(1-sig))./(1-sig)  - alpha*(a2Next_s1s2.^epsil);
    utilityNext2_s2s2 = (c2Next_s2s2.^(1-sig))./(1-sig)  - alpha*(a2Next_s2s2.^epsil);


end;


exp_disc_utility1 = funeval(parexp1,fspace,Gridphi);
exp_disc_utility1_next_s1s1 = funeval(parexp1,fspace,PHINext_s1s1  );
exp_disc_utility1_next_s1s2 = (c1Next_s1s2.^(1-sig))./(1-sig)  - alpha*(a1Next_s1s2.^epsil) ...
    - (c1Next_s1s1.^(1-sig))./(1-sig) +  alpha*(a1Next_s1s1.^epsil) + ...
    funeval(parexp1,fspace,PHINext_s1s2  );
exp_disc_utility1_next_s2s1 = (c1Next_s2s1.^(1-sig))./(1-sig)  - alpha*(a1Next_s2s1.^epsil) ...
    - (c1Next_s1s1.^(1-sig))./(1-sig) +  alpha*(a1Next_s1s1.^epsil) + ...
    funeval(parexp1,fspace,PHINext_s2s1  );
exp_disc_utility1_next_s2s2 =  (c1Next_s2s2.^(1-sig))./(1-sig)  - alpha*(a1Next_s2s2.^epsil) ...
    - (c1Next_s1s1.^(1-sig))./(1-sig) +  alpha*(a1Next_s1s1.^epsil) + ...
    funeval(parexp1,fspace,PHINext_s2s2  );

exp_disc_utility2 = funeval(parexp2,fspace,Gridphi);
exp_disc_utility2_next_s1s1 =  funeval(parexp2,fspace,PHINext_s1s1);
exp_disc_utility2_next_s1s2 = (c2Next_s1s2.^(1-sig))./(1-sig)  - alpha*(a2Next_s1s2.^epsil) ...
    - (c2Next_s1s1.^(1-sig))./(1-sig) + alpha*(a2Next_s1s1.^epsil) + ...
    funeval(parexp2,fspace,PHINext_s1s2  );
exp_disc_utility2_next_s2s1 =  (c2Next_s2s1.^(1-sig))./(1-sig)  - alpha*(a2Next_s2s1.^epsil) ...
    - (c2Next_s1s1.^(1-sig))./(1-sig) + alpha*(a2Next_s1s1.^epsil) + ...
    funeval(parexp2,fspace,PHINext_s2s1 );
exp_disc_utility2_next_s2s2 =  (c2Next_s2s2.^(1-sig))./(1-sig)  - alpha*(a2Next_s2s2.^epsil) ...
    - (c2Next_s1s1.^(1-sig))./(1-sig) + alpha*(a2Next_s1s1.^epsil) + ...
    funeval(parexp2,fspace,PHINext_s2s2  );


exp_planner_disc_utility = funeval(parexp_planner1,fspace,Gridphi);
exp_planner_disc_utility_next_s1s1 = funeval(parexp_planner1,fspace,PHINext_s1s1  );
exp_planner_disc_utility_next_s1s2 =  (c1Next_s1s2.^(1-sig))./(1-sig)  - alpha*(a1Next_s1s2.^epsil) ...
    - (c1Next_s1s1.^(1-sig))./(1-sig) + alpha*(a1Next_s1s1.^epsil) + ...
    PHINext_s1s2.*((c2Next_s1s2.^(1-sig))./(1-sig) - alpha*(a2Next_s1s2.^epsil))...
    - PHINext_s1s1.*((c2Next_s1s1.^(1-sig))./(1-sig) -  alpha*(a2Next_s1s1.^epsil))+ ...
    funeval(parexp_planner1,fspace,PHINext_s1s2  );
exp_planner_disc_utility_next_s2s1 = (c1Next_s2s1.^(1-sig))./(1-sig)  - alpha*(a1Next_s2s1.^epsil) ...
    - (c1Next_s1s1.^(1-sig))./(1-sig) + alpha*(a1Next_s1s1.^epsil) + ...
    PHINext_s2s1.*((c2Next_s2s1.^(1-sig))./(1-sig) - alpha*(a2Next_s2s1.^epsil))...
    - PHINext_s1s1.*((c2Next_s1s1.^(1-sig))./(1-sig) -  alpha*(a2Next_s1s1.^epsil))+ ...
    funeval(parexp_planner1,fspace,PHINext_s2s1  );
exp_planner_disc_utility_next_s2s2 =  (c1Next_s2s2.^(1-sig))./(1-sig)  - alpha*(a1Next_s2s2.^epsil) ...
    - (c1Next_s1s1.^(1-sig))./(1-sig) + alpha*(a1Next_s1s1.^epsil) + ...
    PHINext_s2s2.*((c2Next_s2s2.^(1-sig))./(1-sig) - alpha*(a2Next_s2s2.^epsil))...
    - PHINext_s1s1.*((c2Next_s1s1.^(1-sig))./(1-sig) -  alpha*(a2Next_s1s1.^epsil))+ ...
    funeval(parexp_planner1,fspace,PHINext_s2s2  );




%%%%%%%%
%                       %
% EQUATIONS  %
%                       %
%%%%%%%%

% ICC agent 1
equ1 = epsil*alpha.*(a1.^(epsil-1))  -   prob2_s1.*(betta.*dprob1_s1.*exp_disc_utility1_next_s1s1 ...
    +  betta.*dprob1_s2.*exp_disc_utility1_next_s2s1) - ...
    prob2_s2.*(betta.*dprob1_s1.*exp_disc_utility1_next_s1s2 ...
    +  betta.*dprob1_s2.*exp_disc_utility1_next_s2s2) ;

% ICC agent 2
equ1bis = epsil*alpha.*(a2.^(epsil-1)) -   prob1_s1.*(betta.*dprob2_s1.*exp_disc_utility2_next_s1s1  ...
    +  betta.*dprob2_s2.*exp_disc_utility2_next_s1s2) ...
    -   prob1_s2.*(betta.*dprob2_s1.*exp_disc_utility2_next_s2s1  ...
    +  betta.*dprob2_s2.*exp_disc_utility2_next_s2s2) ;

% dynamic equat. for exp_disc_utility1
if sig==1
    equ2 = exp_disc_utility1 - log(c1) ...
        + alpha.*(a1.^epsil) ...
        -  betta.*prob2_s1.*(prob1_s1.*exp_disc_utility1_next_s1s1 ...
        +  prob1_s2.*exp_disc_utility1_next_s2s1 )...
        -  betta.*prob2_s2.*(prob1_s1.*exp_disc_utility1_next_s1s2 ...
        +  prob1_s2.*exp_disc_utility1_next_s2s2 );

    equ3 = exp_disc_utility2 - log(c2) ...
        + alpha.*(a2.^epsil) ...
        -  betta.*prob2_s1.*(prob1_s1.*exp_disc_utility2_next_s1s1 ...
        + prob1_s2.*exp_disc_utility2_next_s2s1 )...
        -  betta.*prob2_s2.*(prob1_s1.*exp_disc_utility2_next_s1s2 ...
        +  prob1_s2.*exp_disc_utility2_next_s2s2 );

else

    equ2 = exp_disc_utility1 - (c1.^(1-sig))./(1-sig) ...
        + alpha.*(a1.^epsil) ...
        -  betta.*prob2_s1.*(prob1_s1.*exp_disc_utility1_next_s1s1 ...
        +  prob1_s2.*exp_disc_utility1_next_s2s1 )...
        -  betta.*prob2_s2.*(prob1_s1.*exp_disc_utility1_next_s1s2 ...
        +  prob1_s2.*exp_disc_utility1_next_s2s2 );

    equ3 = exp_disc_utility2 - (c2.^(1-sig))./(1-sig) ...
        + alpha.*(a2.^epsil) ...
        -  betta.*prob2_s1.*(prob1_s1.*exp_disc_utility2_next_s1s1 ...
        +  prob1_s2.*exp_disc_utility2_next_s2s1 )...
        -  betta.*prob2_s2.*(prob1_s1.*exp_disc_utility2_next_s1s2 ...
        + prob1_s2.*exp_disc_utility2_next_s2s2 );
end;



incent1 =  -epsil*alpha.*(a1.^(epsil-1)) - epsil*(epsil-1)*alpha.*(a1.^(epsil-2)).*lambda1 ...
    +  betta.*lambda1.*prob2_s1.*(prob1_s1.*dlike_ratio1_s1.*exp_disc_utility1_next_s1s1 ...
    + prob1_s2.*dlike_ratio1_s2.*exp_disc_utility1_next_s2s1) ...
    +  betta.*lambda1.*prob2_s2.*(prob1_s1.*dlike_ratio1_s1.*exp_disc_utility1_next_s1s2 ...
    + prob1_s2.*dlike_ratio1_s2.*exp_disc_utility1_next_s2s2) ;

incent2 =  -epsil*alpha.*(a2.^(epsil-1)) - epsil*(epsil-1)*alpha.*(a2.^(epsil-2)).*lambda2 ...
    +  betta.*lambda2.*prob2_s1.*(prob1_s1.*dlike_ratio2_s1.*exp_disc_utility2_next_s1s1 ...
    + prob1_s2.*dlike_ratio2_s1.*exp_disc_utility2_next_s2s1) ...
    +  betta.*lambda2.*prob2_s2.*(prob1_s1.*dlike_ratio2_s2.*exp_disc_utility2_next_s1s2 ...
    + prob1_s2.*dlike_ratio2_s2.*exp_disc_utility2_next_s2s2) ;


equ4 =  -phi1.*incent2 ...
    - betta.*dprob2_s1.*(prob1_s1.*(1+lambda1.*like_ratio1_s1).*exp_planner_disc_utility_next_s1s1 ...
    +  prob1_s2.*(1+lambda1.*like_ratio1_s2).*exp_planner_disc_utility_next_s2s1 )...
    -  betta.*dprob2_s2.*(prob1_s1.*(1+lambda1.*like_ratio1_s1).*exp_planner_disc_utility_next_s1s2 ...
    +  prob1_s2.*(1+lambda1.*like_ratio1_s2).*exp_planner_disc_utility_next_s2s2 );

equ5 = exp_planner_disc_utility  - (c1.^(1-sig))./(1-sig) + alpha.*(a1.^epsil) + epsil*alpha.*(a1.^(epsil-1)).*lambda1  ...
    - phi1.*((c2.^(1-sig))./(1-sig) - alpha.*(a2.^epsil)  - epsil*alpha.*(a2.^(epsil-1)).*lambda2 ) ...
    -  betta.*prob2_s1.*(prob1_s1.*(1+lambda1.*like_ratio1_s1).*exp_planner_disc_utility_next_s1s1 ...
    +  prob1_s2.*(1+lambda1.*like_ratio1_s2).*exp_planner_disc_utility_next_s2s1 )...
    -  betta.*prob2_s2.*(prob1_s1.*(1+lambda1.*like_ratio1_s1).*exp_planner_disc_utility_next_s1s2 ...
    +  prob1_s2.*(1+lambda1.*like_ratio1_s2).*exp_planner_disc_utility_next_s2s2 );

equ6 = -incent1 ...
    - betta.*prob2_s1.*(dprob1_s1.*(1+lambda1.*like_ratio1_s1).*exp_planner_disc_utility_next_s1s1 ...
    +  dprob1_s2.*(1+lambda1.*like_ratio1_s2).*exp_planner_disc_utility_next_s2s1 )...
    -  betta.*prob2_s2.*(dprob1_s1.*(1+lambda1.*like_ratio1_s1).*exp_planner_disc_utility_next_s1s2 ...
    +  dprob1_s2.*(1+lambda1.*like_ratio1_s2).*exp_planner_disc_utility_next_s2s2 );




equ = [equ1; equ1bis; equ2; equ3;  equ4;equ5; equ6] ;

