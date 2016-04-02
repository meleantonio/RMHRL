function equ = HA_focs(par00, Gridphi,fspace);

global alpha betta sig  epsil nu
global   s output_y


% coefficients
par0 = reshape(par00(1:end),length(par00)/5,5 );

parlambda1 =               par0(:,1);
para1 =                         par0(:,2);
parexp_planner1 =      par0(:,3);
parexp1 =                     par0(:,4);
parc1 =    par0(:,5);


% new names to the grid (just for convenience in writing the equations)
phi1 = Gridphi(:,1);
zeta1 = Gridphi(:,2);

% allocations and Lagrange multipliers
a1 =                     funeval(para1,fspace, Gridphi);
lambda1 = funeval(parlambda1,fspace, Gridphi);
c1 =                     funeval(parc1,fspace, Gridphi);
eta1 =   zeta1./betta + (  1 - phi1.*(c1.^(-sig)))./((-sig).*(c1.^(-sig-1)))   ;

% probabilities
prob1_s1 =  a1.^nu;
prob1_s2 = 1 - prob1_s1;

% derivatives of probabilities
dprob1_s1 =  nu.*(a1.^(nu-1));
dprob1_s2 =  - dprob1_s1;

% likelihood ratios
like_ratio1_s1 =  dprob1_s1./prob1_s1;
like_ratio1_s2 = dprob1_s2./prob1_s2;

% derivatives of likelihood ratios
dlike_ratio1_s1 =  -nu./(a1.^2);%
dlike_ratio1_s2 = (-nu*(nu-1)*(a1.^(nu-2)) - nu*(a1.^(2*nu-2)))./(prob1_s2.^2) ;


% state variables at t+1
phiNext1_s1 = phi1 + lambda1.*like_ratio1_s1;
phiNext1_s2 =phi1 + lambda1.*like_ratio1_s2;
zetaNext1 = eta1;


% consumption at t+1
c1Next_s1 = funeval(parc1,fspace,[phiNext1_s1 zetaNext1]); %
c1Next_s2 = funeval(parc1,fspace,[phiNext1_s2 zetaNext1]);%

% effort at t+1
a1Next_s1 = funeval(para1,fspace,[phiNext1_s1 zetaNext1]);
a1Next_s2 = funeval(para1,fspace,[phiNext1_s2 zetaNext1]);


% agent´s utility at t+1
if sig==1
    utilityNext_s1 = log(c1Next_s1) - alpha*(a1Next_s1.^epsil);
    utilityNext_s2 = log(c1Next_s2) - alpha*(a1Next_s2.^epsil);
else
    utilityNext_s1 = (c1Next_s1.^(1-sig))./(1-sig) - alpha*(a1Next_s1.^epsil);
    utilityNext_s2 = (c1Next_s2.^(1-sig))./(1-sig) - alpha*(a1Next_s2.^epsil);
end;

% planner continuation value at t and t+1
exp_planner_disc_utility1 = funeval(parexp_planner1,fspace,Gridphi);
exp_planner_disc_utility1_next_s1 = funeval(parexp_planner1,fspace,[phiNext1_s1 zetaNext1] );
exp_planner_disc_utility2_next_s2 =  s(1) -s(2) +funeval(parexp_planner1,fspace,[phiNext1_s2 zetaNext1]);

% agent continuation value at t and t+1
exp_disc_utility1 = funeval(parexp1,fspace,Gridphi);
exp_disc_utility1_next_s1 = funeval(parexp1,fspace,[phiNext1_s1  zetaNext1]);
exp_disc_utility1_next_s2 = funeval(parexp1,fspace,[phiNext1_s2   zetaNext1]);


%%%%%%%%
%                       %
% EQUATIONS  %
%                       %
%%%%%%%%


% ICC
equ1 = epsil*alpha.*(a1.^(epsil-1)) -   betta.*(dprob1_s1.*exp_disc_utility1_next_s1 ...
    +  dprob1_s2.*exp_disc_utility1_next_s2)  ;

% Bellman equation for the agent
if sig==1
    equ2 = exp_disc_utility1 - log(c1) ...
        + alpha.*(a1.^epsil) ...
        -  betta.*prob1_s1.*exp_disc_utility1_next_s1 ...
        -  betta.*prob1_s2.*exp_disc_utility1_next_s2;
else
    equ2 = exp_disc_utility1 - (c1.^(1-sig))./(1-sig) ...
        + alpha.*(a1.^epsil) ...
        -  betta.*prob1_s1.*exp_disc_utility1_next_s1 ...
        -  betta.*prob1_s2.*exp_disc_utility1_next_s2;
end;


% FOC a1
equ3 =  -epsil*alpha.*(a1.^(epsil-1)).*phi1 - epsil*(epsil-1)*alpha.*(a1.^(epsil-2)).*lambda1 ...
    + ( betta.*lambda1.*(prob1_s1.*dlike_ratio1_s1.*utilityNext_s1 ...
    + prob1_s2.*dlike_ratio1_s2.*utilityNext_s2) ...
    +   betta.*dprob1_s1.*exp_planner_disc_utility1_next_s1 ...
    +   betta.*dprob1_s2.*exp_planner_disc_utility2_next_s2);

% Bellman equation for the planner
if sig==1
    equ4 = exp_planner_disc_utility1 - (output_y - c1 - s(1) )   ...
        - phi1.*( log(c1)  - alpha.*(a1.^epsil)  ) ...
        +  epsil*alpha.*(a1.^(epsil-1)).*lambda1 - (eta1 - zeta1./betta).*(c1.^(-sig))...
        - (betta.*prob1_s1.*exp_planner_disc_utility1_next_s1 ...
        + betta.*prob1_s2.*exp_planner_disc_utility2_next_s2);
else
    equ4 = exp_planner_disc_utility1 - (output_y - c1 - s(1) )   ...
        - phi1.*(  (c1.^(1-sig))./(1-sig)  - alpha.*(a1.^epsil)  ) ...
        +  epsil*alpha.*(a1.^(epsil-1)).*lambda1 - (eta1 - zeta1./betta).*(c1.^(-sig))...
        - (betta.*prob1_s1.*exp_planner_disc_utility1_next_s1 ...
        + betta.*prob1_s2.*exp_planner_disc_utility2_next_s2);
end;

% Euler equation
equ5 = (c1.^(-sig))- prob1_s1.*(c1Next_s1.^(-sig)) - prob1_s2.*(c1Next_s2.^(-sig));

% all equations together
equ = [equ1; equ2 ; equ3; equ4;equ5];

% avoid the solution is strange
if(any(a1<0)) || (any(c1<1e-6))  ||  (any(a1>1))   || (any( phiNext1_s1<0)) || (any( phiNext1_s2<0))
    equ(1) = 1e100;
end;
