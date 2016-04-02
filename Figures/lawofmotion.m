% generate law of motion for phi and zeta
% in the verification algorithm


function z = lawofmotion(par00, fspace,Gridphi)

global betta sig  nu ;


par0 = reshape(par00(1:end),length(par00)/5,5 );

parlambda1 =               par0(:,1);
para1 =                         par0(:,2);
% parexp_planner1 =      par0(:,3);
% parexp1 =                     par0(:,4);
parc1 =    par0(:,5);

a1 =                     funeval(para1,fspace, Gridphi);
lambda1 = funeval(parlambda1,fspace, Gridphi);
c1 =                     funeval(parc1,fspace, Gridphi);
phi1 = Gridphi(:,1);
zeta1 = Gridphi(:,2);
eta1 =   zeta1./betta + (  1 - phi1.*(c1.^(-sig)))./((-sig).*(c1.^(-sig-1)))   ;   %funeval(pareta1,fspace, Gridphi); %

prob1_s1 =  a1.^nu;
prob1_s2 = 1 - prob1_s1;

dprob1_s1 =  nu.*(a1.^(nu-1));
dprob1_s2 =  - dprob1_s1;


% likelihood ratios
like_ratio1_s1 =  dprob1_s1./prob1_s1;
like_ratio1_s2 = dprob1_s2./prob1_s2;

% state variables
phiNext1_s1 = phi1 + lambda1.*like_ratio1_s1;
phiNext1_s2 =phi1 + lambda1.*like_ratio1_s2;
zetaNext1 = eta1;

z = [phiNext1_s1 phiNext1_s2 zetaNext1 c1];
