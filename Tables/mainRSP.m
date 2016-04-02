%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code solve a two-agents risk sharing problem
% in a production economy with two-sided moral hazard
% with the recursive Lagrangean approach
% and collocation method over Lagrangean FOCs
%
%
%                    Antonio Mele, March 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

for oo = 1:rounds_approx

    RoundAppr = oo;

    phibar = omega2/omega1;

    phi_sd =omega2/omega1;

    phi_min = phibar- .5*phi_sd;
    phi_max = phibar+  .5*phi_sd;

    k_ss = 2.12;
    k_min =3;
    k_max =1.52*k_ss;


    % Range on which we approximate the solution:
    LowerBound = [phi_min; k_min; k_min ];
    UpperBound = [phi_max;  k_max; k_max];

    Order = Order_vector(:,oo);
    if oo >= 2
        fspace_old = fspace;
    end



    disp(sprintf('  RoundAppr = %d ',RoundAppr));


    approxtype = 'cheb';
    %                 approxtype = 'spli';
    %                 splineorder =[];%1;%  1;
    if(strcmp(approxtype,'spli'))
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,splineorder);
    else
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,[]);
    end;
    nodes = funnode(fspace);
    Gridphi = gridmake(nodes);





    if (RoundAppr == 1)

        fspace1 = fspace;
        Gridphi1 = Gridphi;
        load initial_conditions_RSP park1 park2 para1 para2 parlambda1 parlambda2 parexp1  parexp2 parexp_planner1 fspace;%


        mmm = length(Gridphi);
        ggg = 1:1:mmm;
        nnn = mmm +1 - ggg';
        Gridphi_inv = Gridphi(nnn,1);

        a1 =funeval(para1,fspace, Gridphi);
        a2 = funeval(para2,fspace, Gridphi)  ;
        lambda1 =funeval(parlambda1 , fspace, Gridphi)  ;
        lambda2 =funeval(parlambda2 , fspace, Gridphi)  ;
        exp_disc_utility1 = funeval(parexp1,fspace,Gridphi);
        exp_disc_utility2 = funeval(parexp2,fspace,Gridphi);
        exp_planner_disc_utility  = funeval(parexp_planner1,fspace,Gridphi);
        k1_next = funeval(park1,fspace,Gridphi);
        k2_next = funeval(park2,fspace,Gridphi);


        fspace=fspace1;
        Gridphi = Gridphi1;




    else

        a1 =funeval(para1,fspace_old, Gridphi);
        a2 = funeval(para2,fspace_old, Gridphi)  ;
        lambda1 =funeval(parlambda1 , fspace_old, Gridphi)  ;
        lambda2 =funeval(parlambda2 , fspace_old, Gridphi)  ;
        exp_disc_utility1 = funeval(parexp1,fspace_old,Gridphi);
        exp_disc_utility2 = funeval(parexp2,fspace_old,Gridphi);
        exp_planner_disc_utility  = funeval(parexp_planner1,fspace_old,Gridphi);
        k1_next = funeval(park1,fspace_old,Gridphi);
        k2_next = funeval(park2,fspace_old,Gridphi);

    end;


    % generate basis functions Chebbasis at Gridphi where ValueVec was defined:
    Basis = funbas(fspace,Gridphi);%myfunbas(fspace,Gridphi);


    phi1 = Gridphi(:,1);

    parlambda1 =  Basis\lambda1;
    parlambda2 =  Basis\lambda2;
    para1 =  Basis\a1;
    para2 =  Basis\a2;
    parexp_planner1 =  Basis\exp_planner_disc_utility;
    parexp1 =  Basis\exp_disc_utility1 ;
    parexp2 =  Basis\exp_disc_utility2 ;
    park1 = Basis\k1_next;
    park2 = Basis\k2_next;

    parpolicy = [parlambda1; para1 ; parexp1 ; parexp2; parlambda2; para2; ...
        parexp_planner1;park1 ;park2  ];


    %%%%%%%%%%%%%%%
    %                                                  %
    %               MAIN LOOP                %
    %                                                  %
    %%%%%%%%%%%%%%%

    [opt_vec,info] =  broydn('RSP_focs',parpolicy,1e-6,0,1,Gridphi,fspace);
    disp(sprintf(' info = %d',info));
    disp(sprintf('    '));

    par00 = opt_vec;

    par0 = reshape(par00,length(Gridphi),9 );

    parlambda1 =               par0(:,1);
    para1 =                     par0(:,2);
    parexp_planner1 =      par0(:,7);
    parexp1 =                     par0(:,3);
    parexp2 =                     par0(:,4);
    parlambda2 =               par0(:,5);
    para2 =                     par0(:,6);
    park1 =                           par0(:,8);
    park2 =                     par0(:,9);



    a1 =                     funeval(para1,fspace, Gridphi);
    a2 =                     funeval(para2,fspace, Gridphi);
    lambda1 = funeval(parlambda1,fspace, Gridphi);
    lambda2 = funeval(parlambda2,fspace, Gridphi);
    exp_planner_disc_utility = funeval(parexp_planner1,fspace,Gridphi);
    exp_disc_utility1 = funeval(parexp1,fspace,Gridphi);
    exp_disc_utility2 = funeval(parexp2,fspace,Gridphi);
    k1_next = funeval(park1,fspace,Gridphi);
    k2_next = funeval(park2,fspace,Gridphi);

    aggr_output = zeros(length(Gridphi),4);
    aggr_output(:,1) = - k1_next - k2_next +  (1-delta1).*Gridphi(:,2) + (1-delta2).*Gridphi(:,3)  +...
        s1(2).*(Gridphi(:,2).^elast1) + s2(2).*(Gridphi(:,3).^elast2);
    aggr_output(:,2) = - k1_next - k2_next +  (1-delta1).*Gridphi(:,2) + (1-delta2).*Gridphi(:,3)  +...
        s1(1).*(Gridphi(:,2).^elast1) + s2(2).*(Gridphi(:,3).^elast2);
    aggr_output(:,3) = - k1_next - k2_next +  (1-delta1).*Gridphi(:,2) + (1-delta2).*Gridphi(:,3)  +...
        s1(2).*(Gridphi(:,2).^elast1) + s2(1).*(Gridphi(:,3).^elast2);
    aggr_output(:,4) = - k1_next - k2_next +  (1-delta1).*Gridphi(:,2) + (1-delta2).*Gridphi(:,3)  +...
        s1(1).*(Gridphi(:,2).^elast1) + s2(1).*(Gridphi(:,3).^elast2);

    c1 = aggr_output(:,4)./( 1 + phi1.^(1/sig)); % ((omega1 + phi1)./eta).^sig; %

    c2 = aggr_output(:,4) - c1;


    parpolicy = opt_vec;


end;

toc;

time_computation(uu) = toc;

%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                %
%                               TESTING                                  %
%                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%


% Testing the approximation: Euler residuals at points off the grid:

phi_test= linspace(LowerBound(1),UpperBound(1),nphi)';%
k1_test= linspace(LowerBound(2),UpperBound(2),nphi)';%
k2_test= linspace(LowerBound(3),UpperBound(3),nphi)';%
GridTest = gridmake(phi_test,k1_test,k2_test);

test_residuals =  RSP_focs(parpolicy,GridTest,fspace);

norm_test  = norm(test_residuals);

max_test   = max(abs(test_residuals));



testing_max(uu) = max_test(end);
testing_norm(uu) = norm_test(end);

