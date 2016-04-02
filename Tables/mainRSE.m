%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code solve a two-agents risk sharing problem
% in an endowment economy with two-sided moral hazard
% with the recursive Lagrangean approach
% and collocation method over Lagrangean FOCs
%
%
%                    Antonio Mele, March 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic;

for oo = 1: rounds_approx

    RoundAppr = oo;



    phibar = omega2/omega1;
    %         phi_sd =omega2/omega1;

    % without homogen
    phi_min = phibar  -.25;
    phi_max = phibar+ .25;



    aggr_output = 2*output_y - [s1(2)+s2(2); s1(1)+s2(2); s1(2)+s2(1); s1(1)+s2(1)];

    %  Range on which we approximate the solution:
    LowerBound = [phi_min ];
    UpperBound = [phi_max];

    Order = Order_vector(oo);
    if oo >= 2
        fspace_old = fspace;
    end



    disp(sprintf('  RoundAppr = %d ',RoundAppr));

    % approxtype = 'cheb';
    approxtype = 'spli';
    splineorder = [];% 1;%
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
        load initial_conditions_RSE para1 parlambda1 parexp1  parexp2 parexp_planner1 ...
            Gridphi fspace;%

        Gridphi_inv = [Gridphi(10); Gridphi(9);  Gridphi(8);  Gridphi(7);  Gridphi(6);  ...
            Gridphi(5);  Gridphi(4);  Gridphi(3);  Gridphi(2); Gridphi(1) ];

        % agent 1
        a1 =funeval(para1,fspace,[Gridphi(:,1)]);
        lambda1 =funeval(parlambda1 , fspace, [Gridphi(:,1) ]);
        exp_disc_utility1 = funeval(parexp1,fspace,[Gridphi(:,1)]);

        % agent 2
        exp_disc_utility2 = funeval(parexp2,fspace,[Gridphi_inv ]);
        lambda2 = funeval(parlambda1 , fspace, [Gridphi_inv]);
        a2 =funeval(para1,fspace,[Gridphi_inv]);%
        exp_planner_disc_utility  =exp_disc_utility1 + [Gridphi(:,1) ].*exp_disc_utility2 -...
            lambda1.*(a1.^2)./(1-betta) - [Gridphi(:,1)].*(lambda2.*(a2.^2))./(1-betta) ;

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
    end;


    % generate basis functions Basis
    Basis = funbas(fspace,Gridphi);


    phi1 = Gridphi(:,1);


    parlambda1 =  Basis\lambda1;
    parlambda2 =  Basis\lambda2;
    para1 =  Basis\a1;
    para2 =  Basis\a2;
    parexp_planner1 =  Basis\exp_planner_disc_utility;
    parexp1 =  Basis\exp_disc_utility1 ;
    parexp2 =  Basis\exp_disc_utility2 ;
    parpolicy = [parlambda1; para1 ; parexp1 ; parexp2; parlambda2; ...
        para2; parexp_planner1 ];


    %%%%%%%%%%%%%%%
    %
    %               MAIN LOOP
    %
    %%%%%%%%%%%%%%%


    [opt_vec,info] =  broydn('RSE_focs',parpolicy,1e-6,0,1,Gridphi,fspace);
    disp(sprintf(' info = %d',info));
    disp(sprintf('    '));

    par00 = opt_vec;

    par0 = reshape(par00,length(Gridphi),7 );

    parlambda1 =               par0(:,1);
    para1 =                     par0(:,2);
    parexp_planner1 =      par0(:,7);
    parexp1 =                     par0(:,3);
    parexp2 =                     par0(:,4);
    parlambda2 =               par0(:,5);
    para2 =                     par0(:,6);


    a1 =                     funeval(para1,fspace, Gridphi);
    a2 =                     funeval(para2,fspace, Gridphi);
    lambda1 = funeval(parlambda1,fspace, Gridphi);
    lambda2 = funeval(parlambda2,fspace, Gridphi);
    exp_planner_disc_utility = funeval(parexp_planner1,fspace,Gridphi);
    exp_disc_utility1 = funeval(parexp1,fspace,Gridphi);
    exp_disc_utility2 = funeval(parexp2,fspace,Gridphi);
    c1 = aggr_output(4)./( 1 + phi1.^(1/sig));
    c2 = aggr_output(4) - c1;


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

grid_nodes_test = linspace(LowerBound,UpperBound,nphi)';%
GridTest = gridmake(grid_nodes_test);

test_residuals =  RSE_focs(parpolicy,GridTest,fspace);

norm_test  = norm(test_residuals);

max_test   = max(abs(test_residuals));



testing_max(uu) = max_test(end);
testing_norm(uu) = norm_test(end);


