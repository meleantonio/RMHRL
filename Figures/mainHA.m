%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %       
% This code solve a dynamic agency problem  with hidden borrowing and lending      %
% with the recursive Lagrangean approach            %
% and collocation method over Lagrangean FOCs      %
%                                                   %
%                                                   %
%                    Antonio Mele, March 2010     %
%                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;


for oo = 1: rounds_approx;

    RoundAppr = oo;

    

    % Display grid extrema
    disp(sprintf('    '));
    disp(sprintf('[%e , %e]',phi_min,phi_max));
    disp(sprintf('[%e , %e]',zeta_min,zeta_max));
    disp(sprintf('    '));


    % Range on which we approximate the solution:
    LowerBound = [phi_min; zeta_min   ];
    UpperBound = [phi_max; zeta_max  ];

    % number of grid points in each dimension of the state space
    Order = Order_vector(:,oo);
    if oo >= 2
        fspace_old = fspace;
    end

    disp(sprintf('  RoundAppr = %d ',RoundAppr));
    
    
    % Generate the approximation space with Miranda-Fackler Compecon
    % toolbox
    
%              approxtype = 'cheb';
    approxtype = 'spli';
    splineorder = [];
    if(strcmp(approxtype,'spli'))
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,splineorder);
    else
        fspace = fundefn(approxtype,Order,LowerBound,UpperBound,[]);
    end;
    %  fspace1= fspace;
    
    % Generate the grid
    % if RoundAppr <=10
        nodes=funnode(fspace);
        Gridphi = gridmake(nodes);
    % else
    %     disp(sprintf('clusterizing the grid...'));
    %     Gridphi_old = Gridphi;
    %     fspace = fspace_old;
    %     simul_HA; 
    %     [index_clusters,Gridphi] = kmeans([reshape(phi1_simu_real(:,1:end-1), ...
    %         random_generations*periods_simulations ,1) ...
    %         reshape( zeta_simu_real(:,1:end-1), ...
    %         random_generations*periods_simulations ,1)],nClust,...
    %         'emptyaction','drop','start','cluster');
    % end
        
%    fspace = fspace1;

    % Initial conditions for coefficients of the approximated policy
    % functions
    if (RoundAppr == 1)    
        
        load initial_conditions_HA.mat para1  parlambda1 parexp_planner1 parexp1
        OrderIC = 10;
        approxtypeIC = 'spli';
        splineorderIC = [];
        if(strcmp(approxtypeIC,'spli'))
            fspaceIC = fundefn(approxtypeIC,OrderIC,LowerBound(1),UpperBound(1),splineorderIC);
        else
            fspaceIC = fundefn(approxtypeIC,OrderIC,LowerBound(1),UpperBound(1),[]);
        end;
            
        c1 =sqrt(Gridphi(:,1)+Gridphi(:,2)./betta);
        a1 =funeval(para1,fspaceIC, Gridphi(:,1));
        lambda1 =funeval(parlambda1 , fspaceIC, Gridphi(:,1))  ;
        exp_disc_utility1 = funeval(parexp1,fspaceIC,Gridphi(:,1));
        exp_planner_disc_utility1  =  funeval(parexp_planner1,fspaceIC,Gridphi(:,1));
        exp_planner_disc_utility2  =s(1) -s(2) + exp_planner_disc_utility1 ;
    else
        c1 =funeval(parc1,fspace_old, Gridphi);
        a1 =funeval(para1,fspace_old, Gridphi);
        lambda1 =funeval(parlambda1 , fspace_old, Gridphi)  ;
        exp_disc_utility1 = funeval(parexp1,fspace_old,Gridphi);
        exp_planner_disc_utility1  = funeval(parexp_planner1,fspace_old,Gridphi);
        exp_planner_disc_utility2  =s(1) -s(2) + exp_planner_disc_utility1 ;
        eta1 = Gridphi(:,2)./betta + (  1 - Gridphi(:,1).*(c1.^(-sig)))./((-sig).*(c1.^(-sig-1))) ;
    end;


    % generate basis functions Chebbasis at Gridphi where ValueVec was defined:
    Basis = funbas(fspace,Gridphi);

    % give new names to the grid (just for convenience9
    phi1 = Gridphi(:,1);
    zeta1 = Gridphi(:,2);


    % obtain coefficients
    parlambda1 =  Basis\lambda1;
    para1 =  Basis\a1;
    parc1 =  Basis\c1;
    parexp_planner1 =  Basis\exp_planner_disc_utility1;
    parexp1 =  Basis\exp_disc_utility1 ;
    
    % coefficients' vector
    parpolicy = [parlambda1; para1 ; parexp_planner1 ;parexp1 ; parc1];

    %%%%%%%%%%%%%%%
    %                                                  %
    %               MAIN LOOP                %
    %                                                  %
    %%%%%%%%%%%%%%%

   
    % solve FOCs with Broyden's algorithm
    [opt_vec,info] =  broydn('HA_focs',parpolicy,1e-6,0,1,Gridphi,fspace);
    disp(sprintf(' info = %d',info));
    disp(sprintf('    '));

    % give a new name to the solution vector
    par00 = opt_vec;



    % reshape the solution vector and give new names to coefficients
    par0 = reshape(par00,length(Gridphi),5 );

    parlambda1 =               par0(:,1);
    para1 =                     par0(:,2);
    parexp_planner1 =      par0(:,3);
    parexp1 =                     par0(:,4);
    parc1 =    par0(:,5);
    
    
    % calculate policy functions on the grid with solution coefficients
    a1 =                     funeval(para1,fspace, Gridphi);
    lambda1 = funeval(parlambda1,fspace, Gridphi);
    exp_planner_disc_utility1 = funeval(parexp_planner1,fspace,Gridphi);
    exp_planner_disc_utility2 = s(1) -s(2) + exp_planner_disc_utility1 ;
    exp_disc_utility1 = funeval(parexp1,fspace,Gridphi);
    c1 =                     funeval(parc1,fspace, Gridphi);


end;

% update coefficients' vector
parpolicy = opt_vec;


toc;
time_computation = toc;

% Calculate time needed for solution
time_hours = fix(toc/3600);
time_minutes = fix((toc - time_hours*3600)/60);
time_seconds = fix(toc -  time_hours*3600 - time_minutes*60);

disp(sprintf('    '));disp(sprintf('Time for solving the model was    '));
disp(sprintf('%d hours     %d  minutes %d  seconds',time_hours,time_minutes,time_seconds));



%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                %
%                               TESTING                                  %
%                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%


% Testing the approximation: Euler residuals at points off the grid:

phigrid_nodes_test = linspace(phi_min,phi_max,nphi)';%
zetagrid_nodes_test = linspace(zeta_min,zeta_max,nphi)';%

GridTest = gridmake(phigrid_nodes_test,zetagrid_nodes_test);

test_residuals =  HA_focs(parpolicy,GridTest,fspace);

norm_test  = norm(test_residuals);

max_test   = max(abs(test_residuals));


testing_max = max_test(end);
testing_norm = norm_test(end);

