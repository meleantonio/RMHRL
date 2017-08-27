% Verification procedure
% (using an idea from Abraham and Pavoni (2008))

% Range for bonds
bond_min =-10;
bond_max = 10;

%  Range on which we approximate the solution:
LowerBoundB = [bond_min ; phi_min; zeta_min  ];%
UpperBoundB = [bond_max; phi_max; zeta_max  ]; %

OrderB = [10; 10; 10]; % fix number of gridpoints for state variables
nk = 2000; % number of gridpoints for effort

% we use linear interpolation of linear splines
approxtypeB = 'lin'; 
splineorderB = [];

if(strcmp(approxtypeB,'spli'))
    fspaceB = fundefn(approxtypeB,OrderB,LowerBoundB,UpperBoundB,splineorderB);
else
    fspaceB = fundefn(approxtypeB,OrderB,LowerBoundB,UpperBoundB,[]);
end;
nodesB = funnode(fspaceB);
bond_grid = gridmake(nodesB);


% generate basis functions BasisB at bond_grid:
BasisB = funbas(fspaceB,bond_grid);

% evaluate continuation value for the agent using the solution from the Lagrangean approach
exp_disc_utility1B = funeval(parexp1,fspace,bond_grid(:,2:3));

% calculate initial conditions for the interpolation coefficients
parexp1B =  BasisB\exp_disc_utility1B ;
parpolicyB = [ parexp1B ];

U = ones(length(bond_grid),1); % initial guess for the value function

for it = 1:maxits  % value function iteration
    parexp1B_old = parexp1B;
    value_old = U;

    % maximize value function for current basis
    % coefficients and other parameters
    value = value_max_verification(parexp1B,parpolicy,...
        bond_grid, Gridphi,fspaceB,fspace);
    [U,opt_vecB] = max(value,[],2);

    % homotopy step for stability
    parexp1B = .5*(BasisB\U) + .5*parexp1B_old ; %BasisB\U;%

    % display convergence metrics
    disp(sprintf('iteration = %d Diff = %e norm = %e ',...
        it,(max(abs((U-value_old)./value_old))),norm(abs((parexp1B-parexp1B_old)./parexp1B_old))));
    if norm(abs((parexp1B-parexp1B_old)./parexp1B_old))<1e-6, break, end;

end

% Testing the approximation: Euler residuals at points off the grid:
nB = 100;
phigrid_nodes_testB = linspace(phi_min,phi_max,nB)';%
zetagrid_nodes_testB = linspace(zeta_min,zeta_max,nB)';%
bond_nodes_testB = linspace(bond_min,bond_max,nB)';%

GridTest = gridmake(bond_nodes_testB,phigrid_nodes_testB,zetagrid_nodes_testB);


test_index = [1:1:length(GridTest)]';
test_index = reshape(test_index,length(GridTest)/100,100);
max_testB = -Inf;
for ff=1:size(test_index,2)

    testB = funeval(parexp1B,fspaceB,GridTest(ff,:)) - max(value_max_verification(parexp1B,parpolicy,...
        GridTest(ff,:), Gridphi,fspaceB,fspace),[],2); % U; %
    max_testB = max(abs(testB), max_testB);

end

disp(' ');
disp(sprintf('test verification = %g',max_testB));


% calculate the welfare change wrt the solution obtained with the
% Lagrangean approach
BOND0 = -betta*(funeval(parexp_planner1,fspace, [gam1 0])   - gam1*funeval(parexp1,fspace,[gam1 0])) ;%
exp_disc_utilityBOND0 = funeval(parexp1B,fspaceB,[ BOND0 gam1 0]);

delta_welfare = (exp_disc_utilityBOND0 - funeval(parexp1,fspace,[gam1 0]))...
    /abs(funeval(parexp1,fspace,[gam1 0]));



disp(sprintf('    '));
disp(sprintf(' The gain in welfare is  %g',delta_welfare));
disp(sprintf('    '));
if exp_disc_utilityBOND0 < funeval(parexp1,fspace,[gam1 0])
    disp(sprintf('Looks like the FOA is valid!'));
else
    disp(sprintf('Oops... FOA is not valid!'));
end