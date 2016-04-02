% verification procedure


bond_min =0;
bond_max = 1;

%  Range on which we approximate the solution:
LowerBoundB = [bond_min ; phi_min; zeta_min  ];%
UpperBoundB = [bond_max; phi_max; zeta_max  ]; %

OrderB = [10; 10; 10]; % fix number of gridpoints for state variables
nk = 1000; % number of gridpoints for effort

% we use linear interpolation of linear splines
approxtypeB = 'lin'; %'spli'; %
splineorderB = 1; %[];% 1;%

if(strcmp(approxtypeB,'spli'))
    fspaceB = fundefn(approxtypeB,OrderB,LowerBoundB,UpperBoundB,splineorderB);
else
    fspaceB = fundefn(approxtypeB,OrderB,LowerBoundB,UpperBoundB,[]);
end;
nodesB = funnode(fspaceB);
bond_grid = gridmake(nodesB);


% generate basis functions Chebbasis at Gridphi where ValueVec was defined:
BasisB = funbas(fspaceB,bond_grid);



exp_disc_utility1B = funeval(parexp1,fspace,bond_grid(:,2:3));

parexp1B =  BasisB\exp_disc_utility1B ;
parpolicyB = [ parexp1B ];

U = ones(length(bond_grid),1);

for it = 1:maxits  % value function iteration
    parexp1B_old = parexp1B;
    value_old = U;

    % maximize value function for current basis
    % coefficients and other parameters
    value = value_max_verification(parexp1B,parpolicy,...
        bond_grid, Gridphi,fspaceB,fspace);
    [U,opt_vecB] = max(value,[],2);

    parexp1B = BasisB\U;

    disp(sprintf('iteration = %d Diff = %e norm = %e ',...
        it,(max(abs((U-value_old)./value_old))),norm(parexp1B-parexp1B_old)));
    if norm(parexp1B-parexp1B_old)<1e-6, break, end;

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
