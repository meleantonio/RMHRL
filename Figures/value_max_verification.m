function equ = value_max_verification(par00B,par00,bond_grid, Gridphi,fspaceB,fspace)
global alpha betta sig epsil nu 
global nk

parexp1B =                     par00B;

a1 =  linspace(.001,.9999,nk);%                 


bond1Next = repmat(unique(bond_grid(:,1)),1,length(bond_grid(:,1)))';

states = lawofmotion(par00, fspace,bond_grid(:,2:3));

phiNext1_s1 = states(:,1);
phiNext1_s2 = states(:,2);
zetaNext1 = states(:,3);
transf_income = states(:,4);

prob1_s1 =  a1.^nu;
prob1_s2 = 1 - prob1_s1;


c1= repmat(transf_income,1,length(unique(bond_grid(:,1)))) - bond1Next ...
    + repmat(bond_grid(:,1),1,length(unique(bond_grid(:,1))))./betta;


[brow,bcol] = size(bond1Next);

for i = 1: bcol
    exp_disc_utility1_next_s1(:,i) = funeval(parexp1B,fspaceB,[ bond1Next(:,i) phiNext1_s1  zetaNext1]);%
    exp_disc_utility1_next_s2(:,i) = funeval(parexp1B,fspaceB,[ bond1Next(:,i) phiNext1_s2   zetaNext1]);%
end


c1(find(c1<=0)) = NaN;

clear bond1Next bond1


%%%%%%%%%%
%                           %
% EQUATIONS   %
%                           %
%%%%%%%%%%

    u1 = (repmat(c1,1,nk).^(1-sig))./(1-sig);
    u2 = u1 - alpha.*kron(a1,ones(brow,bcol)).*epsil ; clear u1;
    u3 = u2 + betta.*kron(prob1_s1,exp_disc_utility1_next_s1) ; clear u2;
    u4 = u3 + betta.*kron(prob1_s2,exp_disc_utility1_next_s2); clear u3;
    equ2 = u4; clear u4;

equ2(find(isnan(equ2))) = -inf;


equ = equ2;
