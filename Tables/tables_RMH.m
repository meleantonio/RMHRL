%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  program to analyze computational speed
% and accuracy for the RMH solved with
% recursive Lagrangeans
%
%
%
%                       Antonio Mele
%                       February 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;


global alpha betta sig phibar   epsil  nu nu2 gam1;
global phi_min phi_max s  ;
global   phi1  RoundAppr sig2
global  output_y perst_weight rounds_approx OrderVector Order_vector 

% PARAMETERS:

nba = 2;                        % number of states of the world
alpha =.50;                     % disutility of effort is alpha*v(a)
betta =  .95;                   % discount factor
sig = 2; %                  % CRRA utility parameter (c^(1-sig))/(1-sig)
epsil = 2;%                        % v(a) = a^epsil
output_y =1;                    % endowment in good times (in bad times is zero)
nphi = 100000;                     % number of gridpoints used for testing accuracy
perst_weight =0;                % parameter that governs the persistent of endowment: 
                                          %  if 0 no persistence, if 1 there is persistence

nu   = .1;                   % parameter of the prob(y_t+1= y^H| a_t) = a^nu;
nu2 = .1;                    % I allow for the possibility that transition probabilities are
                                % different in different states of the world

% outside option: The initial Pareto weight for the agent
gam1 =.5955;

% States of nature 
s = zeros(nba,1);
s(1) =  0;
s(2) =output_y;

% create a grid for costates around the first pareto weight for agent
phibar = gam1;
phi_min = .3;
phi_max = 1.5;

OrderVector = [10; 15; 20; 30; 50; 100];
time_computation_vec = zeros(1,length(OrderVector));
testing_max_vec = zeros(1,length(OrderVector));
testing_norm_vec = zeros(1,length(OrderVector));

for uu=1:length(OrderVector)
    
    rounds_approx = 2;
    Order_vector = [10 OrderVector(uu)];
    mainRMH;
    time_computation_vec(uu) = time_computation;
    testing_max_vec(uu) =  testing_max;
    testing_norm_vec(uu) = testing_norm;

    
end


%% Print all and save table
table_final_RMH = [OrderVector'; time_computation_vec; testing_max_vec; testing_norm_vec]'
save table_final_RMH table_final_RMH -ascii -double -tabs
