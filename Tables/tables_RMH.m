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
global   phi1  RoundAppr 
global  gam1_2 output_y perst_weight rounds_approx OrderVector Order_vector 

% PARAMETERS:

nba = 2;                        % number of states of the world
alpha =.50;                     % disutility of effort is alpha*v(a)
betta =  .95;                   % discount factor
sig = 2;                        % CRRA utility parameter (c^(1-sig))/(1-sig)
epsil = 2;                      % v(a) = a^epsil
output_y =1;                    % endowment in good times (in bad times is zero)
nphi = 100000;                     % number of gridpoints used for testing accuracy
perst_weight =0;                % parameter that governs the persistent of endowment: if 0 no persistence, if 1 there is persistence

nu   =.5;                       % parameter of the prob(y_t+1= y^H| a_t) = a^nu;
nu2 = .5;                       % I allow for the possibility that transition probabilities are
                                % different in different states of the world

% outside option: The initial Pareto weight for the agent
gam1 =.5955;
% gam1_2 = [gam1; gam1];

% States of nature 
s = zeros(nba,1);
s(1) =  0;
s(2) =output_y;



OrderVector = [10; 15; 20; 30; 50; 100];
time_computation = zeros(1,length(OrderVector));
testing_max = zeros(1,length(OrderVector));
testing_norm = zeros(1,length(OrderVector));

for uu=1:length(OrderVector)
    
    rounds_approx = 2;
    Order_vector = [10 OrderVector(uu)];
    mainRMH;
    
end


disp(sprintf('%g    %g      %g      %g  ',OrderVector', time_computation, testing_max, testing_norm));


table_final_RMH = [OrderVector'; time_computation; testing_max; testing_norm]
save table_final_RMH table_final_RMH -ascii -double -tabs
