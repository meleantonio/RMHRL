%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  program to analyze computational speed
%  and accuracy for the RMH model with hidden assets
%  solved with recursive Lagrangeans
%
%
%
%                       Antonio Mele
%                       March 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all;

global alpha betta sig epsil nu nu2 gam1;
global   s output_y  rounds_approx OrderVector Order_vector nk
global zeta_min zeta_max phi_min phi_max

% PARAMETERS:
maxits = 10000;
alpha =.5;
betta =  .95;
sig =2;
epsil = 2;

output_y =1;

nu   = .1;
nu2 =.1;

% outside option: The initial Pareto weight for the agent
gam1 =.5955;

nba = 2;                              % number of states of the world
nphi = round(sqrt(100000)); % gridpoints for testing accuracy



% States of nature (= income y)
s = zeros(nba,1);
s(1) = 0;
s(2) =output_y;

% Create grid extrema for Pareto weights
phibar = gam1;
phi_sd =.3*gam1;
phi_min = phibar-phi_sd;
phi_max = phibar+phi_sd;

% Create grid extrema for costate EE
zetabar= .422025;
zeta_min = 0;
zeta_max = zetabar;


OrderVector = [ 4 6 8 10 12 15 20;  4 6 8 10 12 15 20 ];
time_computation_vec = zeros(1,length(OrderVector));
testing_max_vec = zeros(1,length(OrderVector));
testing_norm_vec = zeros(1,length(OrderVector));

for uu=1:length(OrderVector)
    
    rounds_approx = 2;
    if uu < 3
        Order_vector = [[4;4] OrderVector(:,uu)];
    elseif uu < 6 
        rounds_approx = 3;
        Order_vector = [[4;4] [8;8] OrderVector(:,uu)];
    else
        rounds_approx = 4;
        Order_vector = [[4;4] [8;8] [12;12] OrderVector(:,uu)];
    end
    mainHA;
    time_computation_vec(uu) = time_computation;
    testing_max_vec(uu) =  testing_max;
    testing_norm_vec(uu) = testing_norm;
    
end

%% Print all and save table
table_final_HA = [OrderVector(1,:); time_computation_vec; testing_max_vec; testing_norm_vec]'
save table_final_HA table_final_HA -ascii -double -tabs
