
clear all;

global alpha betta sig epsil omega1 omega2 nu nu2 ;
global output_y s1 s2 elast1 elast2 delta1 delta2
global  rounds_approx OrderVector Order_vector


% PARAMETERS:
alpha =.05;
betta =.95;
sig = 2;
epsil = 2;

omega1 = .5;
omega2 = .5;

delta1 =.06;
delta2 = .06;

elast1 =.3;
elast2 = .3;

output_y =1;

nu   =.1;
nu2 =.1;

nba = 2;                              % number of states of the world
nphi = round(100000^(1/3));


% States of nature (= income y)
s1 = zeros(nba,1);
s1(1) =output_y-.45;
s1(2) = .45;

s2 = zeros(nba,1);
s2(1) =  output_y-.45;
s2(2) = .45;

% grid for phi
phibar = omega2/omega1;
phi_sd =omega2/omega1;
phi_min = phibar- .15*phi_sd;
phi_max = phibar+  .15*phi_sd;

% grid for k
k_ss = 2.12;
k_min =2;
k_max =4; 


OrderVector = [ 2  4  6 7 8 ;  2  4  6 7 8 ;2  4  6 7 8 ];
time_computation = zeros(1,length(OrderVector));
testing_max = zeros(1,length(OrderVector));
testing_norm = zeros(1,length(OrderVector));

for uu=1:length(OrderVector)
    
    
    rounds_approx = 2;
    Order_vector = [[2;2;2] OrderVector(:,uu)];
    
    disp(sprintf(' Grid is %g by %g by %g',Order_vector(1, rounds_approx),...
        Order_vector(2, rounds_approx),Order_vector(3, rounds_approx)));
    mainRSP;
    time_computation_vec(uu) = time_computation;
    testing_max_vec(uu) =  testing_max;
    testing_norm_vec(uu) = testing_norm;
    
end


%% Print table and save it
table_final_RSP = [OrderVector(1,:); time_computation_vec; testing_max_vec; testing_norm_vec]'
save table_final_RSP table_final_RSP -ascii -double -tabs
