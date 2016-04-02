clear all;
tic;

global alpha betta sig  epsil omega1 omega2 nu nu2 gam1 output_y s1 s2
global  rounds_approx OrderVector Order_vector 


% PARAMETERS:
nba = 2;                              % number of states of the world
nphi = 100000;
alpha =.5;
betta =.95;
sig = 2;
epsil = 2;
omega1 = .5;
omega2 = .5;
output_y =1;
nu   =.5;
nu2 =.5;


% States of nature (= income y)
s1 = zeros(nba,1);
s1(1) =.4;
s1(2) =output_y -.4;

s2 = zeros(nba,1);
s2(1) =.4;
s2(2) =output_y -.4;


OrderVector = [10; 15; 20; 30; 50; 100];
time_computation = zeros(1,length(OrderVector));
testing_max = zeros(1,length(OrderVector));
testing_norm = zeros(1,length(OrderVector));

for uu=1:length(OrderVector)
    
    rounds_approx = 2;
    Order_vector = [10 OrderVector(uu)];
    mainRSE;
    
end



table_final_RSE = [OrderVector'; time_computation; testing_max; testing_norm]
save table_final_RSE table_final_RSE -ascii -double -tabs
