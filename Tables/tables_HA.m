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
global   s output_y  rounds_approx OrderVector Order_vector 


% PARAMETERS:
maxits = 10000; 
alpha =.50;
betta =  .95;
sig =2;
epsil = 2;

output_y =1;

nu   =.5;
nu2 =.5;

% outside option: The initial Pareto weight for the agent
gam1 =.5955;

nba = 2;                              % number of states of the world
nphi = round(sqrt(100000)); % gridpoints for testing accuracy



% States of nature (= income y)
s = zeros(nba,1);
s(1) = 0;
s(2) =output_y;




OrderVector = [ 4 6 8 10 12 15 20;  4 6 8 10 12 15 20 ];
time_computation = zeros(1,length(OrderVector));
testing_max = zeros(1,length(OrderVector));
testing_norm = zeros(1,length(OrderVector));

for uu=1:length(OrderVector)
    
    rounds_approx = 2;
    Order_vector = [[4;4] OrderVector(:,uu)];
    mainHA;
    
end


disp(sprintf('%g    %g      %g      %g  ',OrderVector', time_computation, testing_max, testing_norm));


table_final_HA = [OrderVector; time_computation; testing_max; testing_norm]
save table_final_HA table_final_HA -ascii -double -tabs
