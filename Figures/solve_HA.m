%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Code solving and generating figures for the RMH model with hidden assets
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

% parameters for simulations
random_generations = 1000;%            % random_generations multiplied by metit gives
metit = 50;                                           % the total number of draws used; just for
                                                            % technical reasons, we store
                                                            % the draws on big matrices
                                                            % with no more than 1000 rows.
                                                            % For simulating one series,
                                                            % set both to 1
periods_simulations =  200; % number of periods for the simulation



rounds_approx = 2;
Order_vector = [[4;4] [15;10]];%[6;6]

mainHA;
testing_max
% verification_procedure;

figures_HA;

