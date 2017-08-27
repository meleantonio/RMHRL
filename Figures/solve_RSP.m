%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Code solving and generating figures for the risk sharing model with 
%  two-sided moral hazard in a production economy
%  solved with recursive Lagrangeans
%
%                       Antonio Mele
%                       March 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



% parameters for simulations
% we set them such that we simulate a sample path first
random_generations = 1;              	% random_generations multiplied by metit gives
metit = 1;                                          % the total number of draws used; just for
                                                        % technical reasons, we store
                                                        % the draws on big matrices
                                                        % with no more than 1000 rows.
                                                        % For simulating one series,
                                                        % set both to 1
periods_simulations =  200; % number of periods for the simulation

% Rounds of approximation and corresponding order
rounds_approx = 4;
Order_vector = [ [2;2;2] [ 4  ;   4   ; 4   ] [6;6;6] [8;8;8]]; 

mainRSP;
testing_max

% save solution
save paramaters.mat parlambda1 para1  parexp1  parexp2 parlambda2 para2 ...
        parexp_planner1 park1 park2  fspace;
    

% simulations:

% initial conditions for Pareto weights
omega2 = 1; 
omega1 =1;

% initial conditions for capital
k1_zero =  2.1;
k2_zero =  2.1;

avg = 0; % this is counter for figures' number

% parameter for using random shocks or predetermined ones: if 0, random shocks, o/w shocks are given by g1_predet and g2_predet (defined below)
diverge = 0; % 1; %

g1_predet = s1(1).*ones(random_generations, periods_simulations);
g2_predet = s2(1).*ones(random_generations, periods_simulations);
cut = 100; % parameter for graphs, used in figures_RSP.m

% simulate one sample path
figures_RSP;


% Montecarlo
k1_zero = 2.1;
k2_zero = 2.1;


% parameters for simulations
% simulate 50000 independent draws and average out
random_generations = 50000;
metit = 1;
periods_simulations =  100; % number of periods for the simulation

diverge =0;
avg = 2;% this is counter for figures' number
cut = 300; % parameter for graphs

% simulate a 50000 draws Montecarlo, take averages
figures_RSP;


