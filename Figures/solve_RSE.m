%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Code solving and generating figures for the risk sharing model with 
%  two-sided moral hazard in an endowment economy
%  solved with recursive Lagrangeans
%
%                       Antonio Mele
%                       March 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
tic;

global alpha betta sig  epsil omega1 omega2 nu nu2 gam1 output_y s1 s2
global  rounds_approx Order_vector


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


% parameters for simulations
random_generations = 1;               % random_generations multiplied by metit gives
metit = 1;                                           % the total number of draws used; just for 
                                                           % technical reasons, we store
                                                           % the draws on big matrices
                                                           % with no more than 1000 rows.     
                                                           % For simulating one series,
                                                           % set both to 1
periods_simulations =  200; % number of periods for the simulation


rounds_approx = 1;
Order_vector = 10;
mainRSE;
testing_max

figures_RSE;
