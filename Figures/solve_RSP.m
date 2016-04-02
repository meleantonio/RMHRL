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




rounds_approx = 1;
Order_vector =  [ 4  ;   4   ; 4   ];

%     rounds_approx = 3;
%     Order_vector =  [ 4  6 8 ;   4  6  8 ; 4  6 8 ];

mainRSP;
testing_max

avg = 0; % this is counter for figures' number

figures_RSP;

% parameters for simulations
% simulate 50000 independent draws and average out
random_generations = 1000;
metit = 50;
periods_simulations =  200; % number of periods for the simulation

avg = 2;% this is counter for figures' number
figures_RSP;
