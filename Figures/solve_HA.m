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

% parameters for simulations
random_generations = 1000; % 50000;%            % random_generations multiplied by metit gives
metit = 1;                                           % the total number of draws used; just for
% technical reasons, we store
% the draws on big matrices
% with no more than 1000 rows.
% For simulating one series,
% set both to 1
periods_simulations =  70; % number of periods for the simulation



rounds_approx = 3;
Order_vector = [[4;4] [6;6]  [10; 10]  [6;4]  [8;4] [10; 4] [20;20]];%[4;4] [5;5]  [8;8] 
nClust = 16; % # clusters for choosing the grid

% Create grid extrema for Pareto weights
phibar = gam1;
phi_sd =.3*gam1;
phi_min = phibar-phi_sd; % .3;%
phi_max = phibar+phi_sd;% .8; 2.5*.25  1.3;%

% Create grid extrema for costate EE
zetabar= .422025;%1;%
zeta_min = 0;
zeta_max = zetabar;


mainHA;
testing_max
verification_procedure;
phizero = phibar;%1;%.7;
figures_HA;
testing_grid;
