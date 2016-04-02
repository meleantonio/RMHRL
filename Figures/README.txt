%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Repeated Moral Hazard and Recursive Lagrangeans
%
% Codes for solving numerical examples in
% section 5 ("Numerical examples") of the 2011 version
% Published version (http://dx.doi.org/10.1016/j.jedc.2014.03.007)
% includes only the last example (with suffix RSP, see below)
%
% 			Antonio Mele, April 2016 (First version: June 2008)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


This folder collects all Matlab files needed to solve the examples included in the paper "Repeated Moral Hazard and Recursive Lagrangeans" (2011 version), and plots all figures included in the paper. Every script has a suffix _INDIC (where INDIC = {RMH, HA, RSE, RSP}) that identifies which example solves. In this way, each examples can be solved separately and the files can be modified for your own purposes. 


The meaning of the suffix:

RMH: standard repeated moral hazard
HA : repeated moral hazard with hidden assets
RSE: dynamic risk sharing with two-sided moral hazard in an endowment economy
RSP: dynamic risk sharing with two-sided moral hazard in a production economy with individual capital accumulation



PRELIMINARY STEPS:

The folder "Libraries" must be added to your Matlab path. It includes two libraries:

1) the most recent 64-bit version of the CompEcon Toolbox by Miranda and Fackler. If you are using a 32-bit machine, you may want to use another version. All version of CompEcon can be freely downloaded at http://www4.ncsu.edu/~pfackler/compecon/toolbox.html

2) A library by Michael Reiter, which has a very efficient version of the Broyden algorithm (function broydn.m). IF YOU HAVE YOUR OWN VERSION OF THE BROYDEN ALGORITHM OR YOU WANT TO USE ANOTHER NONLINEAR EQUATION SOLVER, YOU MUST MODIFY ALL "main_INDIC.m" SCRIPTS (look for the function "broydn" in the MAIN LOOP). 


To solve all examples, run "solve_all_examples.m". This file calls four different scripts called "solve_INDIC.m" where INDIC = {RMH, HA, RSE, RSP}. If you want to solve only one of the examples, look for the solve_INDIC.m file, and run that one. 

Each solve_INDIC.m assigns values to parameters, solves the model calling the file "main_INDIC.m", and generates the figures in the file "figures_INDIC.m". 

The script "main_INDIC.m"creates the grid and the functional approximations (using Miranda-Fackler's Compecon toolbox, which MUST be installed in your pc), and solves the first-order conditions of the Lagrangean (in the file "INDIC_focs.m") with Broyden algorithm.

In the HA example, before generating figures the code must verify if the first-order approach is valid. It does this with the routine "verification_procedure.m". YOU MAY WANT TO CHANGE PARAMETERS IN THIS ROUTINE TO GET A MORE PRECISE VERIFICATION. In particular, you may want to increase the number of grid points used.

The script "figures_INDIC.m" calls two files (except for the HA example and the RSP example, where the first is not called). The first file is "plot_grid_results_INDIC.m", which generates policy functions over a grid (This is done for HA example directly in the file "figures_INDIC.m", while it is not possible to do it for the RSP example). The second is "simul_INDIC.m", which simulates 50000 independent draws of each model for a horizon of 200 periods. 

For questions, bugs reporting or suggestions, open an issue on GitHub or send me an email to meleantonio@gmail.com
