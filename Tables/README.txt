%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tables for section "Computational speed and accuracy"
%
% Antonio Mele, April 2016 (First Version: June 2008)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This folder collects all files needed to generate the tables included in the paper "Repeated Moral Hazard and Recursive Lagrangeans", 2011 version, section 4.3. These tables check the computational speed of the algorithm for different sizes of the grid. 

The main file to run is generate_tables.m. This file generates all tables for all models, by calling four different scripts called "tables_INDIC.m" where INDIC = {RMH, HA, RSE, RSP}. 

Each tables_INDIC.m assigns values to parameters, and solves the model for different grids. In the main loop, it calls a file "main_INDIC.m". The script "main_INDIC.m" creates the grid and the functional approximations (using Miranda-Fackler's Compecon toolbox, which MUST be installed in your pc),  and solves the first-order conditions of the Lagrangean (in the file INDIC_focs.m). I use a version of the Broyden algorithm coded by Michael Reiter, included in the "Libraries" folder. IF YOU HAVE YOUR OWN VERSION OF THE BROYDEN ALGORITHM OR YOU WANT TO USE ANOTHER NONLINEAR EQUATION SOLVER, YOU MUST MODIFY ALL "main_INDIC.m" SCRIPTS (look for the function "broydn" in the MAIN LOOP). 

The file tables_INDIC stores the table in a txt file tables_final_INDIC.txt. 