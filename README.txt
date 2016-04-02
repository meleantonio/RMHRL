%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Repeated Moral Hazard and Recursive Lagrangeans
%
% Matlab codes
%
% 	Antonio Mele, April 2016 (First version: June 2008)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

(Feel free to use these Matlab codes, with appropriate recognition of my work and citation of my paper. If you have any problem in running the codes, for questions, bugs reporting or suggestions, open an issue on GitHub or send me an email to meleantonio@gmail.com)

This project solves four different models of repeated moral hazard. These models are described in the 2011 version of the paper. 

The folder "Paper" includes two versions of the paper. The examples solved by the code are described in the 2011 version, section 4.2, while the published 2014 version only include the last example. 

The folder "Figures" contains codes for solving the examples included in section 4.2 of the 2011 version of the paper. It also generates all pictures of the paper. These codes can be used to solve each model separately, and as a base for applying the same algorithm to new problems.

The folder "Tables" contains codes for creating the tables included in the section 4.3 to test accuracy and computational speed.

The folder "Libraries" contains a few routines used to solve the models.




IMPORTANT PRELIMINARY STEPS (more details in the README.txt files in each folder):


The folder "Libraries" must be added to your Matlab path. It includes two libraries:

1) the most recent 64-bit version of the CompEcon Toolbox by Miranda and Fackler. If you are using a 32-bit machine, you may want to use another version. All version of CompEcon can be freely downloaded at http://www4.ncsu.edu/~pfackler/compecon/toolbox.html

2) A library by Michael Reiter, which has a very efficient version of the Broyden algorithm (function broydn.m). IF YOU HAVE YOUR OWN VERSION OF THE BROYDEN ALGORITHM OR YOU WANT TO USE ANOTHER NONLINEAR EQUATION SOLVER, YOU MUST MODIFY ALL "main_INDIC.m" SCRIPTS (look for the function "broydn" in the MAIN LOOP). 

This code comes with no guarantee. Use or modify it at your own risk. 