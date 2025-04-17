%% EXAMPLE GENERIC CASPER HEADER FOR NEW SCRIPTS/FUNCTIONS (your title should replace this line)
%%% Brief description: %%%
% Calculates the magnetic flux and geometric flling factor in a pick-up coil (solenoid or second order
% gradiometer) with user-provided geometry and wire turns based on a
% publication: https://pubmed.ncbi.nlm.nih.gov/19191458/ (Minimizing sample shape and size effects in all purpose magnetometers)
%%% Author: Andrew J. Winter, Boston University %%%
%%% Code originally written: 12/19/2023 %%%
%%% Inputs: %%%
% - (example) Gradiometer structure containing:
% --- Number of turns (unitless)
% --- Winding separation distance (mm)
% --- Radius of coil (mm)
%%% Outputs: %%%
% - (example) Magnetic flux in fT/sqrt(Hz) and filling factor (unitless) as a function of turn number
%%% Revisions: (GITHUB SHOULD TAKE CARE OF THIS BUT MAYBE GOOD PRACTICE TO KEEP TRACK HERE TOO) %%% 
% 7/16/2024
% - Andrew made this a standalone function instead of an ugly script
% - Input is now a single structure that contains all the necessary input values
% - Added support for double-wound coils (since each loop-pair is approximately same radius from center of the sample, the function just evaluates for
%   one iteration of the loops and doubles the results due to superposition principle / symmetry)
%%% Notes: %%%
% - Unit tests performed to make sure calculations are valid, as described in LaTeX write-up
% - Obviously this can be optimized a little better, but considering the calculations this runs pretty fast in most cases (< 2 min for N = 16  gradiometer)
%%% MATLAB dependencies: %%%
%%% Code written with MATLAB 2022a %%%
% REQUIRED MATLAB TOOLBOXES: 
% - Symbolic Math (for legendreP function)
% SUGGESTED MATLAB TOOLBOXES:
% - Parallel Computing (for parfor loops, speeds up code significantly)





















