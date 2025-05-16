%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 🌱 RICHARDS EQUATION SOLVER (1D Vertical) — MULTI-LAYER SOIL PROFILE
%
% Description:
%   Solves the transient mixed-form Richards Equation for unsaturated flow
%   in a variably saturated soil column with multiple soil types.
%
% Features:
%   ✅ Nonlinear vertical grid refinement near the surface (Hydrus-style)
%   ✅ Spatially variable soil hydraulic properties (Van Genuchten + Ks)
%   ✅ Hydrostatic initial condition
%   ✅ Dirichlet, Neumann, Free, and No-Flow boundary conditions
%   ✅ Outputs pressure head, moisture, and bottom fluxes over time
%
% Developed by: Marcus Nobrega, Ph.D. 🌎
% Last Updated: May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =========================================================================
% 🌱 MAIN RICHARDS MODEL LAUNCHER
% Description:
%   Master script to run the 1D mixed-form Richards equation solver using
%   modular configuration, numerical engine, and post-processing routines.
%
% Folders:
%   📁 model_configuration/           – Input setup and parameter loading
%   📁 numerical_solver/             – Time loop, Newton method, residuals
%   📁 physical_functions_and_utilities/ – Soil functions, mesh, BCs
%   📁 visualization_and_output/     – Plotting and result export
%   📁 test_cases_and_scenarios/     – .mat files with benchmark setups
%   📁 simulation_outputs/           – Output figures, logs, saved results
%
% Author   : Marcus Nóbrega, Ph.D.
% Updated  : May 2025
%% =========================================================================

clearvars; close all; clc;

%% === Add Relevant Paths ==================================================
addpath('model_configuration');                % Input setup & parameters
addpath('numerical_solver');                   % Solver and residual logic
addpath('physical_functions_and_utilities');   % Mesh, θ-K, BC utilities
addpath('visualization_and_output');           % Plotting, exports
addpath('examples');                           % Examples pre-defined


%% === Load Example Case ===================================================
% Options (Uncomment one to test):
% load('test_cases_and_scenarios/Example1_Infiltration.mat');
% load('test_cases_and_scenarios/Celia1990.mat');
load('test_cases_and_scenarios/PP.mat');  % ← Currently active case

%% === Execute Solver ======================================================
Main_Solver_no_line_search;

%% === Run Post-Processing Visualization ===================================
Post_Processing;
