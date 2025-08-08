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
% Last Updated: AUG 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =========================================================================
% 🌱 MAIN RICHARDS MODEL LAUNCHER
% Description:
%   Master script to run the 1D mixed-form Richards equation solver using
%   modular configuration, numerical engine, and post-processing routines.
%
% Folders:
%   📁 Model_Configuration/           – Input setup and parameter loading
%   📁 Numerical_Solver/              – Time loop, Newton method, residuals
%   📁 Physical_Functions_and_Utilities/ – Soil functions, mesh, BCs
%   📁 Visualization_and_Output/      – Plotting and result export
%   📁 Examples/                      – .mat files with benchmark setups
%   📁 Modeling_Results/              – Output figures, logs, saved results
%
% Author   : Marcus Nóbrega, Ph.D.
% Updated  : May 2025
%% =========================================================================

clearvars; close all; clc;

%% === Add Relevant Paths ==================================================
addpath('Model_Configurations');               % Input setup & parameters
addpath('Numerical_Solver');                   % Solver and residual logic
addpath('Physical_Functions_and_Utilities');   % Mesh, θ-K, BC utilities
addpath('Main_Functions');                     % Main Solvers
addpath('Visualization_and_Output');           % Plotting, exports
addpath('Model_Calibration');                  % Files for model calibration
addpath('Forcing');                            % Catchment hydrologic functions
addpath('Examples');                           % Examples pre-defined

%% =========================================================================
% 📂 LOAD SCENARIO
% Choose one of the following options to load or configure your simulation.
% =========================================================================

% 🔹 Option 1: Load Parameters from a Saved .mat File (Pre-generated Workspaces)
% Uncomment ONLY one of the following as needed:
% load('Examples/Celia1990.mat');
% load('Examples/Example1_Infiltration_Sand.mat');
% load('Examples/Example2_Clay_Loam_Soil.mat');
% load('Examples/Example3_CapillaryRise.mat');
% load('Examples/Example4_TopNeumann_Sandy.mat');
% load('Examples/Monitored_PP_Data.mat');
% load('Examples/Monitored_PP_Events_Data.mat');

% 🔹 Option 2: Load City-Based Scenarios

% Permeable Pavement (PP)
% load('Examples/PP_NWC.mat');
% load('Examples/PP_MIA.mat');
% load('Examples/PP_PHX.mat');
% load('Examples/PP_SAN.mat');

% Bioretention
% load('Examples/BIO_NWC.mat');
% load('Examples/BIO_MIA.mat');
% load('Examples/BIO_PHX.mat');  
% load('Examples/BIO_SAN.mat');

% Green Roof
% load('Examples/GR_NWC.mat');
% load('Examples/GR_MIA.mat');
% load('Examples/GR_PHX.mat');
% load('Examples/GR_SAN.mat');

% Pre-development Cases
% load('Examples/PRE_NWC.mat');
% load('Examples/PRE_MIA.mat');
% load('Examples/PRE_PHX.mat');
% load('Examples/PRE_SAN.mat');

% 🔹 Option 3: Generate Scenario from Input Script
% Customize your own scenario by modifying 'Input_Data.m'
% Make sure to read the instructions inside the script

run('Input_Data.m');  % Default configurable script

% ✅ Confirmation
disp('📦 Input parameters loaded.');

%% =========================================================================
% 🚀 EXECUTE SOLVER
% =========================================================================

% 🔹 Simulation Start Time (Used for time-series plotting)
% Format: datetime(YYYY, MM, DD, HH, MM, SS)
start_datetime = datetime(2023, 12, 21, 0, 10, 0);

% 🔹 Empirical Evaporation / Evapotranspiration Coefficients
% These affect the surface water balance.
% Reference values (see manuscript draft for details):
%   - Delta: 0.3 for Permeable Pavement (PP), 1 for others
%   - Gamma: 2 for Permeable Pavement (PP), 0 for others

Delta = 1;   % General default
Gamma = 0;   % General default

%% 🔹 Run Main Solver Script
run('Main_Solver.m');

%% 🔹 Save Results ========================================================
SaveSimulationResults

%% 🔹 Run Post-Processing Visualization ===================================
Post_Processing;
