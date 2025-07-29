%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 🌱 BATCH RUNNER FOR RICHARDS EQUATION SCENARIOS
%
% Description:
%   Loops through predefined .mat input scenarios for:
%   - Permeable Pavement (PP)
%   - Bioretention (BIO)
%   - Green Roof (GR)
%   - Pre-Development (PRE)
%
%   Executes solver, saves results, and runs post-processing.
%
% Author   : Adapted by ChatGPT
% Updated  : July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc;

%% === Add Relevant Paths ==================================================
addpath('Model_Configurations');
addpath('Numerical_Solver');
addpath('Physical_Functions_and_Utilities');
addpath('Main_Functions');
addpath('Visualization_and_Output');
addpath('Model_Calibration');
addpath('Forcing');
addpath('Examples');

%% === Define Scenario List ================================================
scenario_list = {
    % 🌧️ Permeable Pavement (PP)
    'PP_NWC.mat'
    'PP_MIA.mat'
    'PP_PHX.mat'
    'PP_SAN.mat'

    % 🌿 Bioretention
    % 'BIO_NWC.mat'
    % 'BIO_MIA.mat'
    % 'BIO_PHX.mat'
    % 'BIO_SAN.mat'

    % 🪴 Green Roof
    % 'GR_NWC.mat' % Already ran it
    % 'GR_MIA.mat'
    % 'GR_PHX.mat'
    % 'GR_SAN.mat'

    % 🌍 Pre-Development
    % 'PRE_NWC.mat'
    % 'PRE_MIA.mat'
    % 'PRE_PHX.mat'
    % 'PRE_SAN.mat'
};

%% === Loop Through Scenarios ==============================================
for i = 1:length(scenario_list)
    scenario_file = scenario_list{i};
    fprintf('\n🔄 Running scenario [%d/%d]: %s\n', i, length(scenario_list), scenario_file);

    % === Clear variables except persistent configs
    clearvars -except scenario_list i scenario_file

    % === Re-add paths after clear
    addpath('Model_Configurations');
    addpath('Numerical_Solver');
    addpath('Physical_Functions_and_Utilities');
    addpath('Main_Functions');
    addpath('Visualization_and_Output');
    addpath('Model_Calibration');
    addpath('Forcing');
    addpath('Examples');

    %% === Load Scenario
    load(fullfile('Examples', scenario_file));
    disp('📦 Scenario loaded.');

    %% === Set delta and Gamma
    if contains(scenario_file, 'PP_')
        Delta = 2;
        Gamma = 0.3;
    else
        Delta = 0;
        Gamma = 1;
    end

    %% === Solver Parameters (tolerances and timestep)
    params.tol = 1e-6;
    params.dt_min = 1e-6;
    params.dt_max = 24 * 3600;

    %% === Set Start Date
    start_datetime = datetime(1995, 1, 1, 0, 1, 0);

    %% === Run Main Solver
    run('Main_Solver.m');
    disp('✅ Solver completed.');

    %% === Save Results
    SaveSimulationResults;
    disp('💾 Results saved.');

    %% === Post-Processing
    Post_Processing;
    disp('📊 Post-processing complete.');
    close all
end

disp('🏁 All scenarios finished.');
