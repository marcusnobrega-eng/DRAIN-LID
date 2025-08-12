Quickstart
==========

Goal
----
Run a DRAIN-LID example and produce θ, h, q profiles.

Prerequisites
-------------
- MATLAB (R2021a+)
- Repo path added to MATLAB (see :doc:`../install/index`)

Steps
-----
1. In MATLAB, add the repo to your path (if not already):
   ::
      addpath(genpath('path/to/LID_Tool'))

2. Open the Main Model File at your main path
   ::
      Mixed_Richards_Model_1D.m


3. Uncomment the example 1 line
   ::
      load('Examples/Example1_Infiltration_Sand.mat');
4. Make sure to comment all other lines that input pre-defined examples, such that your loading scenarion is like:
   ::
        %% =========================================================================
        % LOAD SCENARIO
        % Choose one of the following options to load or configure your simulation.
        % =========================================================================

        % Option 1: Load Parameters from a Saved .mat File (Pre-generated Workspaces)
        % Uncomment ONLY one of the following as needed:
        % load('Examples/Celia1990.mat');
        load('Examples/Example1_Infiltration_Sand.mat');
        % load('Examples/Example2_Clay_Loam_Soil.mat');
        % load('Examples/Example3_CapillaryRise.mat');
        % load('Examples/Example4_TopNeumann_Sandy.mat');
        % load('Examples/Monitored_PP_Data.mat');
        % load('Examples/Monitored_PP_Events_Data.mat');

        % Option 2: Load City-Based Scenarios

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

        % Option 3: Generate Scenario from Input Script
        % Customize your own scenario by modifying 'Input_Data.m'
        % Make sure to read the instructions inside the script

        % run('Input_Data.m');  % Default configurable script

        % Confirmation
        disp('Input parameters loaded.');

5. Note that you can run any scenario by uncommenting the line where it is logged. If two or more loading scenarios are chosen, the last one will be run.

6. If you want to enter personalized data to your case, you need to modify the function:
   ::
      Model_Configurations/Input_Data.m

7. This function is described further in the UserGuide section.

Results of this scenario are summarized by water content :math:`{\theta}` and pressure :math:`{h}` values

