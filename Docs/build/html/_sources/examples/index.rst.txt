.. _examples:

========
Examples
========

DRAIN-LID ships with ready-to-run **workspace files** (``.mat``) for the paper’s cases and several city scenarios.
To run an example, you only need to **uncomment one `load(...)` line** in the driver and press **Run**.

.. toctree::
   :maxdepth: 1

   example01
   example02
   example03
   example04
   example05

How to run an example
---------------------
1. In MATLAB:
   ::
      addpath(genpath('path/to/LID_Tool'))

2. Open ``Mixed_Richards_Model_1D.m`` and find the **LOAD SCENARIO** block.
   Uncomment **exactly one** line and keep the others commented:
   ::
      % 🔹 Option 1: Load Parameters from a Saved .mat File
      % load('Examples/Celia1990.mat');
      load('Examples/Example1_Infiltration_Sand.mat');   % ← this one enabled
      % load('Examples/Example2_Clay_Loam_Soil.mat');
      % load('Examples/Example3_CapillaryRise.mat');
      % load('Examples/Example4_TopNeumann_Sandy.mat');
      % load('Examples/Monitored_PP_Data.mat');
      % load('Examples/Monitored_PP_Events_Data.mat');

   (You can also try the *city-based* scenarios under **Option 2** in the same block.)

3. Press **Run** (F5). The solver will execute using the parameters contained in the selected ``.mat`` file.

4. Results (figures, data, logs) will be written under ``Modeling_Results/<sim_name>/``.

Available core examples
-----------------------
These correspond to the paper’s four cases (names match the ``.mat`` files):

- **Example1_Infiltration_Sand.mat** — Constant-flux infiltration in sandy media; wetting-front advance; free drainage at bottom.
- **Example2_Clay_Loam_Soil.mat** — Infiltration/redistribution in finer media; sensitivity to retention curve.
- **Example3_CapillaryRise.mat** — Capillary rise from a fixed water table (bottom Dirichlet), approach to steady state.
- **Example4_TopNeumann_Sandy.mat** — Constant top rainfall flux
- **Example5_Continous_PP.mat** — Time-varying rainfall/ET with ponding switch (flux→head) in a real-world permeable pavement.

Validation & monitored-data cases
---------------------------------
- **Celia1990.mat** — Classic benchmark for Richards’ equation (transient infiltration test).
- **Monitored_PP_Data.mat / Monitored_PP_Events_Data.mat** — Permeable pavement runs calibrated to observed events.

City-based scenarios
--------------------
Uncomment one under **Option 2** to load representative climates for each LID type:

- **Permeable Pavement (PP)**: ``PP_NWC.mat``, ``PP_MIA.mat``, ``PP_PHX.mat``, ``PP_SAN.mat``
- **Bioretention (BIO)**: ``BIO_NWC.mat``, ``BIO_MIA.mat``, ``BIO_PHX.mat``, ``BIO_SAN.mat``
- **Green Roof (GR)**: ``GR_NWC.mat``, ``GR_MIA.mat``, ``GR_PHX.mat``, ``GR_SAN.mat``
- **Pre-development (PRE)**: ``PRE_NWC.mat``, ``PRE_MIA.mat``, ``PRE_PHX.mat``, ``PRE_SAN.mat``

Using your own configuration (optional)
---------------------------------------
If you want to build a scenario from scratch instead of using a ``.mat`` file:

1. Open ``Input_Data.m`` and edit layers, BCs, solver controls, and (if using a top **Neumann** BC) the **forcing** path.
2. In ``Mixed_Richards_Model_1D.m``, **comment out** all ``load(...)`` lines and **uncomment**:
   ::
      % 🔹 Option 3: Generate Scenario from Input Script
      run('Input_Data.m');

3. Run the model. See :doc:`../userguide/index` for the Input_Data fields and the
   :doc:`../userguide/parameters` section for the Neumann forcing spreadsheet format.

Notes & troubleshooting
-----------------------
- **Only one** scenario should be active (one ``load(...)`` or the ``run('Input_Data.m')`` line).
- Some Neumann-flux examples read ``Forcing/Catchment_Forcing.xlsx`` via
  ``Catchment_Outputs.m``. Ensure the spreadsheet exists and columns/units match the template
  (Time [min], Rainfall [mm/h], ET0 [mm/h], Inflow [mm/h], and ``flag_manual_inflow``).
- If Newton does not converge, reduce the initial time step (``params.dt``) or increase
  ``params.max_iters``; verify units of BCs and layer parameters.
