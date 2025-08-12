Preparing your input data
----------------------------------------

This section explains how to configure a run by editing **Input_Data.m**.
The driver script must use your Input_Data file (and **not** a pre-packaged scenario).

.. warning::

   In the main driver (e.g., ``Mixed_Richards_Model_1D.m``), make sure a scenario loader is **not** used.
   Use your input file instead:

   .. code-block:: matlab

      % --- Load user inputs
      clear params;
      run('Input_Data.m');     % ✔️ use this

      % load('Scenarios/Example01.mat')  % ❌ comment/remove any scenario loads

1) Simulation name and output folders
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set a unique name for each run; folders are auto-created under ``Modeling_Results/<sim_name>/``:

.. code-block:: matlab

   sim_name = 'Generic_Example';   % 🔁 change per run

   base_output_dir   = fullfile('Modeling_Results', sim_name);
   figures_dir       = fullfile(base_output_dir, 'Figures');
   data_dir          = fullfile(base_output_dir, 'Data');
   mesh_dir          = fullfile(figures_dir, 'Mesh');
   profiles_dir      = fullfile(figures_dir, 'Profiles');
   time_series_dir   = fullfile(figures_dir, 'TimeSeries');
   diagnostics_dir   = fullfile(figures_dir, 'Diagnostics');
   fdc_dir           = fullfile(figures_dir, 'FlowDuration');
   wetting_dir       = fullfile(figures_dir, 'WettingFront');

Tips:

- Use only letters, numbers, underscores in ``sim_name``.
- All folders are created if missing; you do not need to create them manually.

2) Domain and mesh
~~~~~~~~~~~~~~~~~~

Choose total depth and number of vertical nodes; an optional refinement near the surface is supported.

.. code-block:: matlab

   params.Nz = 27;            % nodes [-]
   params.L  = 1.0;           % total depth [m], z positive upward
   nonlin_factor = 1;         % 1 = uniform mesh (increase for surface refinement)
   params.LID_area = 1.0;     % 1D column area [m^2]
   [params.z, params.dz] = generate_nonlinear_mesh(params.Nz, params.L, nonlin_factor, mesh_dir);

3) Time stepping and solver controls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: matlab

   params.Tmax   = 24*3600;   % total time [s]
   params.dt     = 1;         % initial dt [s]
   params.dt_min = 1e-3;      % min dt [s]
   params.dt_max = 5*60;      % max dt [s]

   params.adapt_down = 0.5;   % shrink factor
   params.adapt_up   = 2.0;   % growth factor
   params.n_up       = 5;     % threshold for fast convergence
   params.n_down     = 10;    % threshold for slow convergence

   params.max_iters  = 20;    % Newton max iterations
   params.tol        = 1e-6;  % convergence (head) [m]

   % Output cadence
   params.save_interval_min = 5;  % [min]
   params.save_interval     = max(params.save_interval_min*60, params.dt); % [s]

4) Layers and hydraulic parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Specify thickness (bottom→top) and van Genuchten–Mualem parameters per layer.
**Units**: ``alpha [m^-1]``, ``n [-]``, ``theta_r, theta_s [-]``, ``S_s [m^-1]``, ``Ks [m s^-1]``.

.. code-block:: matlab

   media_thicknesses = [1.0];   % bottom->top [m]; add more values for multilayer
   media_props = struct( ...
       'alpha',   [14.5], ...
       'n',       [2.68], ...
       'theta_r', [0.045], ...
       'theta_s', [0.43], ...
       'S_s',     [1e-5], ...
       'Ks',      [8.25e-5]);

   % Map layers to nodes (auto)
   media_interfaces = [-params.L + cumsum(media_thicknesses)];
   media_interfaces = [-params.L, media_interfaces];
   n_layers = numel(media_thicknesses);
   media_id = zeros(1, params.Nz);
   for i = 1:params.Nz
       zi = params.z(i);
       for j = 1:n_layers
           if zi >= media_interfaces(j) && zi < media_interfaces(j+1)
               media_id(i) = j; break;
           elseif zi == media_interfaces(end)
               media_id(i) = n_layers;
           end
       end
   end

   params.alpha   = media_props.alpha(media_id);
   params.n       = media_props.n(media_id);
   params.m       = 1 - 1 ./ params.n;
   params.theta_r = media_props.theta_r(media_id);
   params.theta_s = media_props.theta_s(media_id);
   params.S_s     = media_props.S_s(media_id);
   params.Ks      = media_props.Ks(media_id);

5) Boundary conditions
~~~~~~~~~~~~~~~~~~~~~~

Top boundary (Dirichlet head)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Specify a fixed pressure head at the surface:

.. code-block:: matlab

   params.top_bc_type  = "dirichlet";
   params.top_bc_value = -0.10;         % head at top [m]; h>=0 allows ponding

Bottom boundary
^^^^^^^^^^^^^^^
Choose one: free drainage, fixed head, specified flux, or no-flow.

.. code-block:: matlab

   params.bottom_bc_type  = "free";     % 'free' | 'dirichlet' | 'neumann' | 'noflow'
   params.bottom_bc_value = 0;          % [m] (ignored for 'free' and 'noflow')

Top boundary (Neumann flux) with forcing file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When using a flux boundary at the surface, provide a forcing spreadsheet
``Catchment_Forcing.xlsx`` and point the solver to it. The sheet must contain
**Time**, **Rainfall**, **ET0**, and an optional **Inflow Hydrograph**, plus a
**flag** that controls whether a catchment rainfall–runoff model is used.

**Required columns and units**
- ``Time (min)`` — minutes from start of simulation.
- ``Rainfall (mm/h)`` — precipitation intensity.
- ``ET0 (mm/h)`` — reference evapotranspiration.
- ``Inflow Hydrograph (mm/h)`` — *optional* exogenous inflow to the LID (see flag below).
- ``flag_manual_inflow`` — scalar flag (cell on the sheet, value 0 or 1).

.. figure:: /static/Catchment_Forcing_Fig.png
   :alt: Example Catchment_Forcing.xlsx layout with Time, Rainfall, ET0, Inflow Hydrograph, and flag_manual_inflow
   :align: center
   :width: 700px

   Example layout for ``Catchment_Forcing.xlsx``. Place the file in your ``Forcing/`` folder.

**How the flag works**

- ``flag_manual_inflow = 1`` → **Manual net rainfall only**.
  The model **does not** run the catchment rainfall–runoff module.
  The surface flux is computed as

  .. math::

     q_{\text{top}}(t) \;=\; \big(\text{Rainfall}(t) - ET_0(t)\big)\;\;[\mathrm{mm/h}]

  and internally converted to :math:`[\mathrm{m\,s^{-1}}]`.
  The “Inflow Hydrograph (mm/h)” columiltration into the columnn is always considered as - Rainfall or Precipitation + Evapotranspiration (positive upwards)

- ``flag_manual_inflow = 0`` → **Use upstream catchment model**.
  The function ``Catchment_Outputs.m`` reads the spreadsheet, computes the upstream
  **runoff hydrograph** from the catchment parameters you configure there (area, CN,
  :math:`\lambda`, structures, etc.), and returns the net surface flux to apply at the LID top:

  .. math::

     q_{\text{top}}(t) \;=\; \big(\text{Rainfall}(t) - ET_0(t)\big) \;+\; q_{\text{upstream}}(t)

  in consistent units. The “Inflow Hydrograph (mm/h)” column in the spreadsheet is **not used** in this mode.

**Hooking the file into the model**

.. code-block:: matlab

   params.top_bc_type = "neumann";                        % flux BC at surface
   forcing_path = fullfile('Forcing','Catchment_Forcing.xlsx');

   % LID area is used to normalize/runoff per unit area if needed
   [surface_flux_time, surface_flux_vals, C_top] = ...
       Catchment_Outputs(params.dt, params.LID_area, forcing_path);

   % Store for interpolation in the solver (units already converted inside Catchment_Outputs)
   params.surface_flux_time = surface_flux_time;          % [s]
   params.surface_flux_vals = surface_flux_vals;          % [m/s]
   % C_top is optional tracer concentration [M/L^3]


**Units and sign convention**

- Spreadsheet inputs are in **mm/h**; the function converts to **m/s** internally.
  Conversion: :math:`q[\mathrm{m\,s^{-1}}] = (\text{mm/h} \times 10^{-3})/3600`.
- Positive :math:`q_{\text{top}}` is **upward**.
  ET increases the net flux (i.e., acts upward/positive).
- Time is converted from minutes to seconds internally.

**Common mistakes**
- The spreadsheet name/path is wrong → set ``forcing_path`` correctly.
- ``flag_manual_inflow`` not set → default mode may surprise you.
- Mixed units → keep Rainfall/ET0/Inflow in **mm/h** only.
- Forgot to choose ``"neumann"`` at the top boundary → the file is ignored if you are in Dirichlet mode.

6) Input data defined inside ``Catchment_Outputs.m`` (used only for top Neumann BC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``params.top_bc_type == "neumann"``, the model calls
``Catchment_Outputs(dt, LID_area, forcing_path)`` to build the surface flux
time series for the LID. **This function contains site/case parameters that are
hard-coded at the top of the file and must be edited for each test case.**

Where it is used
^^^^^^^^^^^^^^^^
- Load path in your input script:
  ::
     forcing_path = fullfile('Forcing','Catchment_Forcing.xlsx');
     [t_forc_s, qtop_ms, C_top] = Catchment_Outputs(params.dt, params.LID_area, forcing_path);
     params.surface_flux_time = t_forc_s;    % [s]
     params.surface_flux_vals = qtop_ms;     % [m/s]
- If the top boundary is **Dirichlet**, this function is **not used**.

Hard-coded watershed parameters (edit per case)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Inside ``Catchment_Outputs.m`` you will see:

.. code-block:: matlab

   %% 1) Watershed physical properties
   width_w      = 10;      % [m]   watershed width
   length_w     = 50;      % [m]   watershed length
   Watershed_Area = width_w * length_w;  % [m^2]

   %% 2) SCS Curve Number setup
   CN_per       = 65;      % [-]   pervious CN
   h0_w         = 0.006;   % [m]   initial abstraction
   n_w          = 0.02;    % [-]   Manning roughness
   Aimp         = 0.8;     % [-]   impervious fraction
   slope_w      = 0.015;   % [-]   bed slope (m/m)
   baseflow_w   = 0;       % [m^3/s] baseflow
   ks_w         = 10;      % [mm/h] soil K for recovery
   kr_w         = 2;       % [mm/h] recovery rate
   A_GI         = 0;       % [-]   green-infra area ratio
   CN_GI        = 65;      % [-]   GI CN
   catch_GI_imp = 0;       % [-]   fraction of impervious draining to GI

**These values are entered directly in the function. Change them to represent your
watershed every time you run a new case.** (You may later refactor them into your
``Input_Data.m`` if you prefer not to edit the function.)


Optional pollutant model (inside the function)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``Catchment_Outputs.m`` can also compute a surface concentration using a
buildup–washoff model. Parameters appear at the end of the function:

.. code-block:: matlab

   % Buildup–washoff parameters (edit per case)
   C1 = 50;      % [kg/ha]
   C2 = 0.3;     % [1/day]
   C3 = 0.02;    % [(mm/h)^(-C4) * h^-1]
   C4 = 1.5;     % [-]
   ADD = 10;     % antecedent dry days [day]

   Area_km2 = Watershed_Area / 1e6;     % [km^2]
   [~, Cpol, ~, ~, ~, ~, ~, ~] = Buildup_Washoff_Model(C1,C2,C3,C4,ADD,Area_km2,Qin,time_watershed);

.. note::

   The watershed parameters and pollutant coefficients above are **hard-coded defaults**
   to keep the function self-contained. If your workflow involves many runs, consider
   passing these values from ``Input_Data.m`` to avoid editing the function each time.


7) Optional drainage structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Orifices and a top spillway can be enabled per node:

.. code-block:: matlab

   % Orifice (per node): Q = K_orifice * (max(h,0))^exp_orifice
   params.K_orifice   = zeros(1, params.Nz);   % [m^(exp)/s]
   params.exp_orifice = 0.5 * ones(1, params.Nz);

   % Example coefficients (disabled by default)
   node_idx   = 20; n_orifices = 0; Cd = 0.6; D = 0.10; g = 9.81;
   Aeff = pi*D.^2/4;
   params.K_orifice(node_idx) = n_orifices * Cd .* Aeff * sqrt(2*g);

   % Spillway at top (Neumann BC only): Q = c_spillway * (h - h_spill)^(exp_spillway)
   params.spillway_enabled = true;
   params.c_spillway       = 0 * 1.8 * 1.5;   % set >0 to enable
   params.h_spill          = 0.05;            % [m]
   params.exp_spillway     = 1.5;

8) Initial condition and sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: matlab

   % Hydrostatic/constant suction profile
   p = -3;                       % [m]
   h = p * ones(1, params.Nz);   % here you can enter a veritical profile of pressure head [m]
   if params.bottom_bc_type == "dirichlet"; h(1)  = params.bottom_bc_value; end
   if params.top_bc_type    == "dirichlet"; h(end)= params.top_bc_value;   end

   % Optional distributed source term S(z,t) [s^-1] (default zero)
   params.source_times   = linspace(0, params.Tmax, params.Nt);
   params.source_profile = zeros(params.Nz, params.Nt);

10) Run and outputs
~~~~~~~~~~~~~~~~~~~

Run the driver; outputs go to ``Modeling_Results/<sim_name>/`` and include figures, CSV/MAT files,
and a detailed log (``Log.txt``) with mesh, parameters, BCs, solver settings, and mass-balance summary.
