.. _example05:

Example 05 — Observed Permeable Pavement (255 days)
===================================================

Goal
----
Simulate an **observed permeable pavement** column continuously for **255 days** and
analyze time series of **pressure head** :math:`h`, **water content** :math:`\theta`,
and **Darcy flux** :math:`q`. This example uses monitored rainfall/ET forcing and
site-specific layer parameters included in the workspace.

How to run
----------
Uncomment **exactly one** line in ``Mixed_Richards_Model_1D.m``:

.. code-block:: matlab

   % load('Examples/Example4_TopNeumann_Sandy.mat');
   load('Examples/Monitored_PP_Data.mat');       % ← enable this
   % load('Examples/Monitored_PP_Events_Data.mat');

Press **Run** (F5). Outputs are written to ``Modeling_Results/<sim_name>/`` (Figures, Data, Log.txt).

Setup (summary)
---------------
- **Duration**: ~255 days (continuous, variable rain/ET).
- **Top boundary**: **Neumann** (observed rain/ET forcing).
  The example workspace provides the time series already converted for the solver; if you
  switch to your own data, keep spreadsheet inputs in **mm/h** and convert to **m/s**
  upstream (see User Guide → Boundary conditions).
- **Bottom boundary**: typically **free drainage** (:math:`\partial h/\partial z = 0`) unless
  otherwise specified in the workspace.
- **Layers**: pavement system (e.g., surface + base/subbase), with van Genuchten–Mualem
  parameters (:math:`\theta_r,\theta_s,\alpha,n,\ell,K_s`) and **specific storage** :math:`S_s`
  set in the workspace.
- **Evaporation**: Hargreaves :math:`ET_0` with Feddes reduction; pavement option
  :math:`E_0=\gamma(ET_0+\delta)` may be enabled (see run log).
- **Structures (optional)**: spillway/orifice disabled by default unless enabled in the workspace.

Outputs
-------
- **Time series** (``Figures/TimeSeries``):
  - Top flux components (rain, ET), **net surface flux**, **ponding depth** (if any), **runoff** (if enabled)
  - **Seepage** at the bottom, storage change, cumulative **mass balance**
- **Profiles** (``Figures/Profiles``): :math:`h(z,t)`, :math:`\theta(z,t)`, :math:`q(z,t)` at selected dates
- **Data** (``Data/``): CSV/MAT exports of the variables above
- **Log** (``Log.txt``): mesh, layers/BCs, solver tolerances, mass-balance summary

Notes & tips
------------
- **Performance for long runs**: consider coarser save cadence (e.g., hourly) to keep figure/data
  counts manageable; this is usually pre-configured in the workspace.
- **Units**: spreadsheet inputs (if you swap them) should be **mm/h**; solver uses **m/s**.
- **Pavement evaporation**: if :math:`E_0=\gamma(ET_0+\delta)` is active, confirm the chosen
  :math:`\gamma` and :math:`\delta` in your workspace match the site calibration.

Troubleshooting
---------------
- *No surface response during storms*: verify top BC is **Neumann** and the forcing series
  is nonzero (and correctly converted to m/s).
- *Slow convergence in wet periods*: reduce initial ``params.dt`` or increase ``params.max_iters``;
  check layer :math:`K_s` and retention parameters.
- *Large mass-balance error*: inspect the time step control (tighten tolerances or reduce
  ``dt_max`` during heavy events).

See also
--------
- Back to :doc:`Examples index <index>`
- :doc:`../userguide/index`
- :doc:`example05` (constant-flux Neumann benchmark)

.. _fig-ex05-timeseries:

.. figure:: /static/examples/example05.png
   :alt: Example 05: Permeable pavement—time series of head, θ, and q over 255 days
   :align: center
   :width: 100%

   **Example 05** — Observed PP over 255 days. Example time series of surface forcing, ponding/runoff,
   and profiles/flux summaries. Replace with your generated figure saved at
   ``_static/examples/TimeSeries/ex05_pp_timeseries.png``.
