.. _example02:

Example 02 — Infiltration & Redistribution in Clay-Loam
=======================================================

Goal
----
Demonstrate infiltration into a finer-textured (clay-loam) medium and subsequent
**redistribution** once surface input stops. Compared with sand, expect a slower
wetting-front advance, higher near-surface water content for the same imposed flux,
and a more gradual drainage.

How to run
----------

Uncomment **exactly one** line in ``Mixed_Richards_Model_1D.m``::

   % load('Examples/Example1_Infiltration_Sand.mat');
   load('Examples/Example2_Clay_Loam_Soil.mat');   % enable this
   % load('Examples/Example3_CapillaryRise.mat');
   % load('Examples/Example4_TopNeumann_Sandy.mat');

Press **Run** (F5). Outputs are written to ``Modeling_Results/<sim_name>/``.

Setup (summary)
---------------
- **Column**: depth :math:`L \approx 1.0` m; refined mesh near the surface (see *Log.txt* for exact mesh).
- **Soil**: clay-loam with van Genuchten–Mualem parameters
  :math:`\{\theta_r,\theta_s,\alpha,n,\ell,K_s,S_s\}` set inside the example workspace
  (exact values printed in *Log.txt*).
- **Boundary conditions**:
  - **Top**: **Neumann** (imposed rain flux) for a period, then **zero flux** (redistribution phase).
  - **Bottom**: **free drainage** (:math:`\partial h/\partial z = 0`).
- **Initial condition**: uniform suction (dry start) or hydrostatic, as stored in the example.
- **ET/evaporation**: disabled in this example to isolate redistribution (unless noted in the workspace).

What to expect
--------------
- **Infiltration phase**: a **slower** wetting-front than sand; higher capillary
  retention leads to a more gradual :math:`\theta(z,t)` profile.
- **Redistribution phase**: after the surface flux goes to zero, gradients relax and
  :math:`h(z,t)` smooths; :math:`q(z,t)` trends toward small downward seepage; storage decreases
  slowly compared with the sandy case.
- **Mass balance**: reported each save step; relative error should remain small.

Outputs
-------
- **Profiles**: :math:`h(z,t)`, :math:`\theta(z,t)`, :math:`q(z,t)` saved under ``Figures/Profiles``.
- **Time series**: surface flux (rain/zero), bottom seepage, storage change, and mass-balance in
  ``Figures/TimeSeries``; CSV/MAT in ``Data``; run details in ``Log.txt``.

Troubleshooting
---------------
- If convergence slows near saturation, reduce the initial time step (``params.dt``) and/or raise
  ``params.max_iters`` (e.g., 15–20).
- Verify flux **units** when customizing: top Neumann flux must be in :math:`\mathrm{m\,s^{-1}}`
  inside the solver (spreadsheets in mm/h are converted upstream).
- If profiles look overly sharp/oscillatory, check mesh spacing near the surface and consider a
  slightly finer refinement.

See also
--------
- Back to :doc:`Examples index <index>`
- :doc:`../userguide/index`

.. _fig-ex02-profiles:

.. figure:: /static/examples/example02.png
   :alt: Example 02 profiles: h(z,t) and θ(z,t) for clay-loam during infiltration and redistribution
   :align: center
   :width: 100%

   **Example 02** — Infiltration (early) and redistribution (late) in clay-loam. Slower wetting-front advance
   and higher near-surface retention than sand. Replace this image with your generated figure.
