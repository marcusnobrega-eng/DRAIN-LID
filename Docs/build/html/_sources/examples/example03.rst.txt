.. _example03:

Example 03 — Capillary Rise from a Water Table
==============================================

Goal
----
Demonstrate **upward wetting** driven by a fixed water table at the column base,
and the transient approach to steady state. This case highlights the retention
effects in :math:`\theta(h)` and the evolution of :math:`h(z,t)` under **zero surface flux**.

How to run
----------
Uncomment **exactly one** line in ``Mixed_Richards_Model_1D.m``:

.. code-block:: matlab

   % load('Examples/Example1_Infiltration_Sand.mat');
   % load('Examples/Example2_Clay_Loam_Soil.mat');
   load('Examples/Example3_CapillaryRise.mat');   % ← enable this
   % load('Examples/Example4_TopNeumann_Sandy.mat');

Press **Run** (F5). Outputs are written to ``Modeling_Results/<sim_name>/`` (Figures, Data, Log.txt).

Setup (summary)
---------------
- **Column:** depth :math:`L \approx 1.0` m with refined mesh near the surface/bottom (see *Log.txt* for exact mesh).
- **Soil:** van Genuchten–Mualem parameters loaded from the example workspace (exact values printed in *Log.txt*).
- **Boundary conditions:**
  - **Top:** **zero flux** (no rain; ET disabled) so that upward capillarity dominates.
  - **Bottom:** **Dirichlet** head :math:`h=0` m (water table at the base).
- **Initial condition:** dry (uniform suction) or user-specified profile stored in the example.
- **Structures / sinks:** off.

What to expect
--------------
- The wetting front rises from the bottom boundary; :math:`\theta(z,t)` increases first near
  the base and progresses upward until a steady profile is reached.
- Darcy flux :math:`q(z,t)` is upward near the advancing front and trends toward zero at steady state.
- Mass-balance diagnostics remain small throughout.

Outputs
-------
- **Profiles**: :math:`h(z,t)`, :math:`\theta(z,t)`, :math:`q(z,t)` under ``Figures/Profiles``.
- **Time series**: bottom boundary exchange, storage change, and mass-balance under ``Figures/TimeSeries``.
- **Data & logs**: numeric outputs in ``Data/`` and a detailed ``Log.txt`` with mesh, parameters, BCs, and tolerances.

Troubleshooting
---------------
- If the front stalls early, check that the bottom **Dirichlet** head is set to **0 m** and the top is **zero flux**.
- If Newton iterations climb, reduce the initial time step (``params.dt``) or raise ``params.max_iters``.
- Ensure your soil parameters are realistic (e.g., :math:`K_s`, :math:`\alpha`, :math:`n`)—extreme values can slow convergence.

See also
--------
- Back to :doc:`Examples index <index>`
- :doc:`../userguide/index`

.. _fig-ex03-caprise:

.. figure:: /static/examples/example03.png
   :alt: Example 03: Capillary rise profiles of h(z,t) and θ(z,t) from a fixed water table
   :align: center
   :width: 100%

   **Example 03** — Capillary rise from a fixed water table (bottom Dirichlet, top zero flux). Replace with
   your generated figure saved at ``_static/examples/Profiles/ex03_capillary_rise.png``.
