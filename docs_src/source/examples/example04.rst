.. _example04:

Example 04 — Constant Rainfall (Neumann 1 m/day) in Sandy Soil
==============================================================

Goal
----
Evaluate the model under a **constant downward surface flux** of **1 m/day**
into a dry **sandy** profile. This is a high-conductivity case with **no surface
ponding**. Results are compared against HYDRUS-1D.

How to run
----------
Uncomment **exactly one** line in ``Mixed_Richards_Model_1D.m``:

.. code-block:: matlab

   % load('Examples/Example1_Infiltration_Sand.mat');
   % load('Examples/Example2_Clay_Loam_Soil.mat');
   % load('Examples/Example3_CapillaryRise.mat');
   load('Examples/Example4_TopNeumann_Sandy.mat');   % ← enable this

Press **Run** (F5). Outputs are written to ``Modeling_Results/<sim_name>/`` (Figures, Data, Log.txt).

Setup (summary)
---------------
- **Same inputs as Example 01** (sandy soil, column depth, mesh), **except** the top BC.
- **Top boundary:** **Neumann** constant flux
  :math:`q_{\text{top}} = 1\,\mathrm{m\,day^{-1}} \approx 1.1574\times10^{-5}\,\mathrm{m\,s^{-1}}`
  (negative = downward).
- **Bottom boundary:** **free drainage** (:math:`\partial h/\partial z=0`).
- **Evaporation/ET:** disabled.
- **Initial condition:** dry suction (as in Example 01).

What to expect
--------------
- A rapid, sharp **wetting front** propagating downward; **no ponding** at the surface.
- **Good agreement** with HYDRUS-1D in both pressure head :math:`h(z,t)` and moisture
  :math:`\theta(z,t)`; **low RMSE** for most of the simulation.
- Slight differences may appear near the surface and the advancing front at later times
  (boundary handling, mesh discretization).

Outputs
-------
- **Profiles**: :math:`h(z,t)`, :math:`\theta(z,t)`, :math:`q(z,t)` saved under ``Figures/Profiles``.
- **Time series**: surface flux (rain/zero), bottom seepage, storage change, and mass-balance in
  ``Figures/TimeSeries``; CSV/MAT in ``Data``; run details in ``Log.txt``.

See also
--------
- Back to :doc:`Examples index <index>`
- :doc:`example01` (same soil/geometry with Dirichlet top)
- :doc:`../userguide/index`

.. _fig-ex04-compare:

.. figure:: /static/examples/example04.png
   :alt: Example 04: Comparison of DRAIN-LID and HYDRUS-1D under constant 1 m/day rainfall
   :align: center
   :width: 100%

   **Example 04** — RMSE of head (top-left) and moisture (bottom-left) over time; right panels show
   profile comparisons at 10 %, 25 %, 50 %, 75 %, and 100 % of the run. Solid = HYDRUS-1D, dotted = DRAIN-LID.
