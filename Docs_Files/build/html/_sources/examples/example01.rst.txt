.. _example01:

Example 01 — Infiltration into Sandy Soil (Dirichlet head)
==========================================================

Goal
----
Demonstrate a sharp wetting front produced by constant **pressure head** infiltration
into a dry **sand** column and compare against HYDRUS-1D.

How to run
----------
Uncomment the line below in ``Mixed_Richards_Model_1D.m`` and run the script::


       % load('Examples/Celia1990.mat');
       load('Examples/Example1_Infiltration_Sand.mat');
       % load('Examples/Example2_Clay_Loam_Soil.mat');
       % ...

Scenario at-a-glance
--------------------
- Column depth: :math:`L = 1.0` m; nodes: :math:`N_z = 41` (uniform). :contentReference[oaicite:1]{index=1}
- Time setup: total time :math:`T_{\max}=3600` s (1 h); initial :math:`\Delta t=1` s; save every 60 s.
- Soil (USDA sand): :math:`\theta_r=0.045`, :math:`\theta_s=0.430`, :math:`\alpha=14.5\,\text{m}^{-1}`, :math:`n=2.68`, :math:`K_s=8.25\times10^{-5}\,\text{m s}^{-1}`, :math:`S_s=10^{-5}\,\text{m}^{-1}`.
- Initial condition: hydrostatic (dry near surface).
- Boundary conditions:
  - **Top**: Dirichlet (fixed head) :math:`h_\text{top}=0` m.
  - **Bottom**: free drainage, :math:`\partial h/\partial z = 0`.
- No drainage structures or internal sinks are active.

What to expect
--------------
- A sharp wetting front initiates at the surface and propagates downward rapidly;
  no surface ponding is expected for this sandy profile. :
- Agreement with HYDRUS-1D is strong for both :math:`h(z,t)` and :math:`\theta(z,t)` (low RMSE),
  with small deviations mainly at very early/late times near the surface.
- The manual includes comparison plots (RMSE and profiles) for this case.

Output
------
Running this example will save time series and vertical profiles (head, water content, flux)
and summary figures into the scenario’s output folders (see *Output and Diagnostics* for
file layout).

.. _fig-ex01-wetting-front:

.. figure:: /static/examples/example01.png
   :alt: Wetting-front advance in sandy soil for Example 01
   :align: center
   :width: 100%

   Wetting-front advance in **Example 01** (sandy soil). Profiles of :math:`h(z,t)` and :math:`\theta(z,t)`
   show a sharp, downward-moving front. See :ref:`example01` for setup.


