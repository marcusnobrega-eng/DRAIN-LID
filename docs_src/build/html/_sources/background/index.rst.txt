.. _background:

==========
Background
==========

Motivation
----------
**DRAIN-LID (Darcy–Richards Infiltration model for Low-Impact Development)** simulates vertical water movement in layered LID systems (e.g., permeable pavements, green roofs, bioretention media) under rainfall and evaporation/ET forcing. The goal is to resolve infiltration, storage, and exfiltration dynamics with a physics-based 1-D solver that is transparent, reproducible, and easy to benchmark (e.g., against HYDRUS-1D) while remaining lightweight enough for long-term scenario runs.

Name
----
**DRAIN-LID** expands to *Darcy–Richards Infiltration model for Low-Impact Development*. It emphasizes that the core engine is the Darcy–Richards equation, specialized for vertical infiltration in LID layers.

Scope and assumptions
---------------------
- 1-D vertical column (no lateral flow within layers).
- Layered, homogeneous properties per layer; hysteresis not modeled (optional in roadmap).
- Isothermal; water phase only (no air entrapment or heat transport).
- Runoff is generated when surface flux exceeds infiltration + storage capacity.
- Suitable for permeable pavements, green roofs, and bioretention media where vertical processes dominate.


Theory (governing equation)
---------------------------
We solve the **mixed-form Richards equation** in 1-D vertical coordinates (:math:`z` positive upward):

.. math::
   :label: eq:richards_mixed

   C(h)\,\frac{\partial h}{\partial t}
   \;=\;
   \frac{\partial}{\partial z}\!\left[
      K(h)\,\Big(\frac{\partial h}{\partial z} + 1\Big)
   \right]
   \;-\; S(h,t),

where :math:`h(z,t)` is pressure head [m], :math:`C(h)=\mathrm{d}\theta/\mathrm{d}h` is the specific moisture capacity, :math:`\theta(h)` is volumetric water content, :math:`K(h)` is unsaturated hydraulic conductivity [m/s], and :math:`S(h,t)` is a sink term (e.g., root uptake; often 0 for hardscape LID).

Constitutive relations (van Genuchten–Mualem, with units)
--------------------------------------------------------

We use van Genuchten–Mualem for water retention and conductivity (1-D, pressure head convention :math:`h<0` in the unsaturated zone so :math:`|h|` appears below):

.. math::

   S_e(h) \;=\; \bigl[\,1+(\alpha\,|h|)^n\,\bigr]^{-m},
   \qquad m \equiv 1 - \frac{1}{n} \quad [-]

.. math::

   \theta(h) \;=\; \theta_r \;+\; S_e(h)\,\bigl(\theta_s-\theta_r\bigr)
   \qquad [\theta] = [-]

.. math::

   K(h) \;=\; K_s \, S_e(h)^{\ell}\,
   \bigl[\,1-\bigl(1-S_e(h)^{1/m}\bigr)^m\,\bigr]^2
   \qquad [K] = [K_s] = \mathrm{m\,s^{-1}}

Parameters (per layer):
:math:`\theta_r, \theta_s` (residual/saturated water content, :math:`[-]`),
:math:`\alpha` (inverse air-entry suction, :math:`\mathrm{m^{-1}}`),
:math:`n>1` (shape, :math:`[-]`), :math:`m=1-1/n` (shape, :math:`[-]`),
:math:`\ell` (pore-connectivity, default :math:`\ell=0.5`, user-configurable, :math:`[-]`),
:math:`K_s` (saturated hydraulic conductivity, :math:`\mathrm{m\,s^{-1}}`).

Storage term with specific storage
----------------------------------
The moisture capacity from retention is :math:`C_\theta(h)=\mathrm{d}\theta/\mathrm{d}h \; [\mathrm{m^{-1}}]`.
With compressibility, we include **specific storage** :math:`S_s \; [\mathrm{m^{-1}}]` and use

.. math::

   C_{\text{tot}}(h) \;=\; C_\theta(h) \;+\; S_s \qquad [\mathrm{m^{-1}}].

In the mixed-form Richards equation we therefore use :math:`C_{\text{tot}}(h)\,\partial h/\partial t` as the storage term.


Boundary and initial conditions
-------------------------------
**Top boundary (rain/ET with ponding switch).** A flux (Neumann) is imposed from rainfall/evaporation:

.. math::

   q_\mathrm{surf}(t) \;=\; q_\mathrm{rain}(t) - E(t),

converted to an infiltration flux :math:`q_0(t)` until infiltration capacity is exceeded. When :math:`q_\mathrm{rain}` would cause **ponding**, the model switches to a **Dirichlet** condition with :math:`h_0 \ge 0` (ponded head), and **runoff** is generated as the residual:

.. math::

   q_\mathrm{runoff}(t) \;=\; \max\!\big(0,\; q_\mathrm{rain}(t) - q_\mathrm{infil}(t)\big).

**Bottom boundary.** Either free drainage :math:`\partial h/\partial z = 0`, a fixed head (e.g., water table), or a specified flux.

**Initial condition.** :math:`h(z,0)` (or :math:`\theta(z,0)`) specified per layer; hydrostatic or user-defined.

**Drainage structures (optional)**

- *Orifice flow* (submerged opening): :math:`Q_o = C_d\,A_o\,\sqrt{2g\,H}`  [m³ s⁻¹],
  where :math:`C_d` [-] is the discharge coefficient, :math:`A_o` [m²] the area, :math:`H` [m] the head.

- *Broad-crested weir / spillway*: :math:`Q_s = C_w\,b\,H^{3/2}`  [m³ s⁻¹],
  with :math:`C_w` [-] coefficient, :math:`b` [m] width, :math:`H` [m] head.

Structure outflows are volume-consistent with the surface store and reduce ponded depth/head accordingly.


Evapotranspiration / Evaporation modeling
----------------------------------------

**Reference ET (Hargreaves).**
The Hargreaves method estimates reference evapotranspiration :math:`ET_0` from temperature and
extraterrestrial radiation :math:`R_a`:

.. math::

   ET_0 \;=\; 0.0023\,\big(T_{\text{mean}} + 17.8\big)\,
   \big(T_{\text{max}} - T_{\text{min}}\big)^{1/2}\, R_a
   \qquad [\mathrm{mm\,day^{-1}}],

where :math:`T_{\text{mean}}`, :math:`T_{\text{max}}`, :math:`T_{\text{min}}` are daily air temperatures [°C],
and :math:`R_a` is extraterrestrial radiation [MJ m^{-2} day^{-1}] (function of latitude and day of year).

Functions to estimate hargreaves reference evapotranspiration are available in the folder Physical_Functions_and_Utilities, although users are required to enter their time-series in the forcing file.

**Moisture limitation (Feddes reduction).**
To account for water availability limits near the surface, the potential vertical flux is reduced by the
Feddes stress function :math:`\alpha(h)\in[0,1]` based on surface pressure head :math:`h` [m]:

.. math::

   q_{\mathrm{ET}} \;=\; \alpha\!\big(h\big)\;ET_0
   \qquad [\mathrm{mm\,day^{-1}}].

A two-threshold, piecewise-linear form is used with user-specified wet and dry limits
:math:`h_{\text{lim,upper}}` and :math:`h_{\text{lim,lower}}` (both in meters of water head):

.. math::

   \alpha(h) \;=\;
   \begin{cases}
     1, & h \ge h_{\text{lim,upper}},\\[6pt]
     \dfrac{\,h - h_{\text{lim,lower}}\,}{\,h_{\text{lim,upper}} - h_{\text{lim,lower}}\,},
       & h_{\text{lim,lower}} < h < h_{\text{lim,upper}},\\[10pt]
     0, & h \le h_{\text{lim,lower}}.
   \end{cases}

**Permeable pavements (empirical potential evaporation).**
For non-vegetated systems such as permeable pavements, observations suggest evaporation is
approximately 30 % of gross pan evaporation. We therefore define an empirical potential evaporation
:math:`E_0` as

.. math::

   E_0 \;=\; \gamma\,\big(ET_0 + \delta\big)
   \qquad [\mathrm{mm\,day^{-1}}],

with :math:`\gamma=0.3` (evaporation ratio) and daily bias :math:`\delta=2` mm day⁻¹.
Under pavement conditions, the effective evaporation flux replaces :math:`ET_0` with :math:`E_0` in the Feddes
relation: :math:`q_{\mathrm{ET}} = \alpha(h)\,E_0`.

**Unit note for the solver.**
Fluxes supplied to the Richards solver use SI units [m s⁻¹]. Convert daily depths via
:math:`q\,[\mathrm{m\,s^{-1}}] = (q\,[\mathrm{mm\,day^{-1}}] \times 10^{-3})/86400`.


Numerical method
----------------
DRAIN-LID uses a fully implicit discretization of :eq:`eq:richards_mixed` with Newton–Raphson iterations (line search) and adaptive time stepping. The Jacobian is assembled from :math:`C(h)` and :math:`K(h)` sensitivities. Tolerances control nonlinear convergence; a **running mass balance** diagnostic is reported at each step.

Catchment hydrology and surface coupling
----------------------------------------

We combine (i) **SCS–CN abstraction** for effective rainfall, (ii) a **lumped, non-linear
reservoir** for runoff generation/routing, and (iii) a **top Neumann coupling** to DRAIN-LID.

1) SCS–CN abstraction (event/continuous)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For an event depth :math:`P_e` [m]:

.. math::

   S \;=\; \left(\frac{25400}{\mathrm{CN}} - 254\right)\,\mathrm{mm}, \qquad
   I_a \;=\; \lambda\,S, \quad 0.05 \le \lambda \le 0.2,

.. math::

   R_e \;=\;
   \begin{cases}
   \dfrac{(P_e - I_a)^2}{P_e - I_a + S}, & P_e > I_a,\\[6pt]
   0, & \text{otherwise,}
   \end{cases}

yielding the event runoff depth :math:`R_e` [m] (convert :math:`S` mm → m when used with SI).
For **continuous series**, apply the same formula to the **cumulative precipitation**
:math:`P_{\text{cum}}(t)` to obtain cumulative runoff :math:`R_{\text{cum}}(t)`, then compute
the **effective rainfall rate** as the time derivative:

.. math::

   P_{ef}(t) \;=\; \frac{d R_{\text{cum}}}{dt}
   \;\;\approx\;\; \frac{R_{\text{cum}}(t_{k}) - R_{\text{cum}}(t_{k-1})}{\Delta t}
   \qquad [\mathrm{m\,s^{-1}}].

2) Lumped reservoir (kinematic outflow)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Treat the watershed as a non-linear reservoir with mean water depth :math:`h^w(t)` [m] and a
weir-like threshold :math:`h_0^w` [m] for initial abstractions. With kinematic-wave slope
:math:`s_0` (≈ bed slope) and Manning roughness :math:`n^w`, the discharge is

.. math::

   q^w(t) \;=\; \frac{1}{n^w}\, w^w\,
   \max\!\big(h^w(t)-h_0^w,\,0\big)^{5/3}\, s_0^{1/2}
   \qquad [\mathrm{m^3\,s^{-1}}],

where :math:`w^w` [m] is an effective width. The mass balance per unit area :math:`a^w` [m²] is

.. math::

   \frac{d h^w}{dt} \;=\; P_{ef}(t) \;-\; \frac{q^w(t)}{a^w}
   \qquad [\mathrm{m\,s^{-1}}].

A forward-Euler step advances :math:`h^w`:

.. math::

   h^w_{k+1} \;=\; h^w_k \;+\; \Delta t
   \Big[P_{ef,k} \;-\; \frac{q^w_k}{a^w}\Big].

3) Coupling to DRAIN-LID (top Neumann BC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The watershed discharge is converted to a depth flux over the receiving LID area
:math:`A_{\text{LID}}` [m²] and combined with local rain–ET to define the **net surface flux**:

.. math::

   q_{\text{top}}(t) \;=\; \big(\mathrm{Rain}(t) - ET_0(t)\big)
   \;-\; \frac{q^w(t)}{A_{\text{LID}}}
   \qquad [\mathrm{m\,s^{-1}}],

with **positive** :math:`q_{\text{top}}` defined **downward** (into the column).

**Manual inflow option.**
When a manual mode is selected (e.g., ``flag_manual_inflow = 1``), bypass the watershed term:
:math:`q_{\text{top}}(t) = \mathrm{Rain}(t) - ET_0(t)` (converted to :math:`\mathrm{m\,s^{-1}}`).

**Units note.**
Spreadsheet inputs are typically **mm/h**; convert with
:math:`q[\mathrm{m\,s^{-1}}] = (\text{mm/h}\times 10^{-3})/3600`.


Water quality (optional): buildup–washoff
----------------------------------------

An empirical pollutant module simulates surface mass dynamics using **buildup** on dry days and **washoff** during runoff.
Let :math:`B(t)` be pollutant mass per unit area [M L⁻²] (e.g., g m⁻²).

**Buildup (dry weather)**

.. math::

   \frac{dB}{dt} \;=\; k_b\,\big(B_{\max} - B\big)
   \qquad
   \big[k_b\big] = \mathrm{s^{-1}},\;\; \big[B_{\max}\big] = \mathrm{M\,L^{-2}}.

**Washoff (during runoff)**

.. math::

   W(t) \;=\; k_w\,Q_{\mathrm{surf}}(t)^{\,p}\,B(t)
   \qquad
   \big[k_w\big] = \mathrm{L^{p-1}\,T^{-1}},\;\; p\ge 1

where :math:`Q_{\mathrm{surf}}(t)` is surface discharge per unit width [L² T⁻¹] or depth-based proxy; :math:`W(t)` is areal washoff rate [M L⁻² T⁻¹].
The surface concentration discharged is :math:`C_{\mathrm{out}} = W/q_{\mathrm{runoff}}` [M L⁻³] when runoff exists. Optional first-order decay or
media filtration can be applied within layers during infiltration.

The current version of this model only generate inflow pollutographs that later must be routed in the porous media through an advection-dispersion or sediment-based pollutant model (under development, not available in the current model).

Numerical method: residual, Jacobian, adaptive time stepping
------------------------------------------------------------

Let :math:`z` [m] be positive upward, :math:`\Delta z_i` [m] the cell height, and :math:`\Delta t^n` [s] the time step.
With **specific storage** :math:`S_s` and moisture capacity :math:`C_\theta(h) = \frac{d\theta}{dh}`,
the total storage is :math:`C_{\text{tot}}(h) = C_\theta(h) + S_s`.
We use a **fully implicit** scheme:

**Interior node residual** (node :math:`i`, time :math:`n+1`):

.. math::

   F_i = C_{\text{tot},i}^{\,n+1} \frac{h_i^{n+1} - h_i^{n}}{\Delta t^n}
   - \frac{1}{\Delta z_i} \left[ K_{i+\frac12}^{\,n+1}
      \left( \frac{h_{i+1}^{n+1} - h_i^{n+1}}{\Delta z_{i+\frac12}} + 1 \right)
      - K_{i-\frac12}^{\,n+1}
      \left( \frac{h_i^{n+1} - h_{i-1}^{n+1}}{\Delta z_{i-\frac12}} + 1 \right) \right]
   + S_i^{\,n+1},

where :math:`K_{i\pm\frac12}` are face conductivities and :math:`S_i` is a sink term.
Boundary rows enforce the chosen BCs (flux or head), including the **ponding switch** at the surface.

**Newton–Raphson solver**

At each time step, solve :math:`F(h^{n+1}) = 0` via:

.. math::

   J(h^k) \, \delta h^k = -F(h^k),
   \qquad h^{k+1} = h^k + \lambda \, \delta h^k,

with line search :math:`0 < \lambda \le 1`.
The **Jacobian** :math:`J = \frac{\partial F}{\partial h}` is **sparse, banded (tridiagonal)**;
entries contain derivatives of :math:`C_{\text{tot}}(h)` and flux terms w.r.t. neighboring heads.
A symbolic 5-node system illustrates the nonzeros: :math:`\{ J_{i,i-1}, J_{i,i}, J_{i,i+1} \}` for interior rows;
boundary rows adjust per BC.

**Adaptive time stepping & convergence**
Time step :math:`\Delta t` grows when Newton converges quickly and shrinks if iterations exceed thresholds or if mass-balance error is high.
Typical tolerances: residual norm and head update norm below user-set values.

**Mass balance (diagnostic)**
For each step we compute cumulative **inflows** (:math:`P`, upstream/structure inflow), **outflows** (:math:`q_{\mathrm{runoff}}`, bottom drainage, ET),
and **storage change** :math:`\Delta S = \sum_i (\theta_i^{n+1}-\theta_i^{n})\,\Delta z_i` [m]. The relative mass-balance error is reported as

.. math::

   \varepsilon_{\mathrm{MB}} \;=\;
   \frac{\big|\;(\mathrm{In} - \mathrm{Out}) - \Delta S\;\big|}{\max(10^{-6}\,\mathrm{m},\,\mathrm{In})} \;\;\;[-].

References
----------
.. [1] Richards, L. A. (1931). Capillary conduction of liquids through porous mediums. *Physics*, 1, 318–333.
.. [2] van Genuchten, M. T. (1980). A closed‐form equation for predicting the hydraulic conductivity of unsaturated soils. *Soil Sci. Soc. Am. J.*, 44, 892–898.
.. [3] Mualem, Y. (1976). A new model for predicting the hydraulic conductivity of unsaturated porous media. *Water Resour. Res.*, 12(3), 513–522.
.. [4] Celia, M. A., Bouloutas, E. T., & Zarba, R. L. (1990). A general mass-conservative numerical solution for the unsaturated flow equation. *Water Resour. Res.*, 26(7), 1483–1496.
.. [5] Šimůnek, J., van Genuchten, M. T., & Šejna, M. (2016). HYDRUS: Model use, calibration, and validation. *Trans. ASABE*, 55, 1261–1274.
.. [6] Hargreaves, G. H., & Samani, Z. A. (1985). Reference crop evapotranspiration from temperature. *Applied Engineering in Agriculture*, 1(2), 96–99.
.. [7] Farthing, M. W., & Ogden, F. L. (2017). Numerical solution of Richards’ equation: A review. *Vadose Zone Journal*, 16(7).
