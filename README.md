# 🌧️ LID_Tool: A Physically Based 1D Richards Equation Model with Catchment Hydrology and Drainage Infrastructure

An integrated computational tool for modeling infiltration, ponding, vadose zone dynamics, and engineered drainage systems in urban and natural settings.

---

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)

---

## 📘 Project Summary

LID_tool is a modular, mass-conservative, and high-resolution MATLAB model that solves the **mixed-form Richards equation** to simulate **vertical water movement** in variably saturated soils. It couples **surface flux forcing**, **drainage structures** (orifice + spillway), and **runoff generation** at the catchment scale.

It is built to support:
- 🌱 Ecohydrological modeling
- 🏙️ Green stormwater infrastructure design
- 🌊 Overflow and underdrain system simulation
- 📉 Soil moisture profile tracking
- 🧪 Experimental validation and benchmark reproduction

---

## 🔍 Model Capabilities

| Feature                     | Description                                                                 |
|----------------------------|-----------------------------------------------------------------------------|
| 🌊 Mixed-form Richards      | Solves ∂θ/∂t + Ss·∂h/∂t = ∂/∂z [K(h)(∂h/∂z + 1)]                            |
| 📦 Modular forcing engine   | Supports rainfall, PET, runoff, ET reduction (Feddes), and hydrographs      |
| 💧 Ponding + Neumann logic  | Automatically detects surface saturation, switches between BC modes         |
| ⛲ Drainage systems         | Supports orifice and spillway drainage with analytic formulas and caps      |
| 🧠 Adaptive time stepping   | Efficient Newton solver with Jacobian updates, line search (optional)       |
| 🧪 Pollutant dynamics       | Empirical buildup/washoff model driven by surface water flux                |
| 📊 Mass balance tracking    | Strict conservation check at every step with cumulative diagnostics         |
| 📽️ Visualization suite     | Profile plots, time series, animations, FDCs — publication-quality outputs  |

---


## 📁 Directory Structure

📦 HydroVadose/
├── Main_Function/              # Main time loop, solver, data orchestration
├── Numerical_Solver/          # Residual & Jacobian computation
├── Physical_Functions_And_Utilities/
│   ├── Van Genuchten functions
│   ├── Mesh generation
│   └── Drainage structure logic
├── Forcing/                   # Rainfall, runoff (SCS-CN), hydrographs, ETP
├── Model_Configurations/      # Input scenarios and parameter scripts
├── Visualization_And_Output/ # Plots, exports, animations
├── Docs/                      # Manual, papers, diagrams
└── README.md

📚 Documentation
Full user manual (LaTeX/PDF/Overleaf) includes:

📖 Governing equations and physics

💻 Numerical methods (Newton, Jacobian)

🧩 Mesh setup and time stepping

⚙️ Input parameters and configuration guide

🧪 Validation benchmarks

🧠 Troubleshooting section

📎 Appendices: VG parameter tables, function index, dimensional consistency

---

Results_<scenario>/
├── Data/SimulationResults.xlsx
├── Figures/Plots/*.png
├── Figures/Animations/*.mp4
├── Log.txt

⬇️ Time series of $h$, $\theta$, $q$ at all nodes

🧱 Profile plots at selected timesteps

🌀 MP4 animations of vertical moisture/pressure flux

📉 Flow duration curves

✅ Mass balance error and cumulative diagnostics

---

🧪 Built-In Case Studies

| Case                        | What it Shows                             |
| --------------------------- | ----------------------------------------- |
| Celia1990.mat               | Validation vs. literature benchmark       |
| Infiltration\_LoamySand.mat | Wetting front behavior in dry soil        |
| PermeablePavement.mat       | Gravel + sand + orifice flow              |
| Bioretention\_Overflow\.mat | Surface ponding, overflow + ET + drainage |
| CapillaryRise.mat           | Upward flow from saturated boundary       |

---

👋 Author
Marcus N. Gomes Jr.
Postdoctoral Researcher II, University of Arizona
📧 marcusnobrega.engcivil@gmail.com

MIT License

Copyright (c) 2025 Marcus N. Gomes Jr.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell   
copies of the Software, and to permit persons to whom the Software is         
furnished to do so, subject to the following conditions:                      

The above copyright notice and this permission notice shall be included in    
all copies or substantial portions of the Software.                           

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN     
THE SOFTWARE.

