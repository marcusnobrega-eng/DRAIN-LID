# 🌧️ LID_Tool: A Physically Based 1D Richards Equation Model with Catchment Hydrology and Drainage Infrastructure

An integrated computational tool for modeling infiltration, ponding, vadose zone dynamics, and engineered drainage systems in urban and natural settings.

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)

---

## 📘 Project Summary

**LID_Tool** is a modular, mass-conservative, and high-resolution MATLAB model that solves the **mixed-form Richards equation** to simulate **vertical water movement** in variably saturated soils. It couples **surface flux forcing**, **drainage structures** (orifice + spillway), and **runoff generation** at the catchment scale.

### Built to support:
- 🌱 Ecohydrological modeling
- 🏙️ Green stormwater infrastructure design
- 🌊 Overflow and underdrain system simulation
- 📉 Soil moisture profile tracking
- 🧪 Experimental validation and benchmark reproduction

---

## 🔍 Model Capabilities

| Feature                   | Description                                                                 |
|--------------------------|-----------------------------------------------------------------------------|
| 🌊 Mixed-form Richards    | Solves ∂θ/∂t + Ss·∂h/∂t = ∂/∂z [K(h)(∂h/∂z + 1)]                            |
| 📦 Modular forcing engine | Supports rainfall, PET, runoff, ET reduction (Feddes), and hydrographs     |
| 💧 Ponding logic          | Automatically detects surface saturation, switches BC modes                |
| ⛲ Drainage systems       | Orifice and spillway with analytic flow formulas and activation thresholds |
| 🧠 Adaptive time stepping | Newton solver with line search (optional)                                  |
| 🧪 Pollutant dynamics     | Empirical buildup/washoff module driven by surface flux                    |
| 📊 Mass balance checks    | Cumulative conservation diagnostics at each time step                      |
| 📽️ Visualization suite   | Time series, vertical profiles, animations, FDCs                           |

---

## 📁 Directory Structure

```bash
📦 LID_Tool/
├── Main_Function/               # Main time loop, solver
├── Numerical_Solver/           # Residual & Jacobian computation
├── Physical_Functions_And_Utilities/
│   ├── Van Genuchten functions
│   ├── Mesh generation
│   └── Drainage structure logic
├── Forcing/                    # Rainfall, runoff (SCS-CN), PET, pollutants
├── Model_Configurations/       # Input scenarios and parameter scripts
├── Visualization_And_Output/  # Plots, exports, animations
├── Docs/                       # Manual, papers, diagrams
└── README.md
```markdown

---

## 📚 Documentation

A full user manual (LaTeX, PDF, or Overleaf) includes:

- 📖 Governing equations and boundary conditions
- 💻 Numerical methods (Newton, Jacobian, mass balance)
- 🧩 Mesh generation and time stepping
- ⚙️ Input parameters and scenario configuration
- 🧪 Validation benchmarks and outputs
- 🧠 Troubleshooting guidelines
- 📎 Appendices with VG parameters, units, and function index

---

## 📂 Output Directory Example

```bash
Results_<Scenario>/
├── Data/SimulationResults.xlsx
├── Figures/Plots/*.png
├── Figures/Animations/*.mp4
├── Log.txt


### Includes:
- ⬇️ Time series of $h$, $\theta$, $q$ at all nodes  
- 🧱 Profile plots at selected timesteps  
- 🌀 MP4 animations of vertical moisture/pressure flux  
- 📉 Flow duration curves  
- ✅ Mass balance error and cumulative diagnostics  

---

## 🧪 Built-In Case Studies

| Case                         | What it Shows                                  |
|-----------------------------|-------------------------------------------------|
| Celia1990.mat               | Validation vs. analytical benchmark             |
| Infiltration_LoamySand.mat  | Sharp wetting front in dry sandy soil           |
| PermeablePavement.mat       | Gravel + sand profile with orifice drainage     |
| Bioretention_Overflow.mat   | Surface ponding, overflow, ET, and underdrain   |
| CapillaryRise.mat           | Saturated bottom boundary and upward moisture   |

---

## 👤 Author

**Marcus N. Gomes Jr.**  
Postdoctoral Researcher II, University of Arizona  
📧 Email: [marcusnobrega.engcivil@gmail.com](mailto:marcusnobrega.engcivil@gmail.com)  
🌐 Website: [marcusnobrega-eng.github.io/profile](https://marcusnobrega-eng.github.io/profile)  
📄 CV: [Download PDF](https://marcusnobrega-eng.github.io/profile/files/CV___Marcus_N__Gomes_Jr_.pdf)  
🧪 ORCID: [0000-0002-8250-8195](https://orcid.org/0000-0002-8250-8195)  
🐙 GitHub: [@marcusnobrega-eng](https://github.com/marcusnobrega-eng)

---

## 📜 License

This project is licensed under the [MIT License](./LICENSE).

> You are free to use, modify, and distribute the model with attribution.
