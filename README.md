# 🌧️ DRAIN-GI: Darcy–Richards Analysis of Infiltration in Nature-based Low Impact Development

A MATLAB-based hydrological modeling framework for simulating infiltration, evaporation, and drainage in green infrastructure systems such as permeable pavements, green roofs, and bioretention cells.

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)

---

## 📘 Project Summary

**DRAIN-GI** is a modular, mass-conservative, and high-resolution hydrological model that solves the **1D mixed-form Richards equation** to simulate **vertical water movement** in variably saturated soils. It is designed for system-level analysis of **nature-based Low Impact Development (LID)** systems under real-world climatic and hydrological conditions.

The model couples **infiltration modeling**, **surface fluxes**, **evaporation and evapotranspiration**, and **engineered drainage infrastructure** into a unified and extensible MATLAB-based platform.

### Designed for:
- 🌱 Green infrastructure and sustainable stormwater planning through design events or continuous simulations
- 🧪 Richards-based unstatured zone modeling adapted to infiltration and saturation excess overland flow generation mechanisms
- 🏙️ Urban LID system design and analysis (permeable, green roofs, bioretention systems, or any other infiltration-based LID system)
- 📈 Climate-resilient performance evaluation of GI over multi-decadal timescales through long-term continuous simulations to obtain hydrologic signatures

---

## 🔍 Model Capabilities

| Feature                   | Description                                                                 |
|--------------------------|------------------------------------------------------------------------------|
| 🌊 Mixed-form Richards    | Solves ∂θ/∂t + Ss·∂h/∂t = ∂/∂z [K(h)(∂h/∂z + 1)]                            |
| 📦 Modular forcing engine | Supports rainfall, PET (Hargreaves or user-specified input ETP), runoff (SCS-CN for LID systems that receive upstream catchment runoff), and user-defined hydrograph or rainfall forcing at LID surface |
| 💧 Ponding logic          | Automatically detects saturation and switches between Neumann/Dirichlet BCs |
| ⛲ Drainage structures    | Orifice and spillway modules with activation                                |
| 🧠 Adaptive time stepping | Newton–Raphson solver with mass balance convergence and line search         |
| 🌫️ Evaporation modeling   | Feddes function for ET reduction based on pressure head                     |
| 🧪 Pollutant module       | Empirical buildup and washoff routines driven by surface runoff             |
| 🧮 Catchment runoff       | Lumped kinematic flow model with SCS-CN and Manning’s-based outflow         |
| 📊 Mass balance checks    | Cumulative conservation diagnostics at every time step                      |
| 📽️ Visualization suite    | Animations, profiles, time series, flow duration curves (FDCs)              |

---

## 📁 Directory Structure

bash
📦 DRAIN-GI/
├── Main_Function/               # Main solver and time integration loop
├── Numerical_Solver/           # Residual, Jacobian, and Newton solver functions
├── Physical_Functions_And_Utilities/
│   ├── Van Genuchten functions
│   ├── Mesh generator
│   └── Drainage and ponding logic
├── Forcing/                    # Rainfall, runoff (SCS-CN), PET, pollutants
├── Model_Configurations/       # Input parameter sets and predefined test cases
├── Visualization_And_Output/  # Plots, animations, diagnostics
├── Docs/                       # Technical Manual and Scientific Paper Draft (currently under review)
└── README.md

## 📚 Documentation
A full user manual is included (Docs/Manual_DRAIN_GI.pdf), covering:

📖 Mathematical formulation and boundary condition types

💻 Numerical methods (discretization, Jacobian, Newton solver)

🧩 Mesh design and adaptive time stepping

🧪 Built-in benchmarks vs. Hydrus-1D and experimental data

⚙️ Parameter configuration for LID systems (green roofs, permeable pavements, bioretention)

📈 Output interpretation and performance metrics

💡 Common issues and troubleshooting

## Results_<Scenario>/
bash
├── Data/SimulationResults.xlsx (Excel spreadsheet with all model states and outputs for all nodes in the domain)
├── Figures/Plots/*.png 
├── Figures/Animations/*.mp4
├── Log.txt

## Key outputs:
⬇️ Time series of pressure head ($h$), water content ($\theta$), and flux ($q$)

📊 Profile plots and snapshots at selected time steps

📽️ MP4 animations of infiltration and vertical flux

📈 Flow duration and evaporation duration curves

✅ Cumulative water balance diagnostics and error checks

## 🧪 Case Studies and Applications

| Example | Scenario ID                    | Description                                       | Soil Type   | Duration | Boundary Conditions (Top / Bottom) | Validation       |
|---------|--------------------------------|--------------------------------------------------|-------------|----------|-------------------------------------|------------------|
| 1️⃣      | `Example1_Infiltration_Sand`   | Sharp wetting front from constant infiltration   | Sand        | 1 hr     | Dirichlet / Free drainage           | HYDRUS-1D ✅    |
| 2️⃣      | `Example2_Clay_Loamy_Soil`     | Diffuse front in low-K soil                      | Clay Loam   | 10 hrs   | Dirichlet / Free drainage           | HYDRUS-1D ✅    |
| 3️⃣      | `Example3_CapillaryRise`       | Upward capillary rise from saturated bottom      | Sandy Loam  | 10 hrs   | Dirichlet / Dirichlet               | HYDRUS-1D ✅    |
| 4️⃣      | `Example4_Rainfall_Sand`       | Infiltration under constant rainfall flux        | Sand        | 1 hr     | Neumann (1 m/day) / Free drainage   | HYDRUS-1D ✅    |

Several other examples are available in the \Examples folder. To run a particular example, please just load the scenario in the main code function.


👤 Developer
**Marcus N. Gomes Jr.**  
Postdoctoral Researcher II, University of Arizona  
📧 Email: [marcusnobrega.engcivil@gmail.com](mailto:marcusnobrega.engcivil@gmail.com)  
🌐 Website: [marcusnobrega-eng.github.io/profile](https://marcusnobrega-eng.github.io/profile)  
📄 CV: [Download PDF](https://marcusnobrega-eng.github.io/profile/files/CV___Marcus_N__Gomes_Jr_.pdf)  
🧪 ORCID: [0000-0002-8250-8195](https://orcid.org/0000-0002-8250-8195)  
🐙 GitHub: [@marcusnobrega-eng](https://github.com/marcusnobrega-eng)
