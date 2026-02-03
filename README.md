# Global climate risk indicators for outdoor sport (CMIP6)

This repository provides the full, reproducible code used to compute global climate risk indicators for outdoor sport from CMIP6 climate projections and to generate all figures presented in the associated manuscript.

The workflow is based on daily, bias-corrected and downscaled CMIP6 simulations from the NASA NEX-GDDP dataset and produces a coherent set of sport-relevant climate indicators, aggregated at seasonal scale and combined into a composite multi-hazard index.

---

## Scope and objectives

Outdoor sport practice is increasingly constrained by climate change through the combined effects of extreme heat, humidity, heavy rainfall and heatwaves.  
This repository supports a global, physically based assessment of these constraints by providing:

- transparent and reproducible indicator calculations,
- consistent treatment of multiple hazards,
- seasonally aggregated metrics relevant for sport calendars,
- a composite Sport Climate Risk Index (SCRI) synthesising multi-hazard exposure.

The code is designed to be reusable for future sport-, health- or adaptation-oriented studies.

---

## Climate indicators

The following indicators are computed:

- **Wet-Bulb Globe Temperature (WBGT)**  
  - morning, afternoon and evening conditions  
  - threshold: WBGT ≥ 32 °C (extreme risk)

- **Heat Index (HI)**  
  - NOAA formulation  
  - threshold: HI ≥ 40 °C (danger category)

- **Heavy rainfall days**  
  - threshold: precipitation ≥ 20 mm day⁻¹

- **Heatwave days**  
  - at least 3 consecutive days with daily maximum temperature ≥ 35 °C

Indicators are aggregated by month

and for three 30-year periods:
- historical (1991–2020),
- mid-century (2031–2060),
- late-century (2071–2100).

---

## Composite Sport Climate Risk Index (SCRI)

All indicators are normalised using global min–max scaling across all models, scenarios and grid cells, and combined into a composite index:

SCRI = 0.25 × WBGT + 0.25 × HI + 0.25 × rainfall + 0.25 × heatwaves

The SCRI ranges from 0 (low multi-hazard exposure) to 1 (high multi-hazard exposure) and is used to assess the simultaneous climatic constraints on outdoor sport practice.

---

## Data sources

- **Climate data**: NASA NEX-GDDP-CMIP6  
  (bias-corrected and statistically downscaled CMIP6 projections at 0.25° resolution)

- **CMIP6 models** (subset):
  - GFDL-ESM4  
  - IPSL-CM6A-LR  
  - MPI-ESM1-2-HR  
  - MRI-ESM2-0  
  - UKESM1-0-LL  

All results are presented as multi-model ensemble means.

---

## Repository structure

```text
├── code/
│   ├── indicators/        # WBGT, HI, rainfall, heatwave calculations
│   │── figures/           # scripts to reproduce all manuscript figures
│   └── utils/             # helper functions
│
│
├── README.md
├── LICENSE
