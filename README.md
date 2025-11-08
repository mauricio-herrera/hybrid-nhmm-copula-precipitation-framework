# Hybrid Hidden Markov–Copula Framework for Chilean Rainfall (1980–2022)

**Repository purpose.** This repository provides data (compressed) and R scripts to reproduce the analyses and figures of the manuscript:

> **“Characterizing Climate-Driven Shifts in Chilean Rainfall Regimes with a Hybrid Hidden Markov–Copula Framework.”**

It satisfies journal reproducibility requirements by bundling: (i) a compressed dataset; (ii) self-contained R scripts; and (iii) clear instructions to configure file paths, run simulations for 2022 (MJJA), and export figures/tables.

---

## Abstract (from the manuscript)

Chile’s hydroclimate exhibits pronounced meridional gradients and strong interannual variability, posing persistent challenges for regime-aware, probabilistic rainfall prediction. We introduce a hierarchical framework that explicitly separates large-scale regime dynamics from local spatial dependence. The approach integrates: (i) a covariate-driven non-homogeneous Hidden Markov Model (nHMM) to learn synoptic precipitation regimes and their transitions; (ii) Dynamic Time Warping (DTW) clustering to delineate precipitation-coherent climatic zones; and (iii) state-conditional Regular Vine copulas with Generalized Pareto (GPD) tails to model residual spatial dependence and extremes. The analysis employs the 0.05° daily CR2MET precipitation product over continental Chile (462 grid points, May–August 1980–2021) together with large-scale atmospheric covariates including the Southern Oscillation Index (SOI), the Oceanic Niño Index (ONI), and Global Mean Sea-Surface Temperature (GMSST).

Five physically consistent rainfall regimes emerge, spanning from an anticyclonic dry state to a cyclonic wet state, confirmed by composites of mean sea-level pressure, 850-hPa winds, and 500-hPa geopotential height. Mixed-effects inference on the transition matrix reveals a statistically significant decline in wet-state persistence of ~0.34% yr⁻¹ (≈14.5% over 1980–2022), coincident with rising GMSST. Out-of-sample ensembles for 2022 (100 daily members conditioned on Viterbi states) are well calibrated: central 90% prediction intervals achieve near-nominal coverage, low asymmetry, and widths increasing southward with climatological variance.

By disentangling regime timing and drivers from residual spatial co-variability and extremes, the proposed nHMM–DTW–vine–GPD framework yields meteorologically coherent states, spatially consistent probabilistic simulations, and quantitatively validated forecasts. The method is computationally tractable and transferable, offering a principled pathway for regime-conditioned, uncertainty-aware precipitation prediction to support hydroclimate risk management in Chile and other topographically complex regions.

---

## Repository structure

Hybrid_Hidden_Markov–Copula_Framework/
├── Data/
│ └── mapa_estaciones_precipitaciones_clean.csv.zip (~34 MB, compressed; ~244 MB unzipped)
├── Scripts/
│ ├── Hybrid_nHMM_Copula_Simulation_Validation_2022.R
│ └── (other R helpers, if applicable)
├── (generated on first run)
│ ├── fig_main/
│ └── fig_supp/
├── README.md
└── .gitignore

- **Data/** contains a single compressed CSV (`.zip`). The uncompressed file is large (~244 MB) and **is not committed** to the repository.
- **Scripts/** contains the main R script to reproduce the 2022 MJJA validation runs, GoF diagnostics, multivariate scoring, and figure/table exports.

---

## Data dictionary (columns)

The dataset contains the following variables:

- `year`, `month`, `day`, `dates`: date fields (daily)
- `location`: auxiliary location ID
- `precipitation`: daily precipitation (mm)
- `codigo_estacion`: official station ID
- `z`: altitude (m)
- `x`, `y`: longitude, latitude
- `labels`: DTW-cluster/zone label
- `medoids`: binary flag (1 if station is the medoid of its cluster)
- `States`: hidden state for that day (0–4)
- `Year`: year (duplicated, convenience field)
- `year1980`: integer count since first date with records (baseline 1980)

---

## How to run (paths, decompression, outputs)

1. **Clone the repo and open R** (R ≥ 4.2 recommended).
2. **Set the project root** to the repository folder. We recommend using the [`here`](https://CRAN.R-project.org/package=here) package to avoid hard-coded absolute paths.
3. **Ensure the compressed data exists** at `Data/mapa_estaciones_precipitaciones_clean.csv.zip`.
4. **First run** will:
   - Unzip the CSV (≈244 MB) to `Data/mapa_estaciones_precipitaciones_clean.csv` (ignored by Git).
   - Create `fig_main/` and `fig_supp/` if not present.
   - Generate main figures and supplement tables/CSVs into those folders.

### Minimal path configuration (copy into the top of your R script)

```r
# ---------------- Path configuration ----------------
# Use 'here' to anchor paths at the project root:
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

data_dir <- here::here("Data")
zip_path <- file.path(data_dir, "mapa_estaciones_precipitaciones_clean.csv.zip")
csv_path <- file.path(data_dir, "mapa_estaciones_precipitaciones_clean.csv")

# Decompress if needed (CSV is ~244 MB; ZIP is ~34 MB)
if (!file.exists(csv_path)) {
  message("Unzipping dataset... (this may take a moment)")
  utils::unzip(zipfile = zip_path, exdir = data_dir, overwrite = TRUE)
  stopifnot(file.exists(csv_path))
}

# Read CSV (use readr for speed)
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
datos_base_l <- readr::read_csv(csv_path, progress = TRUE)
# ----------------------------------------------------
If you move the dataset elsewhere, simply change data_dir and zip_path / csv_path accordingly.
Do not hard-code absolute paths (e.g., /Users/<user>/...) in the committed script; keep them relative to the repo root.
Software requirements (R)

The scripts use the following packages (install on first use):

Data & pipes: dplyr, tidyr, purrr, tibble, data.table, stringr, lubridate, readr

Plotting: ggplot2, patchwork

Copulas & extremes: VineCopula, copula, evd

Scoring & tests: scoringRules, forecast, goftest, ADGofTest

Parallel: parallel

Utilities: pacman (optional), here (paths)

You can bootstrap with:
pkgs <- c(
  "dplyr","tidyr","purrr","tibble","data.table","stringr","lubridate","readr",
  "ggplot2","patchwork",
  "VineCopula","copula","evd",
  "scoringRules","forecast","goftest","ADGofTest",
  "parallel","here","pacman"
)
to_install <- pkgs[!pkgs %in% installed.packages()[,1]]
if (length(to_install)) install.packages(to_install)
Reproducible workflow

Load data (unzips on first run).

Train margins by state and station (zero-inflated + ECDF body + GPD tail).

Build U-matrices per state (train) and fit R-vine copulas (state-conditional).

Diagnostics: PIT (univariate & Rosenblatt), White tests on Tree-1 pairs.

Simulation (2022 MJJA): state-conditioned daily ensembles (vine vs independence).

Verification: ES/VS, areal-CRPS, Brier events (≥k wet, joint exceedances), DM tests.

Exports: figures to fig_main/ (paper) and fig_supp/ (supplement), CSV tables to fig_supp/.

Outputs

Main figures: fig_main/ (publication panels).

Supplement figures & CSVs: fig_supp/ (PIT panels, ES/VS/CRPSA gains, Brier gains, GoF summaries).

The exact filenames are created by the script; they follow the manuscript’s S-section naming.


