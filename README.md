# Hybrid Hidden Markov–Copula Framework for Chilean Rainfall (1980–2022)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23322615.svg)](https://doi.org/10.5281/zenodo.23322615)

Research software and processed precipitation input for the manuscript:

> **Characterizing Climate-Driven Shifts in Chilean Rainfall Regimes with a Hybrid Hidden Markov–Copula Framework**

This repository provides an R implementation for state-conditional precipitation simulation and out-of-sample validation during May–August 2022. The scientific framework couples non-homogeneous Hidden Markov Model (nHMM) precipitation regimes, Dynamic Time Warping (DTW) climatic zones, zero-inflated empirical/GPD station margins, and regular vine copulas.

**Archived release:** GitHub [v1.0.0](https://github.com/mauricio-herrera/hybrid-nhmm-copula-precipitation-framework/releases/tag/v1.0.0), published on 12 October 2026 and archived by Zenodo at **[DOI: 10.5281/zenodo.23322615](https://doi.org/10.5281/zenodo.23322615)**. This DOI identifies the initial software/data snapshot; subsequent changes to the `main` branch do not alter that archived release.

## Reproducibility scope

The shipped processed CSV includes precomputed regime and zone assignments (`States`, `labels`, `medoids`). The main R script fits conditional station margins and vine copulas on these processed inputs, runs 2022 MJJA simulations, calculates diagnostics and scoring metrics, and exports figures and tables.

**Important limitation:** this is a *processed-data simulation/validation archive*, not a complete raw-data workflow rebuilding the upstream nHMM regime estimation, atmospheric covariate construction, and DTW clustering from their original sources. The script has not been independently rerun in a clean environment as part of this metadata/pre-release preparation. Users should report any runtime or data-compatibility issues.

## Repository structure

~~~text
Data/
  mapa_estaciones_precipitaciones_clean.csv.zip
Scripts/
  Hybrid_nHMM_Copula_Simulation_Validation_2022.R
README.md
CITATION.cff
LICENSE
LICENSES.md
.gitignore
~~~

The input archive is approximately 34 MB compressed. It expands to a large CSV (approximately 244 MB according to the original project documentation); allow adequate disk space and RAM.

Outputs are written to `fig_main/`, `fig_supp/`, and in some cases `supp_figures/`. Generated files and the extracted CSV are ignored by Git.

## How to run

1. Install **R 4.2 or later** and a suitable CRAN mirror. Clone the repository and open a terminal in its root directory.
2. Install the required R packages before running the script (the script also uses `pacman::p_load` to load packages).
3. Run the main script **from the repository root**. Do not set an absolute local path inside the script.

~~~bash
git clone https://github.com/mauricio-herrera/hybrid-nhmm-copula-precipitation-framework.git
cd hybrid-nhmm-copula-precipitation-framework

Rscript -e 'install.packages(c("pacman","dplyr","tidyr","purrr","tibble","data.table","lubridate","ggplot2","stringr","VineCopula","copula","evd","goftest","ADGofTest","scoringRules","forecast","readr","patchwork"))'

Rscript Scripts/Hybrid_nHMM_Copula_Simulation_Validation_2022.R
~~~

The script detects `Data/mapa_estaciones_precipitaciones_clean.csv.zip`, extracts its CSV if needed, and reads `Data/mapa_estaciones_precipitaciones_clean.csv`. Running outside the repository root gives an explicit diagnostic. The full simulation is potentially computationally demanding.

The main script contains multiple analysis sections and repeated initialization to support legacy checks; its numerical implementation is preserved aside from the two input-path initialization blocks. Figures, diagnostics, and tables may require substantial CPU time and memory.

## Input data dictionary

The accompanying CSV is described in the manuscript/project documentation as containing:

| Column | Interpretation |
|---|---|
| `dates`, `year`, `month`, `day`, `Year` | Dates and calendar fields |
| `codigo_estacion`, `location` | Station and auxiliary location identifiers |
| `precipitation` | Daily precipitation (mm) |
| `x`, `y`, `z` | Coordinates and altitude |
| `labels` | Precomputed precipitation-zone assignment |
| `medoids` | DTW-cluster medoid flag |
| `States` | Precomputed hidden rainfall state |
| `year1980` | Historical-time index |

The actual data archive has not been independently decompressed or audited as part of this metadata update; confirm its schema if the script raises a column-name error.

The training window in the script ends **2021-08-31**, the target year is **2022**, and the validation months are **May–August**. The default `selected_zone` is 5 and `M_per_day` is 100, with random seeds specified in the code.

## Source-data attribution

The manuscript identifies **CR2MET** as the gridded precipitation source. Please acknowledge:

Boisier, Juan P. (2023). *CR2MET: A high-resolution precipitation and temperature dataset for the period 1960–2021 in continental Chile* (v2.5) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7529682

CR2 source: https://www.cr2.cl/datos-productos-grillados/

CR2MET is published with Creative Commons 4.0 licensing information. Source-data and downstream derivative rights are distinct from the MIT software license. See [LICENSES.md](LICENSES.md) for usage and attribution notes. Confirm the provenance and redistribution rights for every additional derived input before republishing the data.

## Citation and archival releases

GitHub's **Cite this repository** action uses [CITATION.cff](CITATION.cff). The first public version was archived via the GitHub–Zenodo integration on 12 October 2026.

**Citation:** Herrera-Marín, M. (2026). *Hybrid Hidden Markov–Copula Precipitation Framework* (Version 1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.23322615

For reproducibility of the first release, use its immutable [GitHub tag](https://github.com/mauricio-herrera/hybrid-nhmm-copula-precipitation-framework/tree/v1.0.0) rather than assuming that the latest `main` matches the frozen Zenodo archive. See [docs/ZENODO_RELEASE.md](docs/ZENODO_RELEASE.md) for archival provenance.

## License

Original software: MIT; see [LICENSE](LICENSE). The dataset and any third-party-derived material require their own attribution and license assessment; see [LICENSES.md](LICENSES.md).
