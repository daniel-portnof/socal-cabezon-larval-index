# socal-cabezon-larval-index
Using sdmTMB to derive a larval abundance index from CalCOFI data to be compared with existing federal stock metrics for cabezon in Southern California.

*NOTE:* You will need to download the `Cab_SCS_BC_STAR/` folder, `cabezon_calcofi_data.csv`, and `bottle_data.csv` into your working directory to run the scripts.

### Notes (04-22-2026):
- Currently working from `cabezon_constrained_model.Rmd` and `cabezon_environmental_analysis.Rmd`.
  
- `cabezon_constrained_model.Rmd` roughly runs as follows:
  
  1. Setup and check STAR output plots
  2. CalCOFI data import and QA/QC
  3. Data cuts and visualization of observed patterns
  4. Model construction and comparison
  5. Residual checks and derive abundance index
  6. Test against STAR outputs via CCF
  7. Regression analysis on significant correlation(s)
    
- Currently, there are no major patterns detected. The strongest result is a weak (and biologically implausible) negative relationship between larval abundance and recruitment metrics, in which recruitment metrics *lead* abundance by one year. This is spurious.
- Further investigation also reveals little to no relationship between larval abundance indices and SSB.
- Curious if a finer-grade data source on adult cabezon abundance (e.g., catch if there was a strong enough dataset) would track better with these indices.

- `cabezon_environmental_analysis.Rmd` is checking for how oceanographic data collected via CalCOFI bottle samples (`bottle_data.csv`) relates to larval presence/absence (using the same manta net tow data `cabezon_calcofi_data.csv`) and positive density. It roughly runs as follows:

  1. Setup and data import & QA/QC
  2. Join and standardize bottle and larval data
  3. Run candidate presence/absence GLMMs and compare (station/location as random effect)
  4. Run candidate presence/absence sdm's and compare
  5. Refit best models with hurdle component
  6. Check summaries and residuals
  
## Repository Contents

| File / Folder | Description |
|---|---|
| `cabezon_primary_analysis.Rmd` | All analyses, inclusive of the full and constrained datasets and appendix materials — knit to reproduce report |
| `cabezon_full_model.Rmd` | Analysis on full CalCOFI dataset; all years, all stations — knit to reproduce report |
| `cabezon_constrained_model.Rmd` | Analysis on constrained dataset; spawning season only, coastal shelf stations only — knit to reproduce report |
| `cabezon_environmental_analysis.Rmd` | Analysis of environmental covariates from CalCOFI bottle data and constrained larval dataset; spawning season only, coastal shelf stations only — knit to reproduce report |
| `cabezon_calcofi_data.csv` | CalCOFI larval survey data (1981–2015) |
| `bottle_data.csv` | CalCOFI oceanographic bottle survey data (1954—2021) |
| `Cab_SCS_BC_STAR/` | Stock Synthesis model outputs (STAR panel assessment) |
| `cabezon_appendix.Rmd` | Appendix containing experimental models, code chunks, and other ephemera |
| `cabezon_primary_analysis.R` | Standalone script — **currently outdated** |

## Requirements

**R Version:** 4.4.3 (2025-02-28 ucrt)

### Key packages
- `sdmTMB` 1.0.0
- `tidyverse` 2.0.0 (includes `ggplot2`, `dplyr`, `tidyr`, `readr`, `tibble`, `stringr`, `forcats`, `purrr`)
- `sf` 1.0-23
- `rnaturalearth` / `rnaturalearthdata` 1.1.0 / 1.0.0
- `maps` 3.4.3
- `here` 1.0.2
- `r4ss` 1.44.0
- `lubridate` 1.9.4
