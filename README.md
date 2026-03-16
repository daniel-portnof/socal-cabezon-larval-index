# socal-cabezon-larval-index
Using sdmTMB to derive a larval abundance index from CalCOFI data to be compared with existing federal stock metrics for cabezon (*Scorpaenichthys marmoratus*) in Southern California.

## Overview

This analysis uses sdmTMB to fit several models to CalCOFI ichthyoplankton survey data (1981-2015) and use best/most parsimonious fit to generate a standardized larval abundance index. The index is then evaluated against STAR outputs (specifically recruitment deviations, age-0 recruits, and SSB) to assess utility as an alternative, novel predictor of stock health indicators.

### Notes:
- I am not as familiar with cross-validation techniques as I am with AIC, so I ran the model fit comparisons (sections 6 + 7) on the AIC best fit model (cab_fit1, Delta-Gamma full AR1 spatial + spatiotemporal). I intend to run the CV best fit (cab_fit2, Tweedie full AR1 spatial + spatiotemporal) as well, but am curious what others would suggest.
- On a similar note, these are all model fits that I built myself with the help of an LLM. They may not reflect the best construction, and I look to others for opinions on how to better optimize the models. For example, I had assumed an AR1 fit would be appropriate, but Brice and I also concluded that if spatiotemporal variation is simply not as strong a predictor, perhaps a better approach would be to model a spatial-only fit and extend the time series since temporal effects will be less germane. This would perhaps be a next step if, as the analysis so far would suggest, there is no significant correlation revealed between the novel constructed index and the STAR outputs.

## Repository Contents

| File / Folder | Description |
|---|---|
| `cabezon_primary_analysis.Rmd` | Full analysis with narrative — knit to reproduce report |
| `cabezon_primary_analysis.R` | Same code as .Rmd, for running as a standalone script |
| `cabezon_calcofi_data.csv` | CalCOFI larval survey data (1981–2015) |
| `Cab_SCS_BC_STAR/` | Stock Synthesis model outputs (STAR panel assessment) |

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

Install all at once:
```r
install.packages(c("sdmTMB", "tidyverse", "sf", "rnaturalearth", 
                   "rnaturalearthdata", "maps", "here", "r4ss", "lubridate"))
```

## How to Run

1. Clone or download this repo
2. Set your working directory to the root of the repo (the folder containing the `.R`, `.Rmd`, and `.csv` files, as well as `Cab_SCS_BC_STAR/`)
3. Open either `cabezon_primary_analysis.R` (run as script) or `cabezon_primary_analysis.Rmd` (knit to document)
4. The script reads Stock Synthesis outputs directly from `Cab_SCS_BC_STAR/` using the `r4ss` package. Shouldn't need any additional setup beyond setting the working directory.

## Data Sources

**CalCOFI ichthyoplankton survey:** `cabezon_calcofi_data.csv`
- Thanks to Andrew for pulling these data!

**STAR Stock Assessment:** `Cab_SCS_BC_STAR/` 
- SS outputs from the STAR panel assessment for Southern California.
- Thanks to Julia for pulling this!
