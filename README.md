# socal-cabezon-larval-index
Using sdmTMB to derive a larval abundance index from CalCOFI data to be compared with existing federal stock metrics for cabezon in Southern California.

### Notes (04-15-2026):
- `cabezon_primary_analysis.Rmd` is typically the most up-to-date file.
- The analysis roughly runs as follows
  1. Check STAR output plots
  2. CalCOFI data import and QA/QC
  3. Visualization of patterns in observed data
  4. Model construction and comparison
  5. Derive abundance index
  6. Test against STAR outputs via CCF
  7. Regression analysis on significant correlation(s)
- Currently, there are no major patterns detected. The strongest result is a weak (and biologically implausible) negative relationship between larval abundance and recruitment metrics, in which recruitment metrics *lead* abundance by one year. This is likely spurious.
  - It is also worth noting that the indices used for comparison are drawn only on the discrete CalCOFI station locations rather than the entire interpolated spatial field. While this does artificially reduce uncertainty in the resulting indices, it is an exercise intended only for the sake of initial comparisons. If a relationship is detected, further examination will be needed against the models' full spatial predictions.
- Next steps will be to further investigate any relationship between novel indices and SSB and to validate best-fit models against environmental covariates.

## Repository Contents

| File / Folder | Description |
|---|---|
| `cabezon_primary_analysis.Rmd` | All analyses, inclusive of the full and constrained datasets and appendix materials — knit to reproduce report |
| `cabezon_full_model.Rmd` | Analysis on full CalCOFI dataset; all years, all stations — knit to reproduce report |
| `cabezon_constrained_model.Rmd` | Analysis on constrained dataset; spawning season only, coastal shelf stations only — knit to reproduce report |
| `cabezon_calcofi_data.csv` | CalCOFI larval survey data (1981–2015) |
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
