"""
Restructures cabezon_primary_analysis.Rmd into:
  1. cabezon_primary_analysis.Rmd  (parent, all content, renumbered)
  2. cabezon_full_model.Rmd        (Part I: sections 0-7.2)
  3. cabezon_constrained_model.Rmd (Part II: sections 8-10.6)
  4. cabezon_appendix.Rmd          (Part III: A1-A4)
"""

import re

SRC = "cabezon_primary_analysis.Rmd"

with open(SRC, "r") as f:
    lines = f.readlines()

text = "".join(lines)

# ── Locate section boundaries ─────────────────────────────────────────────────
def find(marker):
    for i, line in enumerate(lines):
        if marker in line:
            return i
    raise ValueError(f"Marker not found: {marker!r}")

i_yaml_end    = find("```{r setup")                    # line 8
i_s0          = find("# 0. LOAD IN AND STAR REPORT")   # line 12
i_s1          = find("# 1. DATA QC + EXPLORATION")     # line 105
i_s2          = find("# 2. DATA EXPLORATION")          # line 168
i_s3          = find("# 3. STAR COMPARISONS & INITIAL CCFs")  # line 358
i_s4          = find("# 4. DATA VISUALIZATION")        # line 419
i_s5          = find("# 5. sdmTMB MODEL FITTING")      # line 661
i_s6          = find("# 6. cab_fit1 full D-G MODEL")   # line 931
i_s65         = find("# 6.5. cab_fit01 D-G MODEL")     # line 1059
i_s7          = find("# 7. STAR COMPARISONS & CCFs TO cab_fit1") # line 1189
i_s75         = find("# 7.5. STAR COMPARISONS")        # line 1242
i_rerun       = find("# RE-RUNNING WITH CONSTRAINED")  # line 1295
i_s10         = find("## 10. DATA PARING")             # line 1307
i_s11         = find("## 11. sdmTMB MODEL FITTING (CONSTRAINED)") # line 1563
i_s12         = find("## 12. ss_cab_fit1 (CONSTRAINED) FULL") # line 1801
i_s125        = find("### 12.5 ss_cab_fit1 STAR")      # line 1929
i_s127        = find("### 12.7 ss_cab_fit1 correlation") # line 1980
i_s13         = find("## 13. ss_cab_fit01")            # line 2198
i_s135        = find("### 13.5. ss_cab_fit01")         # line 2325
i_s137        = find("### 13.7 ss_cab_fit01")          # line 2378
i_experiments = find("# Experiments\n")                # line 2526
i_addl        = find("# Additional model testing")     # line 2697
i_s67         = find("```{r 6.7. cab_fit2}")           # line 2699
i_s69         = find("```{r 6.9. cab_fit03}")          # line 2850
i_extras      = find("# Extras\n")                     # line 3023
i_eof         = len(lines)

# ── Section extractors ────────────────────────────────────────────────────────
def get(start, end):
    return "".join(lines[start:end])

yaml_setup = get(0, i_s0)
s0_raw     = get(i_s0, i_s1)
s1_raw     = get(i_s1, i_s2)
s2_raw     = get(i_s2, i_s3)
s3_raw     = get(i_s3, i_s4)
s4_raw     = get(i_s4, i_s5)
s5_raw     = get(i_s5, i_s6)
s6_raw     = get(i_s6, i_s65)
s65_raw    = get(i_s65, i_s7)
s7_raw     = get(i_s7, i_s75)
s75_raw    = get(i_s75, i_rerun)
rerun_raw  = get(i_rerun, i_s10)
s10_raw    = get(i_s10, i_s11)
s11_raw    = get(i_s11, i_s12)
s12_raw    = get(i_s12, i_s125)
s125_raw   = get(i_s125, i_s127)
s127_raw   = get(i_s127, i_s13)
s13_raw    = get(i_s13, i_s135)
s135_raw   = get(i_s135, i_s137)
s137_raw   = get(i_s137, i_experiments)
exp_raw    = get(i_experiments, i_addl)
s67_raw    = get(i_addl, i_s69)   # includes "# Additional model testing" header + 6.7 chunk
s69_raw    = get(i_s69, i_extras)
extras_raw = get(i_extras, i_eof)

# ── Text transformations ──────────────────────────────────────────────────────
def rpl(s, pairs):
    for old, new in pairs:
        s = s.replace(old, new)
    return s

# YAML
yaml_setup = rpl(yaml_setup, [
    ('title: "Cabezon sdmTMB"', 'title: "Cabezon sdmTMB — Primary Analysis"'),
])

# Section 0: rename header + chunk label, add fmesher + segmented to library block
s0 = rpl(s0_raw, [
    ("# 0. LOAD IN AND STAR REPORT VIEW\n", "# 0. SETUP & STAR REPORT\n"),
    ("```{r STAR report}", "```{r 0. SETUP & STAR REPORT}"),
    ("# ---- 0. LOAD IN AND STAR REPORT VIEW\n\n#Load packages\n",
     "# ---- 0. SETUP & STAR REPORT\n\n# Load packages\n"),
    ("library(sdmTMB)\n",
     "library(sdmTMB)\nlibrary(fmesher)\nlibrary(segmented)\n"),
])

# Section 5: remove library(fmesher) (now in s0), remove RE-RUN sticky note
s5 = rpl(s5_raw, [
    ("library(fmesher)\n", ""),
    ("\n**RE-RUN** models constrained to spawning season and spatially\n\n\n", "\n"),
    ("**RE-RUN** models constrained to spawning season and spatially\n\n\n", ""),
    ("**RE-RUN** models constrained to spawning season and spatially\n", ""),
])

# Section 6 → 6.1
s61 = rpl(s6_raw, [
    ("# 6. cab_fit1 full D-G MODEL DIAGNOSTICS & ABUNDANCE INDEX\n",
     "## 6.1. cab_fit1 — Full Delta-Gamma Model\n"),
    ("```{r 6. cab_fit1 full D-G MODEL DIAGNOSTICS & ABUNDANCE INDEX}",
     "```{r 6.1. cab_fit1 diagnostics}"),
    ("#---- 6. cab_fit1 full D-G MODEL DIAGNOSTICS & ABUNDANCE INDEX\n\n", ""),
    ("#---- 6.a.", "#---- 6.1.a."),
    ("#---- 6.b.", "#---- 6.1.b."),
    ("#---- 6.c.", "#---- 6.1.c."),
    ("#---- 6.d.", "#---- 6.1.d."),
    ("#---- 6.e.", "#---- 6.1.e."),
    ("#---- 6.f.", "#---- 6.1.f."),
    ("but 4.f. (observed station grid)", "but 6.1.f. (observed station grid)"),
])

# Section 6.5 → 6.2
s62 = rpl(s65_raw, [
    ("# 6.5. cab_fit01 D-G MODEL DIAGNOSTICS & ABUNDANCE INDEX\n",
     "## 6.2. cab_fit01 — Fixed Year Effects Delta-Gamma Model\n"),
    ("```{r 6.5. cab_fit01 D-G MODEL DIAGNOSTICS & ABUNDANCE INDEX}",
     "```{r 6.2. cab_fit01 diagnostics}"),
    ("#---- 6.5. cab_fit01 D-G MODEL DIAGNOSTICS & ABUNDANCE INDEX\n\n", ""),
    ("#---- 6.5.a.", "#---- 6.2.a."),
    ("#---- 6.5.b.", "#---- 6.2.b."),
    ("#---- 6.5.c.", "#---- 6.2.c."),
    ("#---- 6.5.d.", "#---- 6.2.d."),
    ("#---- 6.5.e.", "#---- 6.2.e."),
    ("#---- 6.5.f.", "#---- 6.2.f."),
    ("but 4.f. (observed station grid)", "but 6.2.f. (observed station grid)"),
])

# Section 7 → 7.1
s71 = rpl(s7_raw, [
    ("# 7. STAR COMPARISONS & CCFs TO cab_fit1 full D-G MODEL RUN\n",
     "## 7.1. cab_fit1\n"),
    ("```{r 7. STAR COMPARISONS & CCFs TO cab_fit1 full D-G MODEL RUN}",
     "```{r 7.1. STAR comparisons cab_fit1}"),
    ("#---- 7. STAR COMPARISONS & CCFs TO cab_fit1 full D-G MODEL RUN\n\n", ""),
    ("constrained 6.f. abundnance index", "constrained 6.1.f. abundance index"),
])

# Section 7.5 → 7.2
s72 = rpl(s75_raw, [
    ("# 7.5. STAR COMPARISONS & CCFs TO cab_fit01 D-G MODEL RUN\n",
     "## 7.2. cab_fit01\n"),
    ("```{r 7.5. STAR COMPARISONS & CCFs TO cab_fit01 D-G MODEL RUN}",
     "```{r 7.2. STAR comparisons cab_fit01}"),
    ("#---- 7.5. STAR COMPARISONS & CCFs TO cab_fit01 D-G MODEL RUN\n\n", ""),
    ("constrained 6.f. abundnance index", "constrained 6.2.f. abundance index"),
])

# RE-RUNNING prose: keep as Part II intro, fix nested header
rerun = rpl(rerun_raw, [
    ("## 10. DATA PARING AND EXPLORATION\n", ""),  # will be added as # 8 below
])

# Section 10 → 8
s8 = rpl(s10_raw, [
    ("## 10. DATA PARING AND EXPLORATION\n", "# 8. DATA PARING & EXPLORATION (CONSTRAINED)\n"),
    ("```{r 10. DATA PARING AND EXPLORATION}", "```{r 8. DATA PARING AND EXPLORATION (CONSTRAINED)}"),
    ("# ---- 10.b.", "# ---- 8.b."),
    ("# ---- 10.c.", "# ---- 8.c."),
    ("#---- 4.d.", "#---- 8.d."),
    ("#---- 4.e.", "#---- 8.e."),
])

# Section 11 → 9: remove library(fmesher), rename
s9 = rpl(s11_raw, [
    ("## 11. sdmTMB MODEL FITTING (CONSTRAINED)\n", "# 9. sdmTMB MODEL FITTING (CONSTRAINED)\n"),
    ("```{r 11. sdmTMB MODEL FITTING (CONSTRAINED)}", "```{r 9. sdmTMB MODEL FITTING (CONSTRAINED)}"),
    ("#---- 11.a.", "#---- 9.a."),
    ("#---- 11.b.", "#---- 9.b."),
    ("#---- 11.c.", "#---- 9.c."),
    ("#---- 11.d.", "#---- 9.d."),
    ("library(fmesher)\n", ""),
])

# Section 12 → 10.1
s101 = rpl(s12_raw, [
    ("## 12. ss_cab_fit1 (CONSTRAINED) FULL D-G MODEL DIAGNOSTICS & ABUNDANCE INDEX\n",
     "## 10.1. ss_cab_fit1 — Full Delta-Gamma Model\n"),
    ("```{r 12. ss_cab_fit1 (CONSTRAINED) FULL D-G MODEL DIAGNOSTICS & ABUNDANCE INDEX}",
     "```{r 10.1. ss_cab_fit1 diagnostics}"),
    ("#---- 12. ss_cab_fit1 (CONSTRAINED) FULL D-G MODEL DIAGNOSTICS & ABUNDANCE INDEX\n\n", ""),
    ("#---- 12.a.", "#---- 10.1.a."),
    ("#---- 12.b.", "#---- 10.1.b."),
    ("#---- 12.c.", "#---- 10.1.c."),
    ("#---- 12.d.", "#---- 10.1.d."),
    ("#---- 12.e.", "#---- 10.1.e."),
    ("#---- 12.f.", "#---- 10.1.f."),
    ("but 4.f. (observed station grid)", "but 10.1.f. (observed station grid)"),
])

# Section 12.5 → 10.2
s102 = rpl(s125_raw, [
    ("### 12.5 ss_cab_fit1 STAR comparisons and CCF\n",
     "## 10.2. ss_cab_fit1 — STAR Comparisons & CCFs\n"),
    ("```{r 12.5 ss_cab_fit1 STAR comparisons and CCF}",
     "```{r 10.2. ss_cab_fit1 STAR comparisons}"),
    ("#---- 12.5 ss_cab_fit1 STAR comparisons and CCF\n",
     "#---- 10.2. ss_cab_fit1 STAR comparisons and CCF\n"),
    ("constrained 6.f. abundnance index", "constrained 10.1.f. abundance index"),
])

# Section 12.7 → 10.3: remove library(segmented)
s103 = rpl(s127_raw, [
    ("### 12.7 ss_cab_fit1 correlation testing\n",
     "## 10.3. ss_cab_fit1 — Correlation Testing\n"),
    ("```{r 12.7 ss_cab_fit1 correlation testing}",
     "```{r 10.3. ss_cab_fit1 correlation testing}"),
    ("library(segmented)\n\n", ""),
    ("library(segmented)\n", ""),
    ("# 3. Compare AIC to ss_cab_fit1c (Gamma AR1 only) \n# which had cleaner diagnostics in your original analysis\n",
     "# 3. AIC comparison: ss_cab_fit1 vs. ss_cab_fit1c (Gamma AR1 only)\n"),
])

# Section 13 → 10.4
s104 = rpl(s13_raw, [
    ("## 13. ss_cab_fit01 (CONSTRAINED) D-G MODEL FIXED YEAR EFFECTS DIAGNOSTICS & ABUNDANCE INDEX\n",
     "## 10.4. ss_cab_fit01 — Fixed Year Effects Delta-Gamma Model\n"),
    ("```{r 13. ss_cab_fit01 (CONSTRAINED) D-G MODEL FIXED YEAR EFFECTS DIAGNOSTICS & ABUNDANCE INDEX}",
     "```{r 10.4. ss_cab_fit01 diagnostics}"),
    ("#---- 13.a.", "#---- 10.4.a."),
    ("#---- 13.b.", "#---- 10.4.b."),
    ("#---- 13.c.", "#---- 10.4.c."),
    ("#---- 13.d.", "#---- 10.4.d."),
    ("#---- 13.e.", "#---- 10.4.e."),
    ("#---- 13.f.", "#---- 10.4.f."),
])

# Section 13.5 → 10.5: fix ss_cab_fit1 → ss_cab_fit01 in comment
s105 = rpl(s135_raw, [
    ("### 13.5. ss_cab_fit01 STAR comparisons and CCF\n",
     "## 10.5. ss_cab_fit01 — STAR Comparisons & CCFs\n"),
    ("```{r 13.5. STAR comparisons and CCF}",
     "```{r 10.5. ss_cab_fit01 STAR comparisons}"),
    ("#---- 13.5 ss_cab_fit1 STAR comparisons and CCF\n",
     "#---- 10.5. ss_cab_fit01 STAR comparisons and CCF\n"),
    ("constrained 6.f. abundnance index", "constrained 10.4.f. abundance index"),
])

# Section 13.7 → 10.6
s106 = rpl(s137_raw, [
    ("### 13.7 ss_cab_fit01 correlation testing\n",
     "## 10.6. ss_cab_fit01 — Correlation Testing\n"),
    ("```{r 13.7 ss_cab_fit01 correlation testing}",
     "```{r 10.6. ss_cab_fit01 correlation testing}"),
    ("# 3. Cross-model AIC comparison: AR1 vs fixed year effect (constrained dataset)\n",
     "# 3. AIC comparison: AR1 vs. fixed year effect (constrained dataset)\n"),
])

# Experiments → A3: remove library(segmented)
s_a3 = rpl(exp_raw, [
    ("# Experiments\n", "# A3. EXPERIMENTS\n"),
    ("```{r Experiments}", "```{r A3. Experiments}"),
    ("library(segmented)\n\n", ""),
    ("library(segmented)\n", ""),
])

# 6.7 → A1
s_a1 = rpl(s67_raw, [
    ("# Additional model testing\n", "# A1. TWEEDIE MODEL — cab_fit2\n"),
    ("```{r 6.7. cab_fit2}", "```{r A1. cab_fit2 Tweedie}"),
    ("#---- 6.7. cab_fit2 TWEEDIE MODEL DIAGNOSTICS & ABUNDANCE INDEX\n",
     "#---- A1. cab_fit2 Tweedie — Diagnostics & Abundance Index\n"),
    ("#---- 6.7.", "#---- A1."),
    ("#----7.7.", "#---- A1."),
    ("#----7.7.a.", "#---- A1.STAR.a."),
    ("#---- 7.7.b.", "#---- A1.STAR.b."),
])

# 6.9 → A2
s_a2 = rpl(s69_raw, [
    ("```{r 6.9. cab_fit03}", "```{r A2. cab_fit03}"),
    ("#---- 6.9. cab_fit03 DG FIXED YEAR / IID MODEL DIAGNOSTICS & ABUNDANCE INDEX\n",
     "#---- A2. cab_fit03 DG Fixed Year / IID — Diagnostics & Abundance Index\n"),
    ("#---- 6.9.", "#---- A2."),
    ("#---- 7.9.", "#---- A2.STAR."),
    ("#---- 7.9.a.", "#---- A2.STAR.a."),
    ("#---- 7.9.b.", "#---- A2.STAR.b."),
])

# Extras → A4
s_a4 = rpl(extras_raw, [
    ("# Extras\n", "# A4. EXTRAS\n"),
    ("```{r Extras, error = TRUE}", "```{r A4. Extras, error = TRUE}"),
])

# ── Assemble Part dividers ────────────────────────────────────────────────────
part1_header = "\n---\n\n# PART I — FULL DATASET ANALYSIS\n\n---\n\n"
part2_header = "\n---\n\n# PART II — CONSTRAINED ANALYSIS (SPAWNING SEASON + SHELF REGION)\n\n---\n\n"
part3_header = "\n---\n\n# PART III — APPENDIX / ADDITIONAL MODELS\n\n---\n\n"
s6_header    = "\n# 6. MODEL DIAGNOSTICS & ABUNDANCE INDICES\n\n"
s7_header    = "\n# 7. STAR COMPARISONS & CCFs\n\n"
s10_header   = "\n# 10. MODEL DIAGNOSTICS & ABUNDANCE INDICES (CONSTRAINED)\n\n"

# ── Write PARENT file ─────────────────────────────────────────────────────────
parent = "".join([
    yaml_setup,
    part1_header,
    s0, s1_raw, s2_raw, s3_raw, s4_raw, s5,
    s6_header, s61, s62,
    s7_header, s71, s72,
    part2_header,
    rerun, s8, s9,
    s10_header, s101, s102, s103, s104, s105, s106,
    part3_header,
    s_a1, s_a2, s_a3, s_a4,
])

with open("cabezon_primary_analysis.Rmd", "w") as f:
    f.write(parent)
print("✓ cabezon_primary_analysis.Rmd written")

# ── Write FULL MODEL file ─────────────────────────────────────────────────────
full_yaml = yaml_setup.replace(
    'title: "Cabezon sdmTMB — Primary Analysis"',
    'title: "Cabezon sdmTMB — Full Dataset Analysis"'
)

full_model = "".join([
    full_yaml,
    part1_header,
    s0, s1_raw, s2_raw, s3_raw, s4_raw, s5,
    s6_header, s61, s62,
    s7_header, s71, s72,
])

with open("cabezon_full_model.Rmd", "w") as f:
    f.write(full_model)
print("✓ cabezon_full_model.Rmd written")

# ── Write CONSTRAINED MODEL file ──────────────────────────────────────────────
constrained_yaml = yaml_setup.replace(
    'title: "Cabezon sdmTMB — Primary Analysis"',
    'title: "Cabezon sdmTMB — Constrained Analysis"'
)

# Self-contained setup note for Part II
constrained_setup = """\n```{r constrained setup, include=FALSE}
# ---- Setup for constrained analysis
# Packages and STAR-derived objects are carried over from cabezon_full_model.Rmd.
# Re-load everything here so this file runs independently.
library(tidyverse)
library(r4ss)
library(here)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)
library(sdmTMB)
library(fmesher)
library(segmented)

my.seed <- 666
set.seed(my.seed)

wd <- "Cab_SCS_BC_STAR"
pp <- SS_output(wd)
derived_quants <- pp[["derived_quants"]]

devyrs  <- 1970:2018
devdf   <- pp$parameters[28:76, ]
devvalue <- devdf[, 3]; devsd <- devdf[, 11]
STAR_recdevs <- data.frame(year = devyrs, value = devvalue,
                            lo = devvalue - 1.96 * devsd,
                            hi = devvalue + 1.96 * devsd)

age0yrs <- 1970:2018
age0    <- derived_quants[174:222, "Value"]
age0sd  <- derived_quants[174:222, "StdDev"]
STAR_age0 <- data.frame(year = age0yrs, value = age0,
                         lo = age0 - 1.96 * age0sd,
                         hi = age0 + 1.96 * age0sd)

ssbyrs <- 1916:2018
ssb    <- derived_quants[3:105, "Value"]
ssbsd  <- derived_quants[3:105, "StdDev"]
STAR_SSB <- data.frame(year = ssbyrs, value = ssb,
                        lo = ssb - 1.96 * ssbsd,
                        hi = ssb + 1.96 * ssbsd)

cabezon <- read.csv("cabezon_calcofi_data.csv") %>%
  mutate(season = as.factor(season),
         station_id = paste(line, station, sep = "_"))

cabezon <- sdmTMB::add_utm_columns(cabezon, ll_names = c("longitude", "latitude"))
cabezon$fyear <- as.factor(cabezon$year)
```\n
"""

constrained = "".join([
    constrained_yaml,
    constrained_setup,
    part2_header,
    rerun, s8, s9,
    s10_header, s101, s102, s103, s104, s105, s106,
])

with open("cabezon_constrained_model.Rmd", "w") as f:
    f.write(constrained)
print("✓ cabezon_constrained_model.Rmd written")

# ── Write APPENDIX file ───────────────────────────────────────────────────────
appendix_yaml = yaml_setup.replace(
    'title: "Cabezon sdmTMB — Primary Analysis"',
    'title: "Cabezon sdmTMB — Appendix & Additional Models"'
)

appendix_setup = """\n```{r appendix setup, include=FALSE}
# ---- Setup for appendix
# All objects below are originally produced in cabezon_full_model.Rmd.
# Re-load here so this file runs independently.
library(tidyverse)
library(r4ss)
library(here)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)
library(sdmTMB)
library(fmesher)
library(segmented)

my.seed <- 666
set.seed(my.seed)

wd <- "Cab_SCS_BC_STAR"
pp <- SS_output(wd)
derived_quants <- pp[["derived_quants"]]

devyrs  <- 1970:2018
devdf   <- pp$parameters[28:76, ]
devvalue <- devdf[, 3]; devsd <- devdf[, 11]
STAR_recdevs <- data.frame(year = devyrs, value = devvalue,
                            lo = devvalue - 1.96 * devsd,
                            hi = devvalue + 1.96 * devsd)

age0yrs <- 1970:2018
age0    <- derived_quants[174:222, "Value"]
age0sd  <- derived_quants[174:222, "StdDev"]
STAR_age0 <- data.frame(year = age0yrs, value = age0,
                         lo = age0 - 1.96 * age0sd,
                         hi = age0 + 1.96 * age0sd)

ssbyrs <- 1916:2018
ssb    <- derived_quants[3:105, "Value"]
ssbsd  <- derived_quants[3:105, "StdDev"]
STAR_SSB <- data.frame(year = ssbyrs, value = ssb,
                        lo = ssb - 1.96 * ssbsd,
                        hi = ssb + 1.96 * ssbsd)

cabezon <- read.csv("cabezon_calcofi_data.csv") %>%
  mutate(season = as.factor(season),
         station_id = paste(line, station, sep = "_"))

cabezon <- sdmTMB::add_utm_columns(cabezon, ll_names = c("longitude", "latitude"))
cabezon$fyear <- as.factor(cabezon$year)

# Mesh (needed for model predictions)
mesh <- fmesher::fm_mesh_2d(loc = cabezon[, c("X","Y")],
                             cutoff = 20, max.edge = c(75, 150), offset = c(45, 120))
cab_mesh <- make_mesh(data = cabezon, c("X", "Y"), mesh = mesh)

# NOTE: Models cab_fit1, cab_fit2, cab_fit03 and prediction objects (pred_obs,
# combined_obs, etc.) must be loaded or re-fit before running this file.
# Source cabezon_full_model.Rmd or load a saved .RData workspace to proceed.
```\n
"""

appendix = "".join([
    appendix_yaml,
    appendix_setup,
    part3_header,
    s_a1, s_a2, s_a3, s_a4,
])

with open("cabezon_appendix.Rmd", "w") as f:
    f.write(appendix)
print("✓ cabezon_appendix.Rmd written")

print("\nAll done.")
