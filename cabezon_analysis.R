# -==============================================================================-
# ==== Southern California Cabezon Larval Abundance Index — Analysis Pipeline ====
# -==============================================================================-

# End-to-end analysis for the constrained (spawning-season, SCB-region) CalCOFI cabezon larval index: data load, model fitting, diagnostics, abundance-index construction, stock-assessment comparisons, publication figures, and the environmental covariate analysis.

# The script runs top to bottom. Sections 0-8 are the canonical pipeline that produces the reported results and figures. The Appendix holds supplemental, sensitivity, and alternative-model analyses; nothing in Sections 0-8 depends on it.

# =-------------------------------------------------------------------=
## ---- Expected repo layout ----
# =-------------------------------------------------------------------=
#   Data/
#     cabezon_calcofi_data.csv        larval tow records (CalCOFI)
#     bottle_data.csv                 CalCOFI hydrographic bottle casts
#     CalCOFI_75StandardStations.kml  standard-station grid
#   Figures/                          output figures and tables
#   Cab_SCS_BC_STAR/                  Stock Synthesis assessment outputs
#   cabezon_analysis.R                this script

# =-------------------------------------------------------------------=
## ---- Contents ----
# =-------------------------------------------------------------------=
#   0.  Setup & stock-assessment (STAR) reference series
#   1.  Data QC & exploration
#   2.  Data paring (spawning season + SCB) & exploration
#   3.  sdmTMB model fitting (candidate index models)
#   4.  Preferred-model diagnostics & abundance index
#   5.  Stock-assessment comparisons & cross-correlations
#   6.  Regression analyses (index vs. recruitment)
#   7.  Publication figures & tables
#   8.  Environmental covariate analysis
#   A.  Appendix (supplemental / sensitivity / alternative models)
#
# Preferred index model: ss_cab_fit_dln3 (delta-lognormal, spatial field on, IID encounter spatiotemporal field; "DLn-IID").

# =-------------------------------------------------------------------=
## ---- Dependencies ----
# =-------------------------------------------------------------------=
library(tidyverse)
library(r4ss)
library(sdmTMB)
library(fmesher)
library(funtimes)
library(lme4)
library(corrplot)
library(DHARMa)
library(visreg)
library(patchwork)
library(flextable)
library(scales)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(marmap)
library(viridis)

# -==============================================================================-
# ==== 0.  SETUP & STOCK-ASSESSMENT (STAR) REFERENCE SERIES ====
# -==============================================================================-

# ---- 0. SETUP & STAR REPORT (Part II)

# Load packages
library(tidyverse)
library(r4ss)
library(sdmTMB)
library(fmesher)
library(funtimes)

my.seed <- 666
set.seed(my.seed)

wd <- "Cab_SCS_BC_STAR" # Mac

# This reads the output files
pp <- SS_output(wd)

# Pull out derived quantities.  These are things like SSB and age-0 recruits.
derived_quants <- pp[["derived_quants"]]

# Pull out parameter estimates.  Recruitment deviations are estimated parameters.
params <- pp[["parameters"]]

# Plot recruitment deviations for the main estimation period.
devyrs <- 1970:2018
devdf <- pp$parameters[28:76, ]
devvalue <- devdf[,3]
devsd <- devdf[,11]
devlower <- (devvalue-1.96*devsd)
devmed   <- devvalue
devupper <- (devvalue+1.96*devsd)
STAR_recdevs <- data.frame(year=devyrs, value=devmed,lo=devlower,hi=devupper)

STAR_recdevs_plot <- ggplot(data=STAR_recdevs, aes(x=year, ymin=devlower, ymax=devupper)) +
  geom_ribbon(fill="skyblue", alpha=0.5) +
  geom_line(aes(y=devmed)) +
  labs(title="Cabezon SoCal Recruitment Deviations") +
  theme_bw()

STAR_recdevs_plot

# Plot age-0 recruits
age0yrs <- 1970:2018
startrowindex <- which(rownames(derived_quants)=="Recr_1970")
endrowindex <- which(rownames(derived_quants)=="Recr_2018")
age0 <- derived_quants[174:222, "Value"]
age0sd <- derived_quants[174:222, "StdDev"]
age0lower <- (age0-1.96*age0sd)
age0upper <- (age0+1.96*age0sd)

STAR_age0 <- data.frame(year=age0yrs, value=age0,lo=age0lower,hi=age0upper)

STAR_age0_plot <- ggplot(data=STAR_age0, aes(x=year, ymin=age0lower, ymax=age0upper)) +
  geom_ribbon(fill="skyblue", alpha=0.5) +
  geom_line(aes(y=age0)) +
  labs(title="Cabezon SoCal Age 0 Recruits") +
  theme_bw()

STAR_age0_plot

# Plot SSB
ssbyrs <- 1916:2018
startrowindex <- which(rownames(derived_quants)=="SSB_1916")
endrowindex <- which(rownames(derived_quants)=="SSB_2018")
ssb <- derived_quants[3:105, "Value"]
ssbsd <- derived_quants[3:105, "StdDev"]
ssblower <- (ssb-1.96*ssbsd)
ssbupper <- (ssb+1.96*ssbsd)

STAR_SSB <- data.frame(year=ssbyrs, value=ssb,lo=ssblower,hi=ssbupper)

STAR_SSB_plot <- ggplot(data=STAR_SSB, aes(x=year, ymin=ssblower, ymax=ssbupper)) +
  geom_ribbon(fill="skyblue", alpha=0.5) +
  geom_line(aes(y=ssb)) +
  labs(title="Cabezon SoCal SSB") +
  theme_bw()

STAR_SSB_plot


# -==============================================================================-
# ==== 1.  DATA QC & EXPLORATION ====
# -==============================================================================-

#---- DATA QC + EXPLORATION (Part II)

cabezon <- read.csv("data/cabezon_calcofi_data.csv")

str(cabezon)

# Add factor variable for season
cabezon <- cabezon %>%
  mutate(season = as.factor(season))

# ---- 1. DATA QC + EXPLORATION (Part II)

# Check for how many observations per ...

# year
cabezon %>%
  count(year, season, line.station) %>%
  filter(n > 1)
# month
cabezon %>%
  count(year, month, line.station) %>%
  filter(n > 1)
# cruise
cabezon %>%
  count(cruise, line.station) %>%
  filter(n > 1)

# Station coverage across years
cabezon %>%
  group_by(year) %>%
  summarize(n_stations = n_distinct(line.station))
table(cabezon$line.station, cabezon$year)
# Pretty steady coverage, especially post-1989

# Check the spatial grid for a quick visual
ggplot(cabezon, aes(longitude, latitude)) +
  geom_point(alpha = 0.3)

# How much zero inflation do we have?
cabezon %>%
  mutate(zero = larvae_100m3 == 0) %>%
  count(zero)
# Roughly ~85% zeros

# Raw distribution visual
plot(cabezon$volume_sampled, cabezon$larvae_100m3)

pos_vals <- cabezon$larvae_100m3[cabezon$larvae_100m3 > 0]
hist(log(pos_vals), breaks = 50, main = "log(positive larvae_100m3)")
hist(log(cabezon$larvae_100m3[cabezon$larvae_100m3 > 0]))

cabezon %>%
  filter(larvae_100m3 > 0) %>%
  group_by(year) %>%
  summarize(m = mean(larvae_100m3), v = var(larvae_100m3), n = n()) %>%
  filter(n > 5) %>%
  ggplot(aes(log(m), log(v))) + geom_point() + geom_smooth(method = "lm") +
  geom_abline(slope = 2, intercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(title = "Mean-variance relationship (dashed: red=gamma slope 2, blue=Poisson slope 1)")

# Also useful: what fraction are zeros?
mean(cabezon$larvae_100m3 == 0)


# -==============================================================================-
# ==== 2.  DATA PARING (SPAWNING SEASON + SCB) & EXPLORATION ====
# -==============================================================================-


cabezon_shelf_spawn <- cabezon %>%
  filter(
    month %in% c(1, 2, 3),
    station <= 60.0
  )

cat("Rows retained:", nrow(cabezon_shelf_spawn), "of", nrow(cabezon), "\n")
cat("Years covered:", paste(range(cabezon_shelf_spawn$year), collapse = "–"), "\n")
cat("Months present:", paste(sort(unique(cabezon_shelf_spawn$month)), collapse = ", "), "\n")
cat("Station range:", paste(range(cabezon_shelf_spawn$station), collapse = "–"), "\n")
cat("Lines present:", paste(sort(unique(cabezon_shelf_spawn$line)), collapse = ", "), "\n")

cabezon_shelf_spawn %>%
  group_by(year) %>%
  summarize(
    n_tows = n(),
    n_pos = sum(larvae_100m3 > 0),
    pct_pos = round(100 * n_pos / n_tows, 1)
  ) %>%
  print(n = Inf)

# Composite station column (combines line + station)
cabezon_shelf_spawn <- cabezon_shelf_spawn %>%
  mutate(station_id = paste(line, station, sep = "_"))

# Do non-standard hauls have systematically different standard_haul_factor values?
summary(cabezon_shelf_spawn$standard_haul_factor)

# Are any values notably far from 1.0?
hist(cabezon_shelf_spawn$standard_haul_factor)
hist(log(cabezon_shelf_spawn$larvae_100m3))

# Check the spike at larvae_100m3 == 1
cabezon_shelf_spawn %>% count(larvae_100m3 == 1)

# Check how many true zeros are present
cabezon_shelf_spawn %>% count(larvae_100m3 == 0)

# Confirm no -Inf values making it into model inputs
sum(is.infinite(log(cabezon_shelf_spawn$larvae_100m3)))

cabezon_shelf_spawn %>%
  filter(standard_haul_factor != 0.02) %>%
  count(year) 

# Calculating positive tow values for visualization
ss_pos_tows_year <- cabezon_shelf_spawn %>%
  group_by(year) %>%
  summarise(
    total_tows        = n(),
    positive_tows     = sum(larvae_100m3 > 0),
    pct_positive_tows = positive_tows / total_tows * 100
  )

ss_pos_tows_month <- cabezon_shelf_spawn %>%
  group_by(month) %>%
  summarise(
    total_tows        = n(),
    positive_tows     = sum(larvae_100m3 > 0),
    pct_positive_tows = positive_tows / total_tows * 100
  )

ss_pos_tows_season <- cabezon_shelf_spawn %>%
  group_by(season) %>%
  summarise(
    total_tows        = n(),
    positive_tows     = sum(larvae_100m3 > 0),
    pct_positive_tows = positive_tows / total_tows * 100
  )

ss_pos_tows_station <- cabezon_shelf_spawn %>%
  group_by(station_id) %>%
  summarise(
    total_tows        = n(),
    positive_tows     = sum(larvae_100m3 > 0),
    pct_positive_tows = positive_tows / total_tows * 100
  )

# Annual mean larval observations with SE
annual_cab_ss_summary <- cabezon_shelf_spawn %>%
  group_by(year) %>%
  summarize(
    mean_larvae = mean(larvae_100m3, na.rm = T),
    sd_larvae = sd(larvae_100m3, na.rm = T),
    n = n(),
    se = sd_larvae / sqrt(n)
  )

# ---- 8.b. Annual average larval observations

ggplot(annual_cab_ss_summary, aes(x = year, y = mean_larvae)) +
  geom_ribbon(aes(ymin = mean_larvae - se, ymax = mean_larvae + se),
              alpha = 0.5, fill = "skyblue") +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue") +
  labs(
    x = "Year",
    y = expression("Mean larvae per 100 m"^3),
    title = "Annual Mean Larval Observations, Constrained to Spawning Season and Coastal Shelf",
    subtitle = expression("Ribbon shows " %+-% " 1 SE")
  ) +
  theme_bw()

# ---- 8.c. Stacked bar chart: positive vs. total tows by year, month, season

ss_plot_data <- bind_rows(
  ss_pos_tows_year   %>% mutate(group_var = "Year",   level = as.character(year)),
  ss_pos_tows_month  %>% mutate(group_var = "Month",  level = as.character(month)),
  ss_pos_tows_season %>% mutate(group_var = "Season", level = as.character(season))
) %>%
  mutate(negative_tows = total_tows - positive_tows) %>%
  pivot_longer(
    cols      = c(positive_tows, negative_tows),
    names_to  = "tow_type",
    values_to = "count"
  )

ss_plot_data <- ss_plot_data %>%
  mutate(level = case_when(
    group_var == "Month" ~ factor(level, levels = as.character(1:12)),
    group_var == "Season" ~ factor(level, levels = c("winter", "spring", "summer", "fall")),
    T ~ factor(level)
  ))

# Shared elements
shared_fill <- scale_fill_manual(
  values = c("negative_tows" = "grey80", "positive_tows" = "steelblue"),
  labels = c("negative_tows" = "Larvae Absent", "positive_tows" = "Larvae Present")
)

shared_bar_theme <- theme_minimal(base_size = 12) +
  theme(
    strip.text         = element_text(face = "bold", size = 13),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    legend.position    = "bottom",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold", size = 15),
    plot.subtitle      = element_text(color = "grey40", size = 10)
  )

# One plot for years; one plot for months + seasons
ss_plot_data_year <- ss_plot_data %>% filter(group_var == "Year")
ss_plot_data_ms   <- ss_plot_data %>% filter(group_var %in% c("Month", "Season"))

# Years plot
ggplot(ss_plot_data_year, aes(x = level, y = count, fill = tow_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  shared_fill +
  shared_bar_theme +
  labs(
    title    = "Cabezon Larvae: Positive vs. Total Tows",
    subtitle = "Stacked bars show larvae-present (blue) against total sampling effort (gray)",
    x        = NULL,
    y        = "Number of Tows",
    fill     = NULL
  )

# Months + Seasons plot
ggplot(ss_plot_data_ms, aes(x = level, y = count, fill = tow_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~ group_var, scales = "free_x", nrow = 1) +
  shared_fill + shared_bar_theme +
  labs(
    title = "Cabezon Larvae: Positive vs. Total Tows",
    subtitle = "Stacked bar charts show larve-present (blue) against total sampling effort (gray)",
    x        = NULL,
    y        = "Number of Tows",
    fill     = NULL
  )

#---- 8.d. Boxplots: log-larvae by year, month, season (positive tows only)

# Data binning
ss_box_data <- bind_rows(
  cabezon_shelf_spawn %>% mutate(group_var = "Year",   level = as.character(year)),
  cabezon_shelf_spawn %>% mutate(group_var = "Month",  level = as.character(month))
) %>%
  filter(larvae_100m3 > 0) %>%
  mutate(
    log_larvae = log(larvae_100m3),
    group_var  = factor(group_var, levels = c("Year", "Month", "Season")),
    level = case_when(
      group_var == "Month"  ~ factor(level, levels = as.character(1:12)),
      group_var == "Season" ~ factor(level, levels = c("winter", "spring", "summer", "fall")),
      TRUE                  ~ factor(level)
    )
  )

# Shared elements
shared_box_theme <- theme_minimal(base_size = 12) +
  theme(
    strip.text         = element_text(face = "bold", size = 13),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold", size = 15),
    plot.subtitle      = element_text(color = "grey40", size = 10)
  )

# Years plot
ggplot(ss_box_data %>% filter(group_var == "Year"), 
       aes(x = level, y = log_larvae)) +
  geom_boxplot(fill = "steelblue", color = "grey30", alpha = 0.7, outlier.size = 1) +
  shared_box_theme +
  labs(
    title    = "Cabezon Larvae: Average Count by Year (Positive Tows Only)",
    subtitle = "Log-transformed larvae per 100m³; zero-catch tows excluded",
    x        = NULL,
    y        = "log(Larvae per 100m³)"
  )

# Months + Seasons plot
ggplot(ss_box_data %>% filter(group_var %in% c("Month", "Season")),
       aes(x = level, y = log_larvae)) +
  geom_boxplot(fill = "steelblue", color = "grey30", alpha = 0.7, outlier.size = 1) +
  facet_wrap(~ group_var, scales = "free_x", nrow = 1) +
  shared_box_theme +
  labs(
    title    = "Cabezon Larvae: Average Count by Month and Season (Positive Tows Only)",
    subtitle = "Log-transformed larvae per 100m³; zero-catch tows excluded",
    x = NULL, y = "log(Larvae per 100m³)"
  )


#---- 8.e. Heatmaps: mean log-larvae by station x year, month, season

heatmap_theme <- theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

heatmap_fill <- scale_fill_viridis_c(option = "mako", name = "log(Larvae\nper 100m³)", na.value = "white")

# Year plot
cabezon_shelf_spawn %>%
  filter(larvae_100m3 > 0) %>%
  group_by(line.station, year) %>%
  summarise(mean_log_larvae = mean(log(larvae_100m3)), .groups = "drop") %>%
  ggplot(aes(x = year, y = line.station, fill = mean_log_larvae)) +
  geom_tile() + heatmap_fill + scale_y_discrete(limits = rev) + heatmap_theme +
  labs(title    = "Cabezon Larval Observations by Station and Year",
       subtitle = "Positive tows only; white = no positive tows recorded",
       x = NULL, y = "South \u2190 Station \u2192 North")

# Month plot
cabezon_shelf_spawn %>%
  filter(larvae_100m3 > 0) %>%
  group_by(line.station, month) %>%
  summarise(mean_log_larvae = mean(log(larvae_100m3)), .groups = "drop") %>%
  ggplot(aes(x = month, y = line.station, fill = mean_log_larvae)) +
  geom_tile() + heatmap_fill +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_discrete(limits = rev) + heatmap_theme +
  labs(title    = "Cabezon Larval Observations by Station and Month",
       subtitle = "Positive tows only; white = no positive tows recorded",
       x = NULL, y = "South \u2190 Station \u2192 North")

# Season plot
cabezon_shelf_spawn %>%
  filter(larvae_100m3 > 0) %>%
  group_by(line.station, season) %>%
  summarise(mean_log_larvae = mean(log(larvae_100m3)), .groups = "drop") %>%
  ggplot(aes(x = season, y = line.station, fill = mean_log_larvae)) +
  geom_tile() + heatmap_fill + scale_y_discrete(limits = rev) + heatmap_theme +
  labs(title    = "Cabezon Larval Observations by Station and Season",
       subtitle = "Positive tows only; white = no positive tows recorded",
       x = NULL, y = "South \u2190 Station \u2192 North")

ss_pos_vals <- cabezon_shelf_spawn$larvae_100m3[cabezon_shelf_spawn$larvae_100m3 > 0]
hist(log(ss_pos_vals), breaks = 50, main = "log(positive larvae_100m3)")
hist(log(cabezon_shelf_spawn$larvae_100m3[cabezon_shelf_spawn$larvae_100m3 > 0]))

cabezon_shelf_spawn %>%
  filter(larvae_100m3 > 0) %>%
  group_by(year) %>%
  summarize(m = mean(larvae_100m3), v = var(larvae_100m3), n = n()) %>%
  filter(n > 5) %>%
  ggplot(aes(log(m), log(v))) + geom_point() + geom_smooth(method = "lm") +
  geom_abline(slope = 2, intercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(title = "Mean-variance relationship (dashed: red=gamma slope 2, blue=Poisson slope 1)")

# Back out the raw count from density and volume
cabezon_shelf_spawn %>%
  filter(larvae_100m3 > 0) %>%
  mutate(
    raw_count_est = larvae_100m3 * volume_sampled / 100,
    log_density = log(larvae_100m3)
  ) %>%
  ggplot(aes(log_density, fill = factor(round(raw_count_est)))) +
  geom_histogram(bins = 50, position = "stack") +
  scale_fill_viridis_d() +
  labs(fill = "Estimated raw count",
       title = "Does multi-modality decompose into integer count strata?")


# -==============================================================================-
# ==== 3.  sdmTMB MODEL FITTING (CANDIDATE INDEX MODELS) ====
# -==============================================================================-

# Fits null models (3 families), the delta-lognormal candidate set
# (dln0-dln7), and delta-gamma / Tweedie variants used by the Appendix.


set.seed(my.seed)

#---- 3.a. Build mesh
if (!all(c("X", "Y") %in% names(cabezon_shelf_spawn))) {
  cabezon_shelf_spawn <- sdmTMB::add_utm_columns(cabezon_shelf_spawn, ll_names = c("longitude", "latitude"))
}

ss_mesh <- fm_mesh_2d(
  loc = cabezon_shelf_spawn[,c("X","Y")],
  cutoff = 20,
  max.edge = c(75, 150),
  offset = c(45, 120)
)
ss_cab_mesh <- make_mesh(data = cabezon_shelf_spawn, c("X", "Y"), mesh = ss_mesh)
plot(ss_cab_mesh); title("SPDE Mesh for Cabezon CalCOFI analysis (Constrained)")

#---- 3.b. Fit candidate models

# Add year factor variable
cabezon_shelf_spawn$fyear <- as.factor(cabezon_shelf_spawn$year)

# NULL Delta-Gamma model

ss_cab_null <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = delta_gamma(),
  spatial = "off",
  spatiotemporal = "off",
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_null)

# NULL Delta-Lognormal model

ss_cab_null_dln <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = delta_lognormal(),
  spatial = "off",
  spatiotemporal = "off",
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_null_dln)
AIC(ss_cab_null_dln)

# NULL Tweedie model

ss_cab_null_tw <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = tweedie(),
  spatial = "off",
  spatiotemporal = "off",
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_null_tw)

null_AIC_table <- data.frame(
  fit_name = c(
    "ss_cab_null",
    "ss_cab_null_dln",
    "ss_cab_null_tw"
  ),
  model = c(
    "Delta-Gamma NULL (Constrained)",
    "Delta-Lognormal NULL (Constrained)",
    "Tweedie NULL (Constrained)"
  ),
  AIC = c(
    AIC(ss_cab_null),
    AIC(ss_cab_null_dln),
    AIC(ss_cab_null_tw)
  ),
  delta_AIC = NA
)

null_AIC_table$delta_AIC <- null_AIC_table$AIC - min(null_AIC_table$AIC)

null_AIC_table <- null_AIC_table[order(null_AIC_table$AIC, decreasing = F), ]

null_AIC_table

# Delta-Lognormal: Global Intercept, spatial on, spatiotemporal off
ss_cab_fit_dln0 <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = delta_lognormal(),
  spatial = list("on", "on"),
  spatiotemporal = list("off", "off"),
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_fit_dln0)

AIC(ss_cab_fit_dln0)

# Delta-Lognormal: Global Intercept, spatial on, spatiotemporal delta AR1, spatiotemporal lognormal AR1
ss_cab_fit_dln1 <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = delta_lognormal(),
  spatial = list("on", "on"),
  spatiotemporal = list("ar1", "ar1"),
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_fit_dln1)

AIC(ss_cab_fit_dln1)

# Delta-Lognormal: Global Intercept, spatial on, spatiotemporal delta AR1, spatiotemporal lognormal iid
ss_cab_fit_dln2 <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = delta_lognormal(),
  spatial = list("on", "on"),
  spatiotemporal = list("ar1", "iid"),
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_fit_dln2)

AIC(ss_cab_fit_dln2)

# Delta-Lognormal: Global Intercept, spatial on, spatiotemporal delta iid, spatiotemporal lognormal off
ss_cab_fit_dln3 <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = delta_lognormal(),
  spatial = list("on", "on"),
  spatiotemporal = list("iid", "off"),
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_fit_dln3)

AIC(ss_cab_fit_dln3)

# Delta-Lognormal: fixed year effect, spatial on, spatiotemporal off
ss_cab_fit_dln4 <- sdmTMB(
  larvae_100m3 ~ 0 + fyear,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = delta_lognormal(),
  spatial = list("on", "on"),
  spatiotemporal = list("off", "off"),
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_fit_dln4)

AIC(ss_cab_fit_dln4)

# Delta-Lognormal: fixed year effect, spatial on, spatiotemporal iid
ss_cab_fit_dln5 <- sdmTMB(
  larvae_100m3 ~ 0 + fyear,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = delta_lognormal(),
  spatial = list("on", "on"),
  spatiotemporal = list("iid", "iid"),
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_fit_dln5)

AIC(ss_cab_fit_dln5)

# Delta-Lognormal: fixed year effect, spatial on, iid delta spatiotemporal, lognormal spatiotemporal off
ss_cab_fit_dln6 <- sdmTMB(
  larvae_100m3 ~ 0 + fyear,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = delta_lognormal(),
  spatial = list("on", "on"),
  spatiotemporal = list("iid", "off"),
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_fit_dln6)

AIC(ss_cab_fit_dln6)

# Delta-Lognormal: Global Intercept, spatial on, spatiotemporal delta AR1, spatiotemporal lognormal off
ss_cab_fit_dln7 <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = delta_lognormal(),
  spatial = list("on", "on"),
  spatiotemporal = list("ar1", "off"),
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_fit_dln7)

AIC(ss_cab_fit_dln7)


#---- 3.c. Random effect parameter estimates for best models

tidy(ss_cab_fit_dln1, effects = "ran_pars", model = 1); tidy(ss_cab_fit_dln1, effects = "ran_pars", model = 2)
tidy(ss_cab_fit_dln2, effects = "ran_pars", model = 1); tidy(ss_cab_fit_dln2, effects = "ran_pars", model = 2)
tidy(ss_cab_fit_dln3, effects = "ran_pars", model = 1); tidy(ss_cab_fit_dln3, effects = "ran_pars", model = 2)
tidy(ss_cab_fit_dln4, effects = "ran_pars")
tidy(ss_cab_fit_dln5, effects = "ran_pars")

#---- 3.d. Construct model comparison table

AIC_table <- data.frame(
  fit_name = c(
    "ss_cab_fit_dln0",
    "ss_cab_fit_dln1",
    "ss_cab_fit_dln2",
    "ss_cab_fit_dln3",
    "ss_cab_fit_dln4",
    "ss_cab_fit_dln5",
    "ss_cab_fit_dln6",
    "ss_cab_fit_dln7"
  ),
  model = c(
    "Delta-Lognormal Global Intercept, spatial on, spatiotemporal off",
    "Delta-Lognormal Global Intercept + spatial + spatiotemporal delta AR1 + spatiotemporal lognormal AR1",
    "Delta-Lognormal Global Intercept + spatial + spatiotemporal delta AR1 + spatiotemporal lognormal iid",
    "Delta-Lognormal Global Intercept + spatial + spatiotemporal delta iid",
    "Delta-Lognormal fixed year effect + spatial + spatiotemporal off",
    "Delta-Lognormal fixed year effect + spatial + spatiotemporal delta iid + spatiotemporal lognormal iid",
    "Delta-Lognormal fixed year effect + spatial + spatiotemporal delta iid + spatiotemporal lognormal off",
    "Delta-Lognormal Global Intercept + spatial + spatiotemporal delta AR1 + spatiotemporal lognormal off"
  ),
  AIC = c(
    AIC(ss_cab_fit_dln0),
    AIC(ss_cab_fit_dln1),
    AIC(ss_cab_fit_dln2),
    AIC(ss_cab_fit_dln3),
    AIC(ss_cab_fit_dln4),
    AIC(ss_cab_fit_dln5),
    AIC(ss_cab_fit_dln6),
    AIC(ss_cab_fit_dln7)
  ),
  delta_AIC = NA
)

AIC_table$delta_AIC <- AIC_table$AIC - min(AIC_table$AIC)

AIC_table <- AIC_table[order(AIC_table$AIC, decreasing = F), ]

AIC_table

ss_cab_fit_dg1 <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = delta_gamma(),
  spatial = list("on", "on"),
  spatiotemporal = list("ar1", "off"),
  offset = NULL,
  extra_time = 1982
)

sanity(ss_cab_fit_dg1)
AIC(ss_cab_fit_dg1)

tw_iid_fit <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = tweedie(),
  spatial = "on",
  spatiotemporal = "iid",
  offset = NULL,
  extra_time = 1982
)
sanity(tw_iid_fit)

tw_stoff_fit <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = tweedie(),
  spatial = "on",
  spatiotemporal = "off",
  offset = NULL,
  extra_time = 1982
)
sanity(tw_stoff_fit)

tw_ar1_fit <- sdmTMB(
  larvae_100m3 ~ 1,
  data = cabezon_shelf_spawn,
  mesh = ss_cab_mesh,
  time = "year",
  family = tweedie(),
  spatial = "on",
  spatiotemporal = "ar1",
  offset = NULL,
  extra_time = 1982
)
sanity(tw_ar1_fit)

AIC(tw_iid_fit, tw_stoff_fit, tw_ar1_fit)




# -==============================================================================-
# ==== 4.  PREFERRED-MODEL DIAGNOSTICS & ABUNDANCE INDEX ====
# -==============================================================================-

# Preferred model ss_cab_fit_dln3 (DLn-IID). Produces the regular-grid
# index ss_dln3_cab_index (carried downstream) and the observed-station-grid
# variant ss_dln3_cab_index_obs (used only in the Appendix).


summary(ss_cab_fit_dln3)

#---- 4.1.a. Calibration summaries

set.seed(my.seed)
ss_dln3_pred_obs <- predict(ss_cab_fit_dln3, type = "response")

# Overall mean calibration
ss_dln3_pred_obs %>%
  summarise(
    obs_mean = mean(larvae_100m3, na.rm = TRUE),
    pred_mean = mean(est, na.rm = TRUE),
    obs_occ = mean(larvae_100m3 > 0, na.rm = TRUE),
    pred_occ = mean(est1, na.rm = TRUE),
    obs_pos_mean = mean(larvae_100m3[larvae_100m3 > 0], na.rm = TRUE),
    pred_pos_mean = mean(est2, na.rm = TRUE)
  )

# Binomial calibration; i.e. predicted presence probability vs. observed presence rate
ss_dln3_pred_obs$presence <- ss_dln3_pred_obs$larvae_100m3 > 0

ss_dln3_calibration_summary <- ss_dln3_pred_obs %>%
  mutate(prob_bin = cut(est1, breaks = seq(0, 1, by = 0.1))) %>%
  group_by(prob_bin) %>%
  summarise(
    mean_pred = mean(est1),
    observed = mean(presence),
    n = n()
  )

plot(ss_dln3_calibration_summary$mean_pred, ss_dln3_calibration_summary$observed,
     xlab = "Predicted presence probability",
     ylab = "Observed presence frequency")
abline(0, 1, lty = 2)

#---- 4.1.b. Predicted vs. observed plots

# Raw scale comparison
plot(
  x = ss_dln3_pred_obs$est,
  y = cabezon_shelf_spawn$larvae_100m3,
  xlab = "Predicted",
  ylab = "Observed"
)
abline(0, 1, lty = 2)

# Log scale comparison
plot(log1p(ss_dln3_pred_obs$est), log1p(cabezon_shelf_spawn$larvae_100m3),
     xlab = "log(Predicted + 1)", ylab = "log(Observed + 1)")
abline(0, 1, lty = 2)

# Lognormal (positive density) component only
pos_idx <- cabezon_shelf_spawn$larvae_100m3 > 0
plot(log(ss_dln3_pred_obs$est2[pos_idx]), log(cabezon_shelf_spawn$larvae_100m3[pos_idx]),
     xlab = "log(Predicted positive abundance)",
     ylab = "log(Observed positive abundance)")
abline(0, 1, lty = 2)

#---- 4.1.c. Randomized quantile residuals

set.seed(my.seed)
ss_dln3_rq_res <- residuals(ss_cab_fit_dln3, type = "mle-mvn")

# QQ plot
qqnorm(ss_dln3_rq_res)
qqline(ss_dln3_rq_res)

# Spatial residual map by year
ss_dln3_pred_obs$rq_resid <- ss_dln3_rq_res

ggplot(ss_dln3_pred_obs, aes(x = longitude, y = latitude, color = rq_resid)) +
  geom_point() +
  scale_color_gradient2() +
  facet_wrap(~ year)

#---- 4.1.d. Temporal calibration

# Are predictions tracking temporal trends?
ss_dln3_pred_obs %>%
  group_by(year) %>%
  summarize(
    obs = mean(larvae_100m3),
    pred = mean(est)
  ) %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = obs), color = "black") +
  geom_line(aes(y = pred), color = "blue", linetype = "dashed")

#---- 4.1.e. Abundance index construction
# Note! The index uncertainty here was far too wide using a regular spatial grid, probably because interpolation over years with no CalCOFI coverage was causing inflation. I kept that plot in for reference, but 4.1.f. (observed station grid) is probably preferred for diagnosis.

ss_dln3_pred_grid <- expand.grid(
  X = seq(min(cabezon_shelf_spawn$X), max(cabezon_shelf_spawn$X), by = 5),
  Y = seq(min(cabezon_shelf_spawn$Y), max(cabezon_shelf_spawn$Y), by = 5)
) %>%
  tidyr::crossing(year = as.integer(sort(unique(cabezon_shelf_spawn$year)))) %>%
  mutate(fyear = as.factor(year))

ss_dln3_cab_map <- predict(ss_cab_fit_dln3, newdata = ss_dln3_pred_grid, return_tmb_object = TRUE)
ss_dln3_cab_index <- get_index(ss_dln3_cab_map, area = 1, bias_correct = TRUE)

ggplot(ss_dln3_cab_index, aes(x = year, y = est)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  labs(y = "Estimated abundance index", x = "Year",
       title = "Abundance Index: Delta-Lognormal on Global Intercept with IID (encounter) Spatiotemporal Field")

#---- 4.1.f. Abundance index: constrained to observed station grid (probably preferred)

ss_dln3_pred_grid_obs <- cabezon_shelf_spawn %>%
  dplyr::select(X, Y, year) %>%
  distinct()

ss_dln3_cab_map_obs <- predict(ss_cab_fit_dln3,
                               newdata = ss_dln3_pred_grid_obs,
                               return_tmb_object = TRUE)

ss_dln3_cab_index_obs <- get_index(ss_dln3_cab_map_obs,
                                   area = 1,
                                   bias_correct = TRUE)

ggplot(ss_dln3_cab_index_obs, aes(x = year, y = est)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  labs(y = "Estimated abundance index", x = "Year",
       title = "Abundance Index: Observed Station Grid (Delta-Lognormal, IID Encounter)")


# -==============================================================================-
# ==== 5.  STOCK-ASSESSMENT COMPARISONS & CROSS-CORRELATIONS ====
# -==============================================================================-

# Canonical basis: regular-grid index ss_dln3_cab_index. Builds
# ss_dln3_combined / ss_dln3_combined_obs used by Sections 6-8.

#---- 5.1. ss_cab_fit_dln3 STAR comparisons and CCF
# Create one data frame to house all STAR + larval index data
ss_dln3_combined <- ss_dln3_cab_index %>%
  dplyr::select(year, est, lwr, upr) %>%
  rename(larvae_index = est) %>%
  left_join(
    STAR_SSB %>% rename(ssb = value, ssb_lo = lo, ssb_hi = hi),
    by = "year"
  ) %>%
  left_join(
    STAR_recdevs %>% rename(rec_dev = value, rec_dev_lo = lo, rec_dev_hi = hi),
    by = "year"
  ) %>%
  left_join(
    STAR_age0 %>% rename(age0 = value, age0_lo = lo, age0_hi = hi),
    by = "year"
  )
# Filter to observed years
ss_dln3_combined_obs <- ss_dln3_combined %>%
  filter(!is.na(larvae_index))
# CCF: Restricted larval index vs. age-0 recruits
ccf_boot(
  x = ss_dln3_combined_obs$larvae_index,
  y = ss_dln3_combined_obs$age0,
  lag.max = 5,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. Age-0 Recruits (Delta-Lognormal)")
# CCF: Restricted larval index vs. recruitment deviations
ccf_boot(
  x = ss_dln3_combined_obs$larvae_index,
  y = ss_dln3_combined_obs$rec_dev,
  lag.max = 5,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. Recruitment Deviations (Delta-Lognormal)")
# Nothing really :( when we try against SSB, we get some strong signal, but this is probably due to a shared declining trend. We'll try to parse that.
# CCF: Restricted larval index vs. SSB
ccf_boot(
  x = ss_dln3_combined_obs$larvae_index,
  y = ss_dln3_combined_obs$ssb,
  lag.max = 5,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. SSB (Delta-Lognormal)")
# First-differenced
ccf_boot(
  x = diff(ss_dln3_combined_obs$larvae_index),
  y = diff(ss_dln3_combined_obs$ssb),
  lag.max = 5,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Δ Larval Index vs. Δ SSB (Delta-Lognormal)")
# And only differenced SSB on raw larvae
ccf_boot(
  x = ss_dln3_combined_obs$larvae_index,
  y = diff(ss_dln3_combined_obs$ssb),
  lag.max = 5,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Larval Index vs. Δ SSB (Delta-Lognormal)")
# Join differenced SSB back to index
ssb_diff_df_dln3 <- ss_dln3_combined_obs %>%
  arrange(year) %>%
  mutate(ssb_diff = c(NA, diff(ssb))) %>%
  filter(!is.na(ssb_diff))
lm_lag0_dln3 <- lm(ssb_diff ~ larvae_index, data = ssb_diff_df_dln3)
summary(lm_lag0_dln3)
# Visual
ssb_diff_df_dln3 %>%
  ggplot(aes(x = larvae_index, y = ssb_diff, label = year)) +
  geom_point() +
  geom_text(nudge_y = 50, size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Larval Index", y = "ΔSSB (mt)",
       title = "Contemporaneous larval index vs. change in SSB (Delta-Lognormal)") +
  theme_minimal()
# Assuming ecologically that SSB would have the strongest relationship to larval abundance, we'll plot those two time series visually for a qualitative look
scale_factor_ssb_dln3 <- max(ss_dln3_combined_obs$larvae_index, na.rm = TRUE) /
  max(ss_dln3_combined_obs$ssb, na.rm = TRUE)
ggplot(ss_dln3_combined_obs %>% filter(year >= 1984, year <= 2015),
       aes(x = year)) +
  # SSB scaled UP to primary axis units
  geom_ribbon(aes(ymin = ssb_lo * scale_factor_ssb_dln3,
                  ymax = ssb_hi * scale_factor_ssb_dln3),
              fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = ssb * scale_factor_ssb_dln3), color = "steelblue", linewidth = 0.8) +
  # Larvae index on primary axis
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "coral", alpha = 0.3) +
  geom_line(aes(y = larvae_index), color = "coral", linewidth = 0.8) +
  scale_y_continuous(
    name = "Estimated Larval Abundance Index; Delta-Lognormal IID (encounter)",
    sec.axis = sec_axis(~ . / scale_factor_ssb_dln3,
                        name = "SSB (mt)")
  ) +
  scale_x_continuous(breaks = seq(1984, 2015, by = 4)) +
  labs(title = "SoCal Cabezon: Larval Index vs. SSB (1984–2015), Delta-Lognormal",
       x = "Year") +
  theme_bw()

# Let's look closer at SSB to see if different span specifications change the relationship to our index. Specifically, what proportion of variance in our SSB is attributable to a smooth, long-term trend across three span sensitivities?
ssb_smooth_5 <- loess(ssb ~ year, data = ss_dln3_combined, span = 0.5)$fitted
var(ssb_smooth) / var(ss_dln3_combined$ssb)

idx_smooth_5 <- loess(larvae_index ~ year, data = ss_dln3_combined, span = 0.5)$fitted
var(idx_smooth_5) / var(ss_dln3_combined$larvae_index)

ssb_smooth_3 <- loess(ssb ~ year, data = ss_dln3_combined, span = 0.3)$fitted
var(ssb_smooth_3) / var(ss_dln3_combined$ssb)

idx_smooth_3 <- loess(larvae_index ~ year, data = ss_dln3_combined, span = 0.3)$fitted
var(idx_smooth_3) / var(ss_dln3_combined$larvae_index)

ssb_smooth_7 <- loess(ssb ~ year, data = ss_dln3_combined, span = 0.7)$fitted
var(ssb_smooth_7) / var(ss_dln3_combined$ssb)

idx_smooth_7 <- loess(larvae_index ~ year, data = ss_dln3_combined, span = 0.7)$fitted
var(idx_smooth_7) / var(ss_dln3_combined$larvae_index)

# ---- Larval Index vs. Recruitment Deviations
# Rec devs are log-scale and centered at 0, so a multiplicative scale_factor
# doesn't behave nicely. Instead, shift + scale rec devs into the index range.

idx_range_dln3   <- range(ss_dln3_combined_obs$larvae_index, na.rm = TRUE)
recdev_range_dln3 <- range(c(ss_dln3_combined_obs$rec_dev_lo,
                             ss_dln3_combined_obs$rec_dev_hi), na.rm = TRUE)

# Linear map from rec dev space -> index space
slope_rd_dln3 <- diff(idx_range_dln3) / diff(recdev_range_dln3)
intercept_rd_dln3 <- idx_range_dln3[1] - slope_rd_dln3 * recdev_range_dln3[1]

ggplot(ss_dln3_combined_obs %>% filter(year >= 1984, year <= 2015),
       aes(x = year)) +
  # Rec devs scaled to primary axis units
  geom_ribbon(aes(ymin = intercept_rd_dln3 + slope_rd_dln3 * rec_dev_lo,
                  ymax = intercept_rd_dln3 + slope_rd_dln3 * rec_dev_hi),
              fill = "darkgreen", alpha = 0.25) +
  geom_line(aes(y = intercept_rd_dln3 + slope_rd_dln3 * rec_dev),
            color = "darkgreen", linewidth = 0.8) +
  geom_hline(yintercept = intercept_rd_dln3 + slope_rd_dln3 * 0,
             linetype = "dashed", color = "darkgreen", alpha = 0.6) +
  # Larvae index on primary axis
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "coral", alpha = 0.3) +
  geom_line(aes(y = larvae_index), color = "coral", linewidth = 0.8) +
  scale_y_continuous(
    name = "Estimated Larval Abundance Index; Delta-Lognormal IID (encounter)",
    sec.axis = sec_axis(~ (. - intercept_rd_dln3) / slope_rd_dln3,
                        name = "Recruitment Deviation (log-scale)")
  ) +
  scale_x_continuous(breaks = seq(1984, 2015, by = 4)) +
  labs(title = "SoCal Cabezon: Larval Index vs. Recruitment Deviations (1984–2015), Delta-Lognormal",
       x = "Year") +
  theme_bw()


# ---- Larval Index vs. Age-0 Recruits
# Age-0 is strictly positive like SSB, so the multiplicative scale_factor approach works.

scale_factor_age0_dln3 <- max(ss_dln3_combined_obs$larvae_index, na.rm = TRUE) /
  max(ss_dln3_combined_obs$age0, na.rm = TRUE)

ggplot(ss_dln3_combined_obs %>% filter(year >= 1984, year <= 2015),
       aes(x = year)) +
  # Age-0 scaled UP to primary axis units
  geom_ribbon(aes(ymin = age0_lo * scale_factor_age0_dln3,
                  ymax = age0_hi * scale_factor_age0_dln3),
              fill = "purple", alpha = 0.25) +
  geom_line(aes(y = age0 * scale_factor_age0_dln3),
            color = "purple", linewidth = 0.8) +
  # Larvae index on primary axis
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "coral", alpha = 0.3) +
  geom_line(aes(y = larvae_index), color = "coral", linewidth = 0.8) +
  scale_y_continuous(
    name = "Estimated Larval Abundance Index; Delta-Lognormal IID (encounter)",
    sec.axis = sec_axis(~ . / scale_factor_age0_dln3,
                        name = "Age-0 Recruits")
  ) +
  scale_x_continuous(breaks = seq(1984, 2015, by = 4)) +
  labs(title = "SoCal Cabezon: Larval Index vs. Age-0 Recruits (1984–2015), Delta-Lognormal",
       x = "Year") +
  theme_bw()


# -==============================================================================-
# ==== 6.  REGRESSION ANALYSES (INDEX vs. RECRUITMENT) ====
# -==============================================================================-

# RQ residuals for delta-lognormal IID (encounter) spatiotemporal model
ss_dln3_pred_obs %>%
  group_by(year) %>%
  mutate(med_resid = median(rq_resid)) %>%
  ggplot(aes(x = factor(year), y = rq_resid)) +
  geom_boxplot(outlier.size = 0.6, fill = "steelblue", alpha = 0.4) +
  geom_smooth(
    aes(x = as.numeric(factor(year)), y = rq_resid),
    method  = "loess", span = 0.4,
    se      = TRUE, color = "firebrick", linewidth = 0.8
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x     = "Year",
    y     = "RQ Residual",
    title = "Annual RQ Residual Distributions — ss_cab_fit_dln3"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 7))

ss_dln3_pred_obs %>%
  group_by(year) %>%
  summarise(med_resid = median(rq_resid), .groups = "drop") %>%
  arrange(year) %>%
  mutate(cusum = cumsum(med_resid)) %>%
  ggplot(aes(x = year, y = cusum)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x     = "Year",
    y     = "Cumulative Sum of Median RQ Residual",
    title = "CUSUM — ss_cab_fit_dln3"
  ) +
  theme_minimal()

# Is the underprediction concentrated in one delta component?
# Separate RQ residuals by model component
set.seed(my.seed)
ss_dln3_pred_obs <- ss_dln3_pred_obs %>%
  mutate(
    rq_resid1 = residuals(ss_cab_fit_dln3, model = 1, type = "mle-mvn"),
    rq_resid2 = residuals(ss_cab_fit_dln3, model = 2, type = "mle-mvn")
  )

# CUSUM for each component separately
ss_dln3_pred_obs %>%
  group_by(year) %>%
  summarise(
    med_resid1 = median(rq_resid1, na.rm = TRUE),
    med_resid2 = median(rq_resid2, na.rm = TRUE),  # positive tows only
    .groups = "drop"
  ) %>%
  arrange(year) %>%
  mutate(
    cusum1 = cumsum(replace_na(med_resid1, 0)),
    cusum2 = cumsum(replace_na(med_resid2, 0))
  ) %>%
  pivot_longer(cols = c(cusum1, cusum2),
               names_to = "component",
               values_to = "cusum") %>%
  ggplot(aes(x = year, y = cusum, color = component)) +
  geom_line(linewidth = 0.9) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("cusum1" = "steelblue", "cusum2" = "firebrick"),
                     labels = c("Delta (binomial)", "Lognormal")) +
  labs(x = "Year", y = "CUSUM",
       title = "CUSUM by Delta Model Component — ss_cab_fit_dln3") +
  theme_minimal()

# Is the underprediction spatially structured?
ss_dln3_pred_obs %>%
  group_by(station) %>%
  summarise(mean_resid = mean(rq_resid), .groups = "drop") %>%
  ggplot(aes(x = station, y = mean_resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = TRUE) +
  labs(x = "Station", y = "Mean RQ Residual",
       title = "Mean Residual by Station — ss_cab_fit_dln3") +
  theme_minimal()


#---- Simple lagged regression at -1 to check

ss_dln3_combined_lagged <- ss_dln3_combined_obs %>%
  arrange(year) %>%
  left_join(
    ss_dln3_combined_obs %>%
      dplyr::select(year, age0_lead1 = age0, rec_dev_lead1 = rec_dev) %>%
      mutate(year = year - 1),  # shift target year back by 1
    by = "year"
  )

ss_dln3_lm_age0   <- lm(age0_lead1    ~ larvae_index, data = ss_dln3_combined_lagged)
ss_dln3_lm_recdev <- lm(rec_dev_lead1 ~ larvae_index, data = ss_dln3_combined_lagged)

summary(ss_dln3_lm_age0)
summary(ss_dln3_lm_recdev)

ss_dln3_combined_lagged %>%
  filter(!is.na(age0_lead1)) %>%
  ggplot(aes(x = larvae_index, y = age0_lead1, label = year)) +
  geom_point() +
  geom_text(nudge_y = 8, size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Larval Index (year t)", y = "Age-0 Recruits (year t+1)",
       title = "Lagged Larval Index vs. Age-0 Recruits — ss_cab_fit_dln3") +
  theme_minimal()

ss_dln3_combined_lagged_no83 <- ss_dln3_combined_lagged %>% filter(year != 1983)

ss_dln3_lm_age0_no83 <- lm(age0_lead1 ~ larvae_index, data = ss_dln3_combined_lagged_no83)
summary(ss_dln3_lm_age0_no83)

ss_dln3_lm_age0_detrend   <- lm(age0_lead1    ~ larvae_index + year, data = ss_dln3_combined_lagged)
ss_dln3_lm_recdev_detrend <- lm(rec_dev_lead1 ~ larvae_index + year, data = ss_dln3_combined_lagged)

summary(ss_dln3_lm_age0_detrend)
summary(ss_dln3_lm_recdev_detrend)

plot(ss_dln3_lm_age0, which = 4)  # Cook's distance plot

# Remove the three most anomalous years
ss_dln3_combined_lagged_trimmed <- ss_dln3_combined_lagged %>%
  filter(!year %in% c(1983, 1986, 1998))

ss_dln3_combined_lagged_trimmed %>%
  filter(!is.na(age0_lead1)) %>%
  ggplot(aes(x = larvae_index, y = age0_lead1, label = year)) +
  geom_point() +
  geom_text(nudge_y = 8, size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Larval Index (year t)", y = "Age-0 Recruits (year t+1)",
       title = "Lagged Larval Index vs. Age-0 Recruits (trimmed) — ss_cab_fit_dln3") +
  theme_minimal()

ss_dln3_lm_age0_trimmed <- lm(age0_lead1 ~ larvae_index, data = ss_dln3_combined_lagged_trimmed)
summary(ss_dln3_lm_age0_trimmed)

ss_dln3_lm_age0_trimmed_detrend <- lm(age0_lead1 ~ larvae_index + year, data = ss_dln3_combined_lagged_trimmed)
summary(ss_dln3_lm_age0_trimmed_detrend)



# -==============================================================================-
# ==== 7.  PUBLICATION FIGURES & TABLES ====
# -==============================================================================-

## ---- Predicted-abundance maps (preferred model) ------------------------

library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

utm_crs <- 32611

ss_map_pred <- predict(
  ss_cab_fit_dln3,
  newdata = ss_dln3_pred_grid,
  type = "response"
)

pred_plot <- ss_map_pred |>
  mutate(
    expected = est1 * est2,
    X_m = X * 1000,
    Y_m = Y * 1000
  )

ca_coast <- ne_states(
  country = "United States of America",
  returnclass = "sf"
) |>
  dplyr::filter(name == "California") |>
  st_transform(utm_crs)

calcofi_stations <- cabezon_shelf_spawn |>
  mutate(
    X_m = X * 1000,
    Y_m = Y * 1000
  ) |>
  distinct(line, station, X_m, Y_m) |>
  st_as_sf(coords = c("X_m", "Y_m"), crs = utm_crs)

calcofi_lines <- calcofi_stations |>
  arrange(line, station) |>
  group_by(line) |>
  summarize(do_union = FALSE, .groups = "drop") |>
  st_cast("LINESTRING")

pred_plot <- ss_map_pred |>
  mutate(X_m = X * 1000, Y_m = Y * 1000) |>
  summarise(
    est = mean(est),
    est1 = mean(est1),
    est2 = mean(est2),
    .by = c(X_m, Y_m))


# Plot showing larval observations by station over time
obs_plot <- cabezon_shelf_spawn |>
  mutate(
    X_m = X * 1000,
    Y_m = Y * 1000,
    larvae_positive = larvae_100m3 > 0
  )

x_range <- range(pred_plot$X_m, na.rm = TRUE)
y_range <- range(pred_plot$Y_m, na.rm = TRUE)

x_buffer <- diff(x_range) * 0.05
y_buffer <- diff(y_range) * 0.05

xlim_padded <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
ylim_padded <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

ggplot() +
  geom_sf(
    data = ca_coast,
    inherit.aes = FALSE,
    fill = "grey80",
    color = "grey40"
  ) +
  geom_sf(
    data = calcofi_lines,
    inherit.aes = FALSE,
    color = "grey95",
    linewidth = 0.6
  ) +
  geom_point(
    data = obs_plot |> filter(!larvae_positive),
    aes(x = X_m, y = Y_m),
    size = 0.4,
    alpha = 0.25
  ) +
  geom_point(
    data = obs_plot |> filter(larvae_positive),
    aes(x = X_m, y = Y_m, size = larvae_100m3),
    alpha = 0.7
  ) +
  coord_sf(
    crs = st_crs(utm_crs),
    xlim = xlim_padded,
    ylim = ylim_padded,
    expand = FALSE
  ) +
  facet_wrap(~ year) +
  scale_size_continuous(
    name = "Observed larvae / 100 m³",
    range = c(1, 7)
  ) +
  labs(
    title = "Observed cabezon larval abundance by CalCOFI station",
    x = NULL,
    y = NULL
  ) +
  theme_minimal()

# Plot showing only delta (distribution)
ggplot() +
  geom_raster(
    data = pred_plot,
    aes(x = X_m, y = Y_m, fill = log1p(est))
  ) +
  geom_sf(
    data = ca_coast,
    inherit.aes = FALSE,
    fill = "beige",
    color = "black",
    linewidth = 0.2
  ) +
  #  geom_sf(
  #    data = calcofi_lines,
  #    inherit.aes = FALSE,
  #    color = "white",
  #    linewidth = 1.2,
  #    alpha = 0.35
  #  ) +
  #  geom_point(
  #    data = obs_aggregate %>% filter(ever_positive),
  #    aes(x = X_m, y = Y_m),
  #    inherit.aes = FALSE,
  #    size = 2.5,
  #    color = "goldenrod",
  #    alpha = 1
  #  ) +
  coord_sf(
    crs = st_crs(utm_crs),
    xlim = range(pred_plot$X_m, na.rm = TRUE),
    ylim = range(pred_plot$Y_m, na.rm = TRUE),
    expand = FALSE
  ) +
  scale_fill_viridis_c(
    option = "mako",
    name = "Predicted larval abundance"
  ) +
  labs(
    title = "Persistent spatial structure in modeled larval abundance",
    subtitle = "Log-transformed expected abundance from delta-lognormal index model",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "grey65", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12, color = "grey35"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_line(color = "grey90", linewidth = 0.2)
  )

# Plot showing hurdle (abundance); same across years
ggplot() +
  geom_raster(
    data = pred_plot,
    aes(x = X_m, y = Y_m, fill = est2)
  ) +
  geom_sf(data = ca_coast, inherit.aes = FALSE, fill = "grey80") +
  coord_sf(
    crs = st_crs(utm_crs),
    xlim = range(pred_plot$X_m, na.rm = TRUE),
    ylim = range(pred_plot$Y_m, na.rm = TRUE),
    expand = FALSE
  ) +
  scale_fill_viridis_c(option = "mako",
                       name = "Predicted larvae") +
  theme_minimal()

pred_plot_yearly <- ss_map_pred |>
  mutate(X_m = X * 1000, Y_m = Y * 1000)

ggplot() +
  geom_raster(
    data = pred_plot_yearly,
    aes(x = X_m, y = Y_m, fill = log(est))
  ) +
  geom_sf(data = ca_coast, inherit.aes = FALSE, fill = "grey80") +
  coord_sf(
    crs = st_crs(utm_crs),
    xlim = range(pred_plot$X_m, na.rm = TRUE),
    ylim = range(pred_plot$Y_m, na.rm = TRUE),
    expand = FALSE
  ) +
  facet_wrap(~ year) +
  scale_fill_viridis_c(name = "Predicted larvae ln(abun)") +
  theme_minimal()


## ---- Cross-correlation panel figure ------------------------------------

library(dplyr)
library(ggplot2)
library(funtimes)

# ── 1. EXTRACT CCF RESULTS ──────────────────────────────────────────────
# Compute LOESS residuals first
loess_df <- ss_dln3_combined_obs %>%
  arrange(year) %>%
  mutate(
    loess_index = predict(loess(larvae_index ~ year, span = 0.5)),
    loess_ssb   = predict(loess(ssb ~ year, span = 0.5)),
    resid_index = larvae_index - loess_index,
    resid_ssb   = ssb - loess_ssb
  )

ccf_age0   <- ccf_boot(ss_dln3_combined_obs$larvae_index,
                       ss_dln3_combined_obs$age0,
                       lag.max = 10, plot = "Spearman", B = 1000)
ccf_recdev <- ccf_boot(ss_dln3_combined_obs$larvae_index,
                       ss_dln3_combined_obs$rec_dev,
                       lag.max = 10, plot = "Spearman", B = 1000)
ccf_ssb    <- ccf_boot(ss_dln3_combined_obs$larvae_index,
                       ss_dln3_combined_obs$ssb,
                       lag.max = 10, plot = "Spearman", B = 1000)
ccf_dssb   <- ccf_boot(diff(ss_dln3_combined_obs$larvae_index),
                       diff(ss_dln3_combined_obs$ssb),
                       lag.max = 10, plot = "Spearman", B = 1000)


# Add: LOESS-detrended instead of first-differenced
ccf_dssl   <- ccf_boot(loess_df$resid_index,
                       loess_df$resid_ssb,
                       lag.max = 10, plot = "Spearman", B = 1000)

# tidy data frame
tidy_ccf <- function(df, panel_label) {
  df |>
    transmute(
      lag    = Lag,
      r      = r_S,
      lower  = lower_S,
      upper  = upper_S,
      sig    = r < lower | r > upper,
      panel  = panel_label
    )
}

ccf_df <- bind_rows(
  tidy_ccf(ccf_age0,   "Age-0 recruits"),
  tidy_ccf(ccf_recdev, "Recruitment deviations"),
  tidy_ccf(ccf_ssb,    "SSB (raw)"),
  tidy_ccf(ccf_dssl,   "SSB (LOESS-detrended)")
)

ccf_df$panel <- factor(
  ccf_df$panel,
  levels = c("Age-0 recruits", "Recruitment deviations",
             "SSB (raw)", "SSB (LOESS-detrended)")
)

# Figure
p_ccf <- ggplot(ccf_df, aes(x = lag, y = r)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "skyblue", alpha = 0.6) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.5,
               colour = "grey30") +
  geom_point(aes(fill = sig), shape = 21, size = 2.2,
             colour = "grey20", stroke = 0.4) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "black"),
                    labels = c("ns", "outside 95% CI"),
                    name = NULL) +
  facet_wrap(~ panel, ncol = 2) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  labs(
    x        = "Lag (years)",
    y        = "Spearman cross-correlation",
    title    = "Cross-correlations: restricted larval index vs. stock assessment outputs",
    subtitle = "Delta-lognormal IID model; blue band = bootstrap 95% CI (B = 1000)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background   = element_rect(fill = "grey95", colour = NA),
    strip.text         = element_text(face = "bold", size = 10),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(colour = "grey30", size = 9.5),
    legend.position    = "bottom",
    legend.margin      = margin(t = -5)
  )

p_ccf


ggsave("Figures/ccf_panel.png", p_ccf,
       width = 8, height = 6.5, dpi = 300)

## ---- Cross-correlation alternative (LOESS-detrended) figure ------------

library(patchwork)

# --- Compute LOESS smooths and residuals ---
loess_df <- ss_dln3_combined_obs %>%
  filter(year >= 1984, year <= 2015) %>%
  arrange(year) %>%
  mutate(
    loess_index = predict(loess(larvae_index ~ year, span = 0.5)),
    loess_ssb   = predict(loess(ssb ~ year, span = 0.5)),
    resid_index = larvae_index - loess_index,
    resid_ssb   = ssb - loess_ssb
  )

# Scale factor for dual axis
scale_factor <- max(loess_df$larvae_index, na.rm = TRUE) /
  max(loess_df$ssb, na.rm = TRUE)

# Rescale SSB residuals to index space for visual comparison
resid_scale  <- sd(loess_df$resid_index, na.rm = TRUE) /
  sd(loess_df$resid_ssb,   na.rm = TRUE)

# --- Panel A: Raw series + LOESS smooths ---
p_raw <- ggplot(loess_df, aes(x = year)) +
  # SSB ribbon + line
  geom_ribbon(aes(ymin = ssb_lo * scale_factor,
                  ymax = ssb_hi * scale_factor),
              fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = ssb * scale_factor),
            color = "steelblue", linewidth = 0.8, alpha = 0.5) +
  # SSB LOESS smooth
  geom_line(aes(y = loess_ssb * scale_factor),
            color = "steelblue", linewidth = 1.3) +
  # Larval index ribbon + line
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "coral", alpha = 0.2) +
  geom_line(aes(y = larvae_index),
            color = "coral", linewidth = 0.8, alpha = 0.5) +
  # Index LOESS smooth
  geom_line(aes(y = loess_index),
            color = "coral", linewidth = 1.3) +
  scale_y_continuous(
    name = "Larval abundance index",
    sec.axis = sec_axis(~ . / scale_factor, name = "SSB (mt)")
  ) +
  scale_x_continuous(breaks = seq(1984, 2015, by = 4)) +
  labs(title = "A. Raw series with LOESS smooths (span = 0.5)",
       x = NULL) +
  theme_bw() +
  theme(axis.title.y.left  = element_text(color = "coral"),
        axis.title.y.right = element_text(color = "steelblue"))

# --- Panel B: Detrended residuals ---
p_resid <- ggplot(loess_df, aes(x = year)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(y = resid_index),
            color = "coral", linewidth = 0.8) +
  geom_point(aes(y = resid_index),
             color = "coral", size = 1.5) +
  geom_line(aes(y = resid_ssb * resid_scale),
            color = "steelblue", linewidth = 0.8) +
  geom_point(aes(y = resid_ssb * resid_scale),
             color = "steelblue", size = 1.5) +
  scale_y_continuous(
    name = "Larval index residual",
    sec.axis = sec_axis(~ . / resid_scale, name = "SSB residual (mt)")
  ) +
  scale_x_continuous(breaks = seq(1984, 2015, by = 4)) +
  labs(title = "B. Detrended residuals (raw \u2212 LOESS smooth)",
       x = "Year") +
  theme_bw() +
  theme(axis.title.y.left  = element_text(color = "coral"),
        axis.title.y.right = element_text(color = "steelblue"))

# --- Combine ---
p_raw / p_resid

# --- CCF on LOESS-detrended residuals ---

# Compute residuals (reusing loess_df from the plot code)
loess_df <- ss_dln3_combined_obs %>%
  filter(year >= 1984, year <= 2015) %>%
  arrange(year) %>%
  mutate(
    loess_index = predict(loess(larvae_index ~ year, span = 0.5)),
    loess_ssb   = predict(loess(ssb ~ year, span = 0.5)),
    resid_index = larvae_index - loess_index,
    resid_ssb   = ssb - loess_ssb
  )

# CCF with bootstrap CIs on detrended residuals
ccf_boot(
  x       = loess_df$resid_index,
  y       = loess_df$resid_ssb,
  lag.max = 10,
  plot    = "Spearman",
  B       = 1000
); title(main = "CCF: Detrended Larval Index vs. Detrended SSB (LOESS residuals, span = 0.5)")


# --- Windowed CCF on LOESS residuals (1997-2012) ---

loess_window <- ss_dln3_combined_obs %>%
  filter(year >= 1999, year <= 2009) %>%
  arrange(year) %>%
  mutate(
    loess_index = predict(loess(larvae_index ~ year, span = 0.5)),
    loess_ssb   = predict(loess(ssb ~ year, span = 0.5)),
    resid_index = larvae_index - loess_index,
    resid_ssb   = ssb - loess_ssb
  )

# With only ~15 obs, pull lag.max in to 5 to avoid
# burning too many degrees of freedom at the outer lags
ccf_boot(
  x       = loess_window$resid_index,
  y       = loess_window$resid_ssb,
  lag.max = 5,
  plot    = "Spearman",
  B       = 1000
); title(main = "CCF: Detrended Larval Index vs. Detrended SSB\n1997-2012 window (LOESS residuals, span = 0.5)")

## ---- Sampling-effort figures -------------------------------------------

#ss_plot_data_year
#ss_plot_data_ms 

library(scales)

month_data <- ss_plot_data_ms %>%
  filter(group_var == "Month") %>%
  mutate(
    level = recode(as.character(level),
                   "1" = "January", "2" = "February", "3" = "March"),
    level = factor(level, levels = c("January", "February", "March"))
  )

p_month <- ggplot(month_data, aes(x = level, y = count, fill = tow_type)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  shared_fill +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Cabezon Larvae: Positive vs. Total Tows",
    subtitle = "Sampling effort by month",
    y = "Number of Tows"
  ) +
  shared_bar_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Refined fill palette
shared_fill <- scale_fill_manual(
  values = c("negative_tows" = "grey88", "positive_tows" = "#2C5F8D"),
  labels = c("negative_tows" = "Larvae Absent", "positive_tows" = "Larvae Present"),
  name = NULL
)

# Refined shared theme
shared_bar_theme <- theme_minimal(base_size = 12) +
  theme(
    strip.text         = element_text(face = "bold", size = 13, margin = margin(b = 8)),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y        = element_text(size = 10),
    axis.title.y       = element_text(size = 11, margin = margin(r = 10)),
    axis.title.x       = element_blank(),
    legend.position    = "bottom",
    legend.text        = element_text(size = 11),
    legend.key.size    = unit(0.5, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    plot.title         = element_text(face = "bold", size = 15, margin = margin(b = 4)),
    plot.subtitle      = element_text(color = "grey40", size = 10, margin = margin(b = 12)),
    plot.margin        = margin(15, 15, 10, 15)
  )

# For the year plot
p_year <- ggplot(ss_plot_data_year, aes(x = year, y = count, fill = tow_type)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  shared_fill +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Cabezon Larvae: Positive vs. Total Tows",
    subtitle = "Annual sampling effort, 1981\u20132015",
    y = "Number of Tows",
    caption = "Stacked bars show larvae-present tows (blue) against total sampling effort (grey).\nElevated effort in some years (1981, 1984, 1991, 1994, 2009, 2012) reflects expanded grid coverage or event-specific surveys."
  ) +
  shared_bar_theme +
  theme(
    plot.caption = element_text(
      hjust = 0, size = 9, color = "grey30",
      margin = margin(t = 10), lineheight = 1.1
    )
  )

p_month
p_year

#ggsave("Figures/tows-by-month.png", p_month,
#       width = 9.85, height = 6.8, dpi = 300)

#ggsave("Figures/tows-by-year.png", p_year,
#       width = 9.85, height = 6.8, dpi = 300)


obs_aggregate <- obs_plot %>%
  group_by(station_id) %>%
  summarise(
    X_m = mean(X_m, na.rm = TRUE),
    Y_m = mean(Y_m, na.rm = TRUE),
    n_visits = n(),
    n_positive = sum(larvae_100m3 > 0, na.rm = TRUE),
    
    mean_abundance = mean(larvae_100m3, na.rm = TRUE),
    mean_positive_abundance = ifelse(
      sum(larvae_100m3 > 0, na.rm = TRUE) > 0,
      mean(larvae_100m3[larvae_100m3 > 0], na.rm = TRUE),
      NA_real_
    ),
    sum_abundance = sum(larvae_100m3, na.rm = TRUE),
    max_abundance = max(larvae_100m3, na.rm = TRUE),
    
    ever_positive = any(larvae_100m3 > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    plot_abundance = mean_positive_abundance
  )

# Compute padded plot limits

library(rnaturalearth)

# Get both US and Mexico landmass
na_coast <- ne_countries(
  scale = "medium",
  country = c("United States of America", "Mexico"),
  returnclass = "sf"
)

x_range <- range(obs_aggregate$X_m, na.rm = TRUE)
y_range <- range(obs_aggregate$Y_m, na.rm = TRUE)
x_buffer <- diff(x_range) * 0.10
y_buffer <- diff(y_range) * 0.10

xlim_padded <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
ylim_padded <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
library(marmap)
library(ggplot2)
library(sf)
library(dplyr)

# Fetch bathymetry
bathy <- getNOAA.bathy(
  lon1 = -122, lon2 = -116.5,
  lat1 = 31, lat2 = 35.5,
  resolution = 1
)
bathy_df <- fortify.bathy(bathy)

# Color palette
ocean_color <- "#CFE3F2"
land_color  <- "#E8DEC3"
contour_color <- "grey50"

# Region label coordinates
# Define in lon/lat, then transform to the UTM CRS
labels_sf <- st_as_sf(
  data.frame(
    name = c("California", "Mexico"),
    lon = c(-118, -117),
    lat = c(34.8, 32.4)
  ),
  coords = c("lon", "lat"),
  crs = 4326
) |>
  st_transform(st_crs(utm_crs))

label_coords <- st_coordinates(labels_sf)
X_california_label <- label_coords[1, "X"]
Y_california_label <- label_coords[1, "Y"]
X_mexico_label     <- label_coords[2, "X"]
Y_mexico_label     <- label_coords[2, "Y"]

# Transform bathymetry contour to UTM
# Extract just the 200m contour line as sf, transform to UTM
bathy_sf <- bathy_df |>
  filter(z >= -1000 & z <= 0) |>  # trim to relevant depth range for speed
  st_as_sf(coords = c("x", "y"), crs = 4326)

# Simpler approach: draw the contour in lon/lat space and let coord_sf reproject
# We do this by converting bathy to a stars/terra raster, extracting the contour
# as an sf LINESTRING, then transforming.

library(terra)
bathy_rast <- marmap::as.raster(bathy)
bathy_terra <- terra::rast(bathy_rast)
contour_levels <- c(-200, -1000, -2000)
contour_lines <- terra::as.contour(bathy_terra, levels = contour_levels) |>
  sf::st_as_sf() |>
  sf::st_set_crs(4326) |>
  sf::st_transform(st_crs(utm_crs))

# Plot
p_effort_grid <- ggplot() +
  # 200m isobath
  geom_sf(
    data = contour_lines,
    aes(alpha = factor(level)),
    color = contour_color,
    linewidth = 0.3,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  # Land
  geom_sf(
    data = na_coast,
    inherit.aes = FALSE,
    fill = land_color,
    color = "grey40",
    linewidth = 0.3
  ) +
  # CalCOFI grid lines
  geom_sf(
    data = calcofi_lines,
    inherit.aes = FALSE,
    color = "grey75",
    linewidth = 0.25
  ) +
  # Negative tows (small open circles for full sampling footprint)
  geom_point(
    data = obs_aggregate %>% filter(!ever_positive),
    aes(x = X_m, y = Y_m),
    size = 0.6,
    alpha = 0.5,
    color = "grey40",
    shape = 1
  ) +
  # Positive tows (sized filled bubbles)
  geom_point(
    data = obs_aggregate %>% filter(ever_positive),
    aes(x = X_m, y = Y_m, size = plot_abundance),
    alpha = 0.75,
    color = "#2C5F8D"
  ) +
  # Region labels
  annotate("text",
           x = X_california_label, y = Y_california_label,
           label = "California",
           color = "grey25", size = 4, fontface = "italic"
  ) +
  annotate("text",
           x = X_mexico_label, y = Y_mexico_label,
           label = "Mexico",
           color = "grey25", size = 4, fontface = "italic"
  ) +
  coord_sf(
    crs = st_crs(utm_crs),
    xlim = xlim_padded,
    ylim = ylim_padded,
    expand = FALSE
  ) +
  scale_size_continuous(
    name = "Mean larvae / 100 m\u00b3\n(positive tows)",
    range = c(1, 8)
  ) +
  labs(
    title = "Observed cabezon larval abundance by CalCOFI station",
    subtitle = "Mean density from positive tows, aggregated by station, 1981\u20132015",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = ocean_color, color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10, margin = margin(b = 10)),
    legend.position = "right"
  )

ggsave(
  "Figures/core-grid-study-region-effort.png",
  p_effort_grid,
  width = 16,
  height = 11,
  units = "in",
  dpi = 300
)


## ---- Preferred larval index figure -------------------------------------

library(scales)  # for comma formatting

p_index <- ggplot(ss_dln3_cab_index, aes(x = year, y = est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              fill = "#2C5F8D", alpha = 0.2) +
  geom_line(color = "#1A3E5C", linewidth = 0.7) +
  geom_point(color = "#1A3E5C", size = 1.5) +
  scale_x_continuous(breaks = seq(1980, 2015, by = 5)) +
  scale_y_continuous(labels = comma, 
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Estimated cabezon larval abundance index, 1981–2015",
    subtitle = "DLn-IID: Delta-lognormal model, IID encounter spatiotemporal field",
    x = NULL,
    y = "Estimated abundance index",
    caption = "Shaded band: 95% CI."
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(color = "grey40", margin = margin(b = 15)),
    plot.caption = element_text(color = "grey50", size = 11, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey92"),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text = element_text(size = 13),
    plot.margin = margin(15, 20, 15, 15)
  ) 

ggsave(
  "Figures/preferred_index.png",
  plot = p_index,
  width = 12.5,
  height = 5.5,
  units = "in",
  dpi = 300
)

## ---- Study-region map --------------------------------------------------

library(sf)
library(dplyr)
library(stringr)
library(ggplot2)
library(rnaturalearth)
library(tibble)
library(scales)

# 0. Map settings

kml_path <- "Data/CalCOFI_75StandardStations.kml"

shelf_station_max <- 60
concave_ratio     <- 0.4
zoom_pad          <- 0.18

# Increase this if red border lines still hug the coast.
# Decrease if too much of the desired border disappears.
coast_clip_km <- 6

# Remove tiny leftover red line fragments after clipping.
# Increase if you still see small hanging red pieces.
min_border_segment_km <- 8

ocean_col  <- "#CFE3F2"
land_col   <- "#E8DEC3"
grid_col   <- "grey55"
border_col <- "#B23A48"
bubble_col <- "#2C5F8D"

core_lines <- c(76.7, 80.0, 81.8, 83.3, 86.7, 90.0, 93.3)

# 1. Build land layer: US + Mexico
# This fixes the issue where Baja/Mexico reads as water.

# Label coordinates (lon/lat -> UTM)

# Helper: convert a lon/lat point to UTM (matching utm_crs)
to_utm <- function(lon, lat, crs) {
  pt <- st_sfc(st_point(c(lon, lat)), crs = 4326) %>%
    st_transform(crs)
  coords <- st_coordinates(pt)
  list(X_m = coords[1], Y_m = coords[2])
}

ca  <- to_utm(-118, 34.8, utm_crs)   # central coast / LA area
mex <- to_utm(-117, 32.4, utm_crs)   # northern Baja

california_xy <- tibble(
  X_m   = ca$X_m,
  Y_m   = ca$Y_m,
  label = "California"
)

mexico_xy <- tibble(
  X_m   = mex$X_m,
  Y_m   = mex$Y_m,
  label = "Mexico"
)

land_ne <- ne_countries(
  country = c("United States of America", "Mexico"),
  scale = 10,
  returnclass = "sf"
) %>%
  st_make_valid() %>%
  st_transform(utm_crs)

# If you want to keep the old object name:
ca_coast <- land_ne

# 2. Parse + clean CalCOFI KML

stations <- st_read(kml_path, quiet = TRUE) %>%
  st_zm(drop = TRUE) %>%
  mutate(
    nm   = as.character(Name),
    line = suppressWarnings(as.numeric(str_extract(nm, "^[0-9]+\\.[0-9]+"))),
    stn  = suppressWarnings(as.numeric(str_extract(nm, "(?<=\\s)[0-9]+\\.[0-9]+")))
  ) %>%
  filter(!is.na(line), !is.na(stn)) %>%
  distinct(nm, .keep_all = TRUE) %>%
  filter(line %in% core_lines, stn >= 26.7, stn <= 120) %>%
  st_transform(utm_crs)

# 3. CalCOFI transect lines

transects <- stations %>%
  group_by(line) %>%
  filter(n() >= 2) %>%
  arrange(stn, .by_group = TRUE) %>%
  summarise(do_union = FALSE, .groups = "drop") %>%
  st_cast("LINESTRING")

# 4. Study-region outline: open C-shape

outline_spec <- tibble::tribble(
  ~order, ~line_target, ~stn_target,
  1,      76.7,         49.0,
  2,      76.7,         60.0,
  3,      93.3,         60.0,
  4,      93.3,         26.7
)

outline_nodes <- purrr::pmap_dfr(
  outline_spec,
  function(order, line_target, stn_target) {
    stations %>%
      filter(
        dplyr::near(line, line_target, tol = 0.05),
        dplyr::near(stn,  stn_target,  tol = 0.15)
      ) %>%
      slice(1) %>%
      mutate(order = order)
  }
) %>%
  arrange(order)

if (nrow(outline_nodes) != 4) {
  print(outline_nodes %>% st_drop_geometry())
  stop("Did not find exactly 4 outline stations. Check line/station values or tolerances.")
}

outline_coords <- st_coordinates(outline_nodes)

region_outline <- st_sf(
  geometry = st_sfc(
    st_linestring(outline_coords[, 1:2]),
    crs = st_crs(utm_crs)
  )
)

# 5. Map extent from station bounding box

bbox <- st_bbox(stations)

x_range <- bbox["xmax"] - bbox["xmin"]
y_range <- bbox["ymax"] - bbox["ymin"]

xlim <- c(bbox["xmin"] - zoom_pad * x_range, bbox["xmax"] + zoom_pad * x_range)
ylim <- c(bbox["ymin"] - zoom_pad * y_range, bbox["ymax"] + zoom_pad * y_range)

p_grid <- ggplot() +
  geom_sf(
    data = ca_coast,
    inherit.aes = FALSE,
    fill = land_col,
    color = "grey55",
    linewidth = 0.45
  ) +
  
  # CalCOFI transect lines
  geom_sf(
    data = transects,
    inherit.aes = FALSE,
    color = "grey50",
    linewidth = 1.2,
    alpha = 0.9
  ) +
  
  # CalCOFI stations
  geom_sf(
    data = stations,
    inherit.aes = FALSE,
    shape = 21,
    size = 7,
    stroke = 0.65,
    fill = "white",
    color = "grey10"
  ) +
  
  # Shelf-restricted study-region outline: open C-shape
  geom_sf(
    data = region_outline,
    inherit.aes = FALSE,
    color = border_col,
    linewidth = 2,
    lineend = "round",
    linejoin = "round"
  ) +
  
  geom_text(
    data = california_xy,
    aes(x = X_m, y = Y_m, label = label),
    size = 7.0,
    fontface = "italic",
    color = "grey15",
    hjust = 0
  ) +
  
  geom_text(
    data = mexico_xy,
    aes(x = X_m, y = Y_m, label = label),
    size = 7,
    fontface = "italic",
    color = "grey15",
    hjust = 0
  ) +
  
  coord_sf(
    crs = st_crs(utm_crs),
    xlim = xlim,
    ylim = ylim,
    expand = FALSE
  ) +
  
  labs(
    title = "CalCOFI core grid and SCB-restricted study region",
    subtitle = "Core CalCOFI stations shown; study region restricted to nearshore-associated stations",
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(base_size = 18) +
  theme(
    panel.background = element_rect(fill = ocean_col, color = NA),
    panel.grid.major = element_line(color = "grey82", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey75",
      fill = NA,
      linewidth = 1.2
    ),
    
    axis.text = element_text(size = 17, color = "grey20"),
    axis.title = element_blank(),
    
    legend.position = "none",
    
    plot.title = element_text(
      face = "bold",
      size = 28,
      color = "black",
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      color = "grey30",
      size = 20,
      margin = margin(b = 14)
    ),
    plot.margin = margin(16, 18, 16, 18)
  )

p_grid

ggsave(
  "Figures/calcofi-core-grid-study-region-poster.png",
  p_grid,
  width = 16,
  height = 11,
  units = "in",
  dpi = 300
)


## ---- Larval index vs. SSB figure ---------------------------------------

library(dplyr)
library(ggplot2)
library(scales)

plot_df <- ss_dln3_combined_obs %>%
  filter(year >= 1981, year <= 2015)

scale_factor_ssb_dln3 <- max(plot_df$larvae_index, na.rm = TRUE) /
  max(plot_df$ssb, na.rm = TRUE)

p_index_ssb <- ggplot(plot_df, aes(x = year)) +
  # SSB ribbon + line (scaled to primary axis)
  geom_ribbon(
    aes(
      ymin = ssb_lo * scale_factor_ssb_dln3,
      ymax = ssb_hi * scale_factor_ssb_dln3
    ),
    fill = "coral4",
    alpha = 0.18
  ) +
  geom_line(
    aes(y = ssb * scale_factor_ssb_dln3),
    color = "coral4",
    linewidth = 1.0
  ) +
  
  # Larval index ribbon + line
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    fill = "dodgerblue4",
    alpha = 0.20
  ) +
  geom_line(
    aes(y = larvae_index),
    color = "dodgerblue4",
    linewidth = 1.0
  ) +
  
  scale_x_continuous(
    breaks = seq(1984, 2015, by = 4),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    name = "Larval abundance index",
    labels = comma,
    sec.axis = sec_axis(
      ~ . / scale_factor_ssb_dln3,
      name = "Spawning stock biomass (mt)",
      labels = comma
    )
  ) +
  labs(
    title = "Estimated larval abundance index and spawning stock biomass",
    subtitle = "Southern California cabezon, 1984–2015; ribbons show 95% uncertainty intervals",
    x = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      lineheight = 0.95,
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      size = 17,
      color = "grey30",
      linheight = 1.0,
      margin = margin(b = 14)
    ),
    axis.title.y = element_text(
      size = 19,
      face = "bold",
      color = "dodgerblue4",
      margin = margin(r = 12)
    ),
    axis.title.y.right = element_text(
      size = 19,
      face = "bold",
      color = "coral4",
      margin = margin(l = 12)
    ),
    
    axis.text.x = element_text(size = 16, color = "grey25"),
    axis.text.y = element_text(size = 16, color = "dodgerblue4"),
    axis.text.y.right = element_text(size = 16, color = "coral4"),
    
    panel.grid.major = element_line(color = "grey86", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey78",
      fill = NA,
      linewidth = 0.45
    ),
    
    plot.margin = margin(16, 18, 16, 18)
  )

p_index_ssb

ggsave(
  "Figures/idx-vs-ssb.png",
  p_index_ssb,
  width = 12.5,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

## ---- Index model-selection table ---------------------------------------


library(flextable)

aic_table <- data.frame(
  Label = c("DLn0", "DLn1", "DLn2", "DLn3", "DLn4", "DLn5", "DLn6", "DLn7"),
  Intercept = c("Global", "Global", "Global", "Global",
                "Year FE", "Year FE", "Year FE", "Global"),
  Spatial = c("on", "on", "on", "on", "on", "on", "on", "on"),
  ST_binomial = c("off", "AR1", "AR1", "IID", "off", "IID", "IID", "AR1"),
  ST_lognormal = c("off", "AR1", "IID", "off", "off", "IID", "off", "off"),
  AIC = c(
    AIC(ss_cab_fit_dln0),
    AIC(ss_cab_fit_dln1),
    AIC(ss_cab_fit_dln2),
    AIC(ss_cab_fit_dln3),
    AIC(ss_cab_fit_dln4),
    AIC(ss_cab_fit_dln5),
    AIC(ss_cab_fit_dln6),
    AIC(ss_cab_fit_dln7)
  )
)

# Compute ΔAIC and sort
aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC)
aic_table <- aic_table[order(aic_table$AIC), ]

# Round for display
aic_table$AIC <- round(aic_table$AIC, 1)
aic_table$delta_AIC <- round(aic_table$delta_AIC, 1)

# Build the flextable
ft <- flextable(aic_table) |>
  set_header_labels(
    Label = "Model",
    Intercept = "Intercept",
    Spatial = "Spatial RF",
    ST_binomial = "ST RF (binomial)",
    ST_lognormal = "ST RF (lognormal)",
    AIC = "AIC",
    delta_AIC = "\u0394AIC"
  ) |>
  align(j = c("AIC", "delta_AIC"), align = "right", part = "all") |>
  align(j = c("Intercept", "Spatial", "ST_binomial", "ST_lognormal"),
        align = "center", part = "all") |>
  bold(i = 1) |>  # bold the top (best) row
  add_footer_lines(
    values = c(
      "Model: candidate model identifier. Intercept: global intercept (\"Global\") or fixed year effects (\"Year FE\"). Spatial RF: spatial random field (always on both components). ST RF: spatiotemporal random field structure on the binomial (presence) and lognormal (positive density) components. \u0394AIC: difference from best-fit (lowest-AIC) model. All models fit using sdmTMB (Anderson et al., 2024) with delta-lognormal family."
    )
  ) |>
  autofit() |>
  fontsize(size = 10, part = "all") |>
  fontsize(size = 9, part = "footer")

ft

save_as_docx(ft, path = "Figures/model_selection_table.docx")


# -==============================================================================-
# ==== 8.  ENVIRONMENTAL COVARIATE ANALYSIS ====
# -==============================================================================-

# Comes after the index pipeline because the index-vs-environment comparisons depend on ss_dln3_cab_index and ss_dln3_combined_obs (Sec. 4-5).


# 8.0  Bottle data prep, GLMM & sdmTMB covariate models

# Reuses cabezon_shelf_spawn from Section 2 (same spawning-season + shelf
# subset, already carrying UTM coordinates, fyear, and station_id). We only
# add the rounded line/sta join keys and the annual CPUE covariate needed
# for the bottle-cast join below.
set.seed(my.seed)

# Data check
cat("Rows retained:", nrow(cabezon_shelf_spawn), "of", nrow(cabezon), "\n")
cat("Years covered:", paste(range(cabezon_shelf_spawn$year), collapse = "–"), "\n")
cat("Months present:", paste(sort(unique(cabezon_shelf_spawn$month)), collapse = ", "), "\n")
cat("Station range:", paste(range(cabezon_shelf_spawn$station), collapse = "–"), "\n")
cat("Lines present:", paste(sort(unique(cabezon_shelf_spawn$line)), collapse = ", "), "\n")

cabezon_shelf_spawn %>%
  group_by(year) %>%
  summarize(
    n_tows = n(),
    n_pos = sum(larvae_100m3 > 0),
    pct_pos = round(100 * n_pos / n_tows, 1)
  ) %>%
  print(n = Inf)

# Composite station column (combines line + station)
cabezon_shelf_spawn <- cabezon_shelf_spawn %>%
  mutate(
    line = round(line, 1),
    sta = round(station, 1)
  )

# Compute an annual cpue value to test against the response; i.e., how does larval abundance vary at a given station given overall larval abundance in a year?
annual_cpue <- cabezon_shelf_spawn %>%
  group_by(year) %>%
  summarize(mean_cpue = mean(larvae_100m3, na.rm = T))

bottle <- read_csv("Data/bottle_data.csv", locale = locale(encoding = "latin1"))

str(bottle)

# Split out bottle metadata and get unique station ID's to compare against larval data
bottle_clean <- bottle %>%
  separate(Sta_ID, into = c("line", "sta"), sep = " ", remove = F) %>%
  mutate(
    line = as.numeric(line, 1),
    sta = as.numeric(sta, 1)
  ) %>%
  mutate(
    century = as.integer(str_sub(Depth_ID, 1, 2)) * 100L,
    depth_id_core = str_split_fixed(Depth_ID, "-", n = 2)[, 2],
    year = as.integer(str_sub(depth_id_core, 1, 2)) + century,
    month = as.integer(str_sub(depth_id_core, 3, 4))
  ) %>%
  dplyr::select(-century)

bottle_clean <- bottle_clean %>%
  mutate(
    line = round(line, 1),
    sta = round(sta, 1)
  )

# Begin with environmental covariates at the surface
bottle_surface <- bottle_clean %>%
  filter(Depthm <= 10,
         month %in% c(1, 2, 3)) %>%
  group_by(line, sta, year, month) %>%
  summarize(
    temp_surf = mean(T_degC, na.rm = T),
    sal_surf = mean(Salnty, na.rm = T),
    o2_surf = mean(O2ml_L, na.rm = T),
    chla_surf = mean(ChlorA, na.rm = T),
    .groups = "drop"
  )

# Attach cleaned bottle data to larval data in one data table
cabezon_env <- cabezon_shelf_spawn %>%
  left_join(bottle_surface, by = c("line", "sta", "year", "month"))

# What's our coverage looking like?
cabezon_env %>%
  group_by(year) %>%
  summarise(
    n          = n(),
    prop_temp  = mean(!is.na(temp_surf)),
    prop_sal   = mean(!is.na(sal_surf)),
    prop_o2    = mean(!is.na(o2_surf)),
    prop_chla  = mean(!is.na(chla_surf))
  ) %>%
  print(n = Inf)

# Match new data table with associated values from larval and bottle datasets
cabezon_env %>%
  dplyr::select(temp_surf, sal_surf, o2_surf, chla_surf) %>%
  filter(complete.cases(.)) %>%
  cor() %>%
  corrplot(method = "number", type = "upper")

cabezon_env <- cabezon_env %>%
  mutate(cabezon_present = as.integer(larvae_100m3 > 0)) %>%
  mutate(cabezon_present = as.integer(cabezon_present))

table(cabezon_env$present, useNA = "ifany")

cabezon_env <- cabezon_env %>%
  mutate(sta_id = paste(line, sta, sep = "_"))

# Standardize!
cabezon_env <- cabezon_env %>%
  mutate(across(c(temp_surf, sal_surf, o2_surf, chla_surf), scale,
                .names = "{.col}_z"))

cabezon_env <- cabezon_env %>%
  left_join(annual_cpue, by = "year") %>%
  mutate(mean_cpue_z = scale(mean_cpue))

# Constrain to observations where all covariates are present
cabezon_env_complete <- cabezon_env %>%
  filter(complete.cases(temp_surf_z, sal_surf_z, o2_surf_z, mean_cpue_z))

cabezon_env_complete$fyear <- as.factor(cabezon_env_complete$year)

# Candidate models, inclusive of 1981/1983 years with incomplete data.

m1 <- glmer(cabezon_present ~ temp_surf_z + sal_surf_z + (1 | sta_id),
            data = cabezon_env_complete, family = binomial)

m2 <- glmer(cabezon_present ~ temp_surf_z + sal_surf_z + o2_surf_z + (1 | sta_id),
            data = cabezon_env_complete, family = binomial)

m3 <- glmer(cabezon_present ~ temp_surf_z + sal_surf_z + chla_surf_z + (1 | sta_id),
            data = cabezon_env_complete, family = binomial)

m4 <- glmer(cabezon_present ~ sal_surf_z + o2_surf_z + chla_surf_z + (1 | sta_id),
            data = cabezon_env_complete, family = binomial)

m5 <- glmer(cabezon_present ~ temp_surf_z + sal_surf_z + o2_surf_z + chla_surf_z + (1 | sta_id),
            data = cabezon_env_complete, family = binomial)

AIC(m1, m2, m3, m4, m5)

summary(m3)
summary(m5)

summary(m2)

if (!all(c("X", "Y") %in% names(cabezon_env_complete))) {
  cabezon_env_complete <- sdmTMB::add_utm_columns(cabezon_env_complete, ll_names = c("longitude", "latitude"))
}

cab_env_mesh <- fm_mesh_2d(
  loc = cabezon_env_complete[,c("X","Y")],
  cutoff = 20,
  max.edge = c(75, 150),
  offset = c(45, 120)
)
cab_env_mesh2 <- make_mesh(data = cabezon_env_complete, c("X", "Y"), mesh = cab_env_mesh)
plot(cab_env_mesh2); title("SPDE Mesh for Cabezon CalCOFI analysis (Constrained)")

# M1: temp + sal
m1_sdm <- sdmTMB(
  cabezon_present ~ temp_surf_z + sal_surf_z + mean_cpue_z,
  data           = cabezon_env_complete,
  mesh           = cab_env_mesh2,
  family         = binomial(),
  spatial        = "on",
  spatiotemporal = "off"
)

# M2: temp + sal + o2 (preferred inferential model)
m2_sdm <- sdmTMB(
  cabezon_present ~ temp_surf_z + sal_surf_z + o2_surf_z + mean_cpue_z,
  data           = cabezon_env_complete,
  mesh           = cab_env_mesh2,
  family         = binomial(),
  spatial        = "on",
  spatiotemporal = "off"
)

# M3: temp + sal + chla
m3_sdm <- sdmTMB(
  cabezon_present ~ temp_surf_z + sal_surf_z + chla_surf_z + mean_cpue_z,
  data           = cabezon_env_complete,
  mesh           = cab_env_mesh2,
  family         = binomial(),
  spatial        = "on",
  spatiotemporal = "off"
)

# M4: sal + o2 + chla
m4_sdm <- sdmTMB(
  cabezon_present ~ sal_surf_z + o2_surf_z + chla_surf_z + mean_cpue_z,
  data           = cabezon_env_complete,
  mesh           = cab_env_mesh2,
  family         = binomial(),
  spatial        = "on",
  spatiotemporal = "off"
)

# M5: full model
m5_sdm <- sdmTMB(
  cabezon_present ~ temp_surf_z + sal_surf_z + o2_surf_z + chla_surf_z + mean_cpue_z,
  data           = cabezon_env_complete,
  mesh           = cab_env_mesh2,
  family         = binomial(),
  spatial        = "on",
  spatiotemporal = "off"
)

# Compare
AIC(m1_sdm, m2_sdm, m3_sdm, m4_sdm, m5_sdm)

sanity(m2_sdm)
summary(m2_sdm)

cabezon_env_complete$resid_m2 <- residuals(m2_sdm, type = "response")

plot(fitted(m2_sdm), cabezon_env_complete$resid_m2,
     xlab = "Fitted values", ylab = "RQ Residuals",
     main = "Residuals vs. Fitted - m2_sdm")
abline(h = 0, lty = 2)

ggplot(cabezon_env_complete, aes(X, Y, color = resid_m2)) +
  geom_point(size = 2) +
  scale_color_gradient2() +
  facet_wrap(~year) +
  coord_fixed() +
  labs(title = "Spatial residuals by year - m2_sdm")

# What years were lost in the complete case filter?
setdiff(unique(cabezon_env$year), unique(cabezon_env_complete$year))

# Distribution of positive catches
cabezon_env_complete %>%
  filter(cabezon_present == TRUE) %>%
  pull(larvae_100m3) %>%
  summary()

# And a quick histogram
cabezon_env_complete %>%
  filter(cabezon_present == TRUE) %>%
  ggplot(aes(larvae_100m3)) +
  geom_histogram(bins = 30) +
  scale_x_log10() +
  labs(title = "Distribution of positive catches (log scale)")


m2_hurdle_ln <- sdmTMB(
  larvae_100m3 ~ temp_surf_z + sal_surf_z + o2_surf_z + mean_cpue_z,
  data = cabezon_env_complete,
  mesh = cab_env_mesh2,
  family = delta_lognormal(),
  spatial = list("on", "on"),
  spatiotemporal = "off"
)

m2_hurdle_ln_time <- sdmTMB(
  larvae_100m3 ~ temp_surf_z + fyear + sal_surf_z + o2_surf_z,
  data = cabezon_env_complete,
  time = "year",
  mesh = cab_env_mesh2,
  family = delta_lognormal(),
  spatial = list("on", "on"),
  spatiotemporal = "off"
)
sanity(m2_hurdle_ln_time)
summary(m2_hurdle_ln_time)
AIC(m2_hurdle_ln_time)
tidy(m2_hurdle_ln_time, effects = "fixed", model = 1, conf.int = TRUE)
tidy(m2_hurdle_ln_time, effects = "fixed", model = 2, conf.int = TRUE)

m2_hurdle_gm <- sdmTMB(
  larvae_100m3 ~ temp_surf_z + sal_surf_z + o2_surf_z + mean_cpue_z,
  data = cabezon_env_complete,
  mesh = cab_env_mesh2,
  family = delta_gamma(),
  spatial = "on",
  spatiotemporal = "off"
)

m2_tw <- sdmTMB(
  larvae_100m3 ~ temp_surf_z + sal_surf_z + o2_surf_z + mean_cpue_z,
  data = cabezon_env_complete,
  mesh = cab_env_mesh2,
  family = tweedie(),
  spatial = "on",
  spatiotemporal = "off"
)

AIC(m2_hurdle_ln, m2_hurdle_gm, m2_tw)

sanity(m2_hurdle_ln)
summary(m2)
summary(m2_sdm)
summary(m2_hurdle_ln)

tidy(m2_hurdle_ln, effects = "fixed", model = 1, conf.int = TRUE)
tidy(m2_hurdle_ln, effects = "fixed", model = 2, conf.int = TRUE)

# How well is mean_cpue_z (our index of annual larval prevalence) capturing a fixed year effect?
year_fx <- tidy(m2_hurdle_ln_time, model = 1, conf.int = TRUE) |>
  filter(grepl("fyear", term)) |>
  mutate(year = as.integer(gsub("fyear", "", term)))

# Get mean_cpue_z by year (one value per year)
year_cpue <- cabezon_env_complete |>
  group_by(year) |>
  summarize(mean_cpue_z = first(mean_cpue_z))

# Join and plot
year_data <- left_join(year_fx, year_cpue, by = "year")

left_join(year_fx, year_cpue, by = "year") |>
  ggplot(aes(mean_cpue_z, estimate)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_smooth(method = "lm") +
  labs(x = "Annual mean CPUE (z-scored)",
       y = "Year fixed effect (binomial component)")

fit <- lm(estimate ~ mean_cpue_z, data = year_data)
summary(fit)$r.squared

library(DHARMa)

sim_hurdle <- simulate(m2_hurdle_ln, nsim = 500, type = "mle-mvn")

dharma_hurdle <- createDHARMa(
  simulatedResponse = sim_hurdle,
  observedResponse = cabezon_env_complete$larvae_100m3,
  fittedPredictedResponse = predict(m2_hurdle_ln)$est
)

plot(dharma_hurdle)

library(visreg)

# Component 1 (binomial / encounter)
visreg_delta(m2_hurdle_ln, xvar = "temp_surf_z", model = 1, 
             scale = "response", gg = TRUE)
visreg_delta(m2_hurdle_ln, xvar = "o2_surf_z", model = 1, 
             scale = "response", gg = TRUE)

# Component 2 (lognormal / positive density) — for completeness
visreg_delta(m2_hurdle_ln, xvar = "temp_surf_z", model = 2,
             scale = "response", gg = TRUE)

library(dplyr)
library(ggplot2)

# Build a prediction grid: vary one covariate, hold others at zero (= mean, since z-scored)
nd_temp <- data.frame(
  temp_surf_z = seq(-2.5, 2.5, length.out = 100),
  sal_surf_z  = 0,
  o2_surf_z   = 0,
  mean_cpue_z = 0,
  X = mean(cabezon_env_complete$X),   # not used once re_form = NA
  Y = mean(cabezon_env_complete$Y)  # not used once re_form = NA
)

# Predict, excluding random fields so you isolate the fixed-effect response
pred_temp <- predict(m2_hurdle_ln, newdata = nd_temp, 
                     re_form = NA, re_form_iid = NA,
                     se_fit = TRUE, model = 1)  # component 1

# Back-transform from logit to probability
pred_temp$prob      <- plogis(pred_temp$est)
pred_temp$prob_lwr  <- plogis(pred_temp$est - 1.96 * pred_temp$est_se)
pred_temp$prob_upr  <- plogis(pred_temp$est + 1.96 * pred_temp$est_se)

ggplot(pred_temp, aes(temp_surf_z, prob)) +
  geom_ribbon(aes(ymin = prob_lwr, ymax = prob_upr), alpha = 0.3) +
  geom_line(linewidth = 1) +
  labs(x = "Surface temperature (z-scored)",
       y = "Predicted encounter probability",
       title = "Component 1: encounter probability vs. SST")

# 8.1  Constrained-window (1984-2014) refit  [basis for partial-effect figures]

# NOTE: the partial-effect figures below use m2_sdm_sens and cabezon_env_constrained, both defined here. Need to decide whether the reported partial effects should come from this window or the full-window m2_sdm.

# Sensitivity analysis: constrain to 1984-2014

cabezon_env_constrained <- cabezon_env_complete %>%
  filter(year >= 1984, year <= 2014)

cat("Rows retained:", nrow(cabezon_env_constrained), 
    "of", nrow(cabezon_env_complete), "\n")
cat("Years covered:", paste(range(cabezon_env_constrained$year), collapse = "–"), "\n")

# Rebuild mesh on constrained dataset
cab_env_mesh_sens <- fm_mesh_2d(
  loc = cabezon_env_constrained[, c("X", "Y")],
  cutoff = 20,
  max.edge = c(75, 150),
  offset = c(45, 120)
)
cab_env_mesh_sens2 <- make_mesh(data = cabezon_env_constrained, 
                                c("X", "Y"), 
                                mesh = cab_env_mesh_sens)


c_m1 <- glmer(cabezon_present ~ temp_surf_z + sal_surf_z + (1 | sta_id),
              data = cabezon_env_constrained, family = binomial)

c_m2 <- glmer(cabezon_present ~ temp_surf_z + sal_surf_z + o2_surf_z + (1 | sta_id),
              data = cabezon_env_constrained, family = binomial)

c_m3 <- glmer(cabezon_present ~ temp_surf_z + sal_surf_z + chla_surf_z + (1 | sta_id),
              data = cabezon_env_constrained, family = binomial)

c_m4 <- glmer(cabezon_present ~ sal_surf_z + o2_surf_z + chla_surf_z + (1 | sta_id),
              data = cabezon_env_constrained, family = binomial)

c_m5 <- glmer(cabezon_present ~ temp_surf_z + sal_surf_z + o2_surf_z + chla_surf_z + (1 | sta_id),
              data = cabezon_env_constrained, family = binomial)

AIC(c_m1, c_m2, c_m3, c_m4, c_m5)

summary(c_m2)
summary(c_m5)

# M1: temp + sal
c_m1_sdm <- sdmTMB(
  cabezon_present ~ temp_surf_z + sal_surf_z + mean_cpue_z,
  data           = cabezon_env_constrained,
  mesh           = cab_env_mesh_sens2,
  family         = binomial(),
  spatial        = "on",
  spatiotemporal = "off"
)

# M2: temp + sal + o2 (preferred inferential model)
c_m2_sdm <- sdmTMB(
  cabezon_present ~ temp_surf_z + sal_surf_z + o2_surf_z + mean_cpue_z,
  data           = cabezon_env_constrained,
  mesh           = cab_env_mesh_sens2,
  family         = binomial(),
  spatial        = "on",
  spatiotemporal = "off"
)

# M3: temp + sal + chla
c_m3_sdm <- sdmTMB(
  cabezon_present ~ temp_surf_z + sal_surf_z + chla_surf_z + mean_cpue_z,
  data           = cabezon_env_constrained,
  mesh           = cab_env_mesh_sens2,
  family         = binomial(),
  spatial        = "on",
  spatiotemporal = "off"
)

# M4: sal + o2 + chla
c_m4_sdm <- sdmTMB(
  cabezon_present ~ sal_surf_z + o2_surf_z + chla_surf_z + mean_cpue_z,
  data           = cabezon_env_constrained,
  mesh           = cab_env_mesh_sens2,
  family         = binomial(),
  spatial        = "on",
  spatiotemporal = "off"
)

# M5: full model
c_m5_sdm <- sdmTMB(
  cabezon_present ~ temp_surf_z + sal_surf_z + o2_surf_z + chla_surf_z + mean_cpue_z,
  data           = cabezon_env_constrained,
  mesh           = cab_env_mesh_sens2,
  family         = binomial(),
  spatial        = "on",
  spatiotemporal = "off"
)

# Compare
AIC(c_m1_sdm, c_m2_sdm, c_m3_sdm, c_m4_sdm, c_m5_sdm)

# Refit preferred binomial model
m2_sdm_sens <- sdmTMB(
  cabezon_present ~ temp_surf_z + sal_surf_z + o2_surf_z + mean_cpue_z,
  data           = cabezon_env_constrained,
  mesh           = cab_env_mesh_sens2,
  family         = binomial(),
  spatial        = "on",
  spatiotemporal = "off"
)

sanity(m2_sdm_sens)
summary(m2_sdm_sens)

# Refit preferred hurdle model
m2_hurdle_ln_sens <- sdmTMB(
  larvae_100m3 ~ temp_surf_z + sal_surf_z + o2_surf_z + mean_cpue_z,
  data           = cabezon_env_constrained,
  mesh           = cab_env_mesh_sens2,
  family         = delta_lognormal(),
  spatial        = list("on", "on"),
  spatiotemporal = "off"
)

m2_hurdle_ln_sens_spoff <- sdmTMB(
  larvae_100m3 ~ temp_surf_z + sal_surf_z + o2_surf_z + mean_cpue_z,
  data           = cabezon_env_constrained,
  mesh           = cab_env_mesh_sens2,
  family         = delta_lognormal(),
  spatial        = list("on", "off"),
  spatiotemporal = "off"
)

sanity(m2_hurdle_ln_sens)
summary(m2_hurdle_ln_sens)
summary(m2_hurdle_ln_sens_spoff)

tidy(m2_hurdle_ln_sens, effects = "fixed", model = 1, conf.int = TRUE)
tidy(m2_hurdle_ln_sens, effects = "fixed", model = 2, conf.int = TRUE)
tidy(m2_hurdle_ln_sens, effects = "ran_pars", model = 1, conf.int = TRUE)
tidy(m2_hurdle_ln_sens, effects = "ran_pars", model = 2, conf.int = TRUE)
tidy(m2_hurdle_ln_sens_spoff, effects = "ran_pars", model = 1, conf.int = TRUE)
tidy(m2_hurdle_ln_sens_spoff, effects = "ran_pars", model = 2, conf.int = TRUE)

## 8.2 ENVIRONMENTAL ANALYSIS FIGURES ----

# Index vs. environmental covariates (SST, SSO)

library(tidyverse)

# Summarize temp_surf by year (or year-month if you want finer resolution)
env_summary <- cabezon_env_complete %>%
  dplyr::select(year, temp_surf, sal_surf, o2_surf, chla_surf) %>%
  pivot_longer(
    cols      = c(temp_surf, sal_surf, o2_surf, chla_surf),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  filter(!is.na(value) & !is.nan(value)) %>%
  group_by(year, variable) %>%
  summarise(
    mean  = mean(value),
    sd    = sd(value),
    n     = n(),
    se    = sd / sqrt(n),
    ci_lo = mean - 1.96 * se,
    ci_hi = mean + 1.96 * se,
    .groups = "drop"
  )

env_summary %>%
  filter(variable == "temp_surf") %>%
  ggplot(aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "steelblue", alpha = 0.15) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(size = 2, color = "gray20") +
  scale_x_continuous(breaks = seq(1984, 2020, by = 4)) +
  labs(
    title = "Surface Temperature (C) Over Time",
    subtitle = "Annual means ± 1 SD",
    x = "Year",
    y = "Surface Temperature (C)",
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray50", size = 10)
  )

env_summary %>%
  filter(variable == "o2_surf") %>%
  ggplot(aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "steelblue", alpha = 0.15) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(size = 2, color = "gray20") +
  scale_x_continuous(breaks = seq(1984, 2020, by = 4)) +
  labs(
    title = "Surface Oxygen Over Time",
    subtitle = "Annual means ± 1 SD",
    x = "Year",
    y = "Surface Oxygen",
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray50", size = 10)
  )


ggplot(ss_dln3_cab_index, aes(x = year, y = est)) +
  geom_line(color = "coral", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "coral", alpha = 0.15) +
  geom_point(size = 2, color = "gray20") +
  scale_x_continuous(breaks = seq(1984, 2020, by = 4)) +
  labs(y = "Estimated abundance index", x = "Year", title = "Abundance Index: Delta-Lognormal on Global Intercept with AR1 and off Spatiotemporal Fields") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray50", size = 10)
  )


# Match oxygen's range to the index's range
o2_data <- env_summary %>% filter(variable == "o2_surf")

index_min <- min(ss_dln3_cab_index$est, na.rm = TRUE)
index_max <- max(ss_dln3_cab_index$est, na.rm = TRUE)
o2_min    <- min(o2_data$mean, na.rm = TRUE)
o2_max    <- max(o2_data$mean, na.rm = TRUE)

# Linear rescale: maps o2 range onto index range
rescale_o2  <- function(x) (x - o2_min) / (o2_max - o2_min) * (index_max - index_min) + index_min
unscale_o2  <- function(x) (x - index_min) / (index_max - index_min) * (o2_max - o2_min) + o2_min

ggplot() +
  # Oxygen (rescaled to index range)
  geom_ribbon(data = o2_data,
              aes(x = year, ymin = rescale_o2(mean - sd),
                  ymax = rescale_o2(mean + sd)),
              fill = "steelblue", alpha = 0.15) +
  geom_line(data = o2_data,
            aes(x = year, y = rescale_o2(mean)),
            color = "steelblue", linewidth = 0.8) +
  geom_point(data = o2_data,
             aes(x = year, y = rescale_o2(mean)),
             size = 2, color = "steelblue4") +
  # Larval index (primary axis)
  geom_ribbon(data = ss_dln3_cab_index,
              aes(x = year, ymin = lwr, ymax = upr),
              fill = "coral", alpha = 0.15) +
  geom_line(data = ss_dln3_cab_index,
            aes(x = year, y = est),
            color = "coral", linewidth = 0.8) +
  geom_point(data = ss_dln3_cab_index,
             aes(x = year, y = est),
             size = 2, color = "coral4") +
  scale_y_continuous(
    name = "Estimated Larval Abundance Index (Delta-Lognormal AR1)",
    sec.axis = sec_axis(~ unscale_o2(.), name = "Surface Oxygen (mean ± 1 SD)")
  ) +
  scale_x_continuous(breaks = seq(1984, 2020, by = 4)) +
  labs(
    title    = "Larval Abundance Index vs. Surface Oxygen",
    subtitle = "Coral: DG-AR1 larval index  |  Blue: surface O₂ annual mean ± 1 SD",
    x = "Year"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(color = "gray50", size = 10),
    axis.title.y.right = element_text(color = "steelblue"),
    axis.text.y.right  = element_text(color = "steelblue"),
    axis.title.y.left  = element_text(color = "coral4"),
    axis.text.y.left   = element_text(color = "coral4")
  )


temp_data <- env_summary %>% filter(variable == "temp_surf")

index_min  <- min(ss_dln3_cab_index$est, na.rm = TRUE)
index_max  <- max(ss_dln3_cab_index$est, na.rm = TRUE)
temp_min   <- min(temp_data$mean, na.rm = TRUE)
temp_max   <- max(temp_data$mean, na.rm = TRUE)

rescale_temp <- function(x) (x - temp_min) / (temp_max - temp_min) * (index_max - index_min) + index_min
unscale_temp <- function(x) (x - index_min) / (index_max - index_min) * (temp_max - temp_min) + temp_min

ggplot() +
  # Temperature (rescaled to index range)
  geom_ribbon(data = temp_data,
              aes(x = year, ymin = rescale_temp(mean - sd),
                  ymax = rescale_temp(mean + sd)),
              fill = "steelblue", alpha = 0.15) +
  geom_line(data = temp_data,
            aes(x = year, y = rescale_temp(mean)),
            color = "steelblue", linewidth = 0.8) +
  geom_point(data = temp_data,
             aes(x = year, y = rescale_temp(mean)),
             size = 2, color = "steelblue4") +
  # Larval index (primary axis)
  geom_ribbon(data = ss_dln3_cab_index,
              aes(x = year, ymin = lwr, ymax = upr),
              fill = "coral", alpha = 0.15) +
  geom_line(data = ss_dln3_cab_index,
            aes(x = year, y = est),
            color = "coral", linewidth = 0.8) +
  geom_point(data = ss_dln3_cab_index,
             aes(x = year, y = est),
             size = 2, color = "coral4") +
  scale_y_continuous(
    name = "Estimated Larval Abundance Index (Full AR1 Delta-Gamma Model)",
    sec.axis = sec_axis(~ unscale_temp(.), name = "Surface Temperature (°C, mean ± 1 SD)")
  ) +
  scale_x_continuous(breaks = seq(1984, 2020, by = 4)) +
  labs(
    title    = "Larval Abundance Index vs. Surface Temperature",
    subtitle = "Coral: DLn-AR1 larval index  |  Blue: surface temperature annual mean ± 1 SD",
    x = "Year"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(color = "gray50", size = 10),
    axis.title.y.right = element_text(color = "steelblue"),
    axis.text.y.right  = element_text(color = "steelblue"),
    axis.title.y.left  = element_text(color = "coral4"),
    axis.text.y.left   = element_text(color = "coral4")
  )

ccf_boot(
  x = ss_dln3_combined_obs %>% filter(year >= 1984, year <= 2015) %>% pull(larvae_index),
  y = o2_data$mean,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
)

ccf_boot(
  x = ss_dln3_combined_obs %>% filter(year >= 1984, year <= 2015) %>% pull(larvae_index),
  y = temp_data$mean,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
)

ccf_boot(
  x = diff(ss_dln3_combined_obs %>% filter(year >= 1984, year <= 2015) %>% pull(larvae_index)),
  y = diff(o2_data$mean),
  lag.max = 10,
  plot = "Spearman",
  B = 1000
)

ggplot(cabezon_env_complete, aes(x = o2_surf_z, y = temp_surf_z, color = cabezon_present)) +
  geom_point()

ggplot(cabezon_env_complete, aes(x = o2_surf, y = temp_surf, color = cabezon_present)) +
  geom_point()



# Partial-effect figures (SST)

library(dplyr)
library(ggplot2)

# ─── 1. Build prediction grid for SST ───
# Sequence across observed range, hold other covariates at their means
sst_grid <- data.frame(
  temp_surf_z = seq(min(cabezon_env_constrained$temp_surf_z, na.rm = TRUE),
                    max(cabezon_env_constrained$temp_surf_z, na.rm = TRUE),
                    length.out = 100),
  sal_surf_z  = mean(cabezon_env_constrained$sal_surf_z, na.rm = TRUE),
  o2_surf_z   = mean(cabezon_env_constrained$o2_surf_z, na.rm = TRUE),
  mean_cpue_z = mean(cabezon_env_constrained$mean_cpue_z, na.rm = TRUE),
  # spatial coords needed by sdmTMB even if marginalized — pick any valid point
  X = mean(cabezon_env_constrained$X, na.rm = TRUE),
  Y = mean(cabezon_env_constrained$Y, na.rm = TRUE)
)

# ─── 2. Predict, marginalizing over spatial field ───
sst_pred <- predict(m2_sdm_sens, 
                    newdata = sst_grid,
                    re_form = NA,            # exclude spatial random field
                    se_fit = TRUE)

# ─── 3. Convert to response scale (binomial link is logit) ───
sst_pred <- sst_pred %>%
  mutate(
    fit_resp = plogis(est),
    lwr_resp = plogis(est - 1.96 * est_se),
    upr_resp = plogis(est + 1.96 * est_se)
  )

# ─── 4. Back-transform x-axis to original units (so the plot is interpretable) ───
sst_mean <- mean(cabezon_env_constrained$temp_surf, na.rm = TRUE)
sst_sd   <- sd(cabezon_env_constrained$temp_surf, na.rm = TRUE)
sst_pred$temp_surf <- sst_pred$temp_surf_z * sst_sd + sst_mean

# ─── 5. Plot ───
sst_partial <- ggplot(sst_pred, aes(x = temp_surf, y = fit_resp)) +
  geom_ribbon(aes(ymin = lwr_resp, ymax = upr_resp), 
              fill = "dodgerblue4", alpha = 0.25) +
  geom_line(color = "dodgerblue4", linewidth = 0.9) +
  labs(x = "Sea surface temperature (°C)",
       y = "Predicted probability of larval presence",
       title = "Partial effect of SST") +
  theme_minimal()

sst_partial

library(ggplot2)
library(scales)

sst_partial <- ggplot(sst_pred, aes(x = temp_surf, y = fit_resp)) +
  geom_ribbon(
    aes(ymin = lwr_resp, ymax = upr_resp),
    fill = "dodgerblue4",
    alpha = 0.22
  ) +
  geom_line(
    color = "dodgerblue4",
    linewidth = 1.6
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Sea surface temperature and larval occurrence",
    subtitle = "Partial effect from presence/absence model; other covariates held at their means",
    x = "Sea surface temperature (°C)",
    y = "Predicted probability of larval presence"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 26,
      color = "black",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 18,
      color = "grey30",
      margin = margin(b = 14)
    ),
    
    axis.title.x = element_text(
      size = 20,
      face = "bold",
      margin = margin(t = 12)
    ),
    axis.title.y = element_text(
      size = 20,
      face = "bold",
      margin = margin(r = 12)
    ),
    axis.text = element_text(
      size = 17,
      color = "grey20"
    ),
    
    panel.grid.major = element_line(
      color = "grey86",
      linewidth = 0.45
    ),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey78",
      fill = NA,
      linewidth = 0.45
    ),
    
    plot.margin = margin(16, 18, 16, 18)
  )

sst_partial

ggsave(
  "Figures/sst-partial-effect-poster.png",
  sst_partial,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)

# Partial-effect figures (SSO)

library(dplyr)
library(ggplot2)

# ─── 1. Build prediction grid for SST ───
# Sequence across observed range, hold other covariates at their means
sso_grid <- data.frame(
  o2_surf_z = seq(min(cabezon_env_constrained$o2_surf_z, na.rm = TRUE),
                  max(cabezon_env_constrained$o2_surf_z, na.rm = TRUE),
                  length.out = 100),
  sal_surf_z  = mean(cabezon_env_constrained$sal_surf_z, na.rm = TRUE),
  temp_surf_z   = mean(cabezon_env_constrained$temp_surf_z, na.rm = TRUE),
  mean_cpue_z = mean(cabezon_env_constrained$mean_cpue_z, na.rm = TRUE),
  # spatial coords needed by sdmTMB even if marginalized — pick any valid point
  X = mean(cabezon_env_constrained$X, na.rm = TRUE),
  Y = mean(cabezon_env_constrained$Y, na.rm = TRUE)
)

# ─── 2. Predict, marginalizing over spatial field ───
sso_pred <- predict(m2_sdm_sens, 
                    newdata = sso_grid,
                    re_form = NA,            # exclude spatial random field
                    se_fit = TRUE)

# ─── 3. Convert to response scale (binomial link is logit) ───
sso_pred <- sso_pred %>%
  mutate(
    fit_resp = plogis(est),
    lwr_resp = plogis(est - 1.96 * est_se),
    upr_resp = plogis(est + 1.96 * est_se)
  )

# ─── 4. Back-transform x-axis to original units (so the plot is interpretable) ───
sso_mean <- mean(cabezon_env_constrained$o2_surf, na.rm = TRUE)
sso_sd   <- sd(cabezon_env_constrained$o2_surf, na.rm = TRUE)
sso_pred$o2_surf <- sso_pred$o2_surf_z * sso_sd + sso_mean

# ─── 5. Plot ───
sso_partial <- ggplot(sso_pred, aes(x = o2_surf, y = fit_resp)) +
  geom_ribbon(aes(ymin = lwr_resp, ymax = upr_resp), 
              fill = "darkolivegreen", alpha = 0.25) +
  geom_line(color = "darkolivegreen", linewidth = 0.9) +
  labs(x = "Sea surface dissolved oxygen",
       y = "Predicted probability of larval presence",
       title = "Partial effect of SSO") +
  theme_minimal()

sso_partial


# Partial-effect figures (combined, paper)

library(patchwork)
library(ggplot2)
library(scales)

# ─── SST partial plot ───
sst_partial <- ggplot(sst_pred, aes(x = temp_surf, y = fit_resp)) +
  geom_ribbon(aes(ymin = lwr_resp, ymax = upr_resp),
              fill = "dodgerblue4", alpha = 0.25) +
  geom_line(color = "dodgerblue4", linewidth = 0.9) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.85),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    x = "Sea surface temperature (°C)",
    y = "Predicted probability of larval presence",
    title = "Partial effect of SST"
  ) +
  theme_minimal()

# ─── SSO partial plot ───
sso_partial <- ggplot(sso_pred, aes(x = o2_surf, y = fit_resp)) +
  geom_ribbon(aes(ymin = lwr_resp, ymax = upr_resp),
              fill = "darkolivegreen", alpha = 0.25) +
  geom_line(color = "darkolivegreen", linewidth = 0.9) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.85),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    x = expression("Sea surface dissolved oxygen (mL L"^{-1}*")"),
    y = "Predicted probability of larval presence",
    title = "Partial effect of SSO"
  ) +
  theme_minimal()

# ─── Combine with patchwork ───
library(patchwork)
partial_effects <- sst_partial + sso_partial +
  plot_annotation(
    title   = "Partial effects of environmental covariates on cabezon larval presence",
    caption = "Predicted probability of larval detection across the observed range of each covariate, with all other covariates held at their means and spatial random field set to zero."
  )

partial_effects

ggsave("Figures/partial_effects_paper.png", partial_effects,
       width = 8, height = 4, units = "in", dpi = 300)

# Partial-effect figures (SSO poster)


sso_partial <- ggplot(sso_pred, aes(x = o2_surf, y = fit_resp)) +
  geom_ribbon(
    aes(ymin = lwr_resp, ymax = upr_resp),
    fill = "goldenrod4",
    alpha = 0.22
  ) +
  geom_line(
    color = "goldenrod4",
    linewidth = 1.6
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Dissolved oxygen and larval occurrence",
    subtitle = "Partial effect from presence/absence model; other covariates held at their means",
    x = "Sea surface dissolved oxygen",
    y = "Predicted probability of larval presence"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 26,
      color = "black",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 18,
      color = "grey30",
      margin = margin(b = 14)
    ),
    
    axis.title.x = element_text(
      size = 20,
      face = "bold",
      margin = margin(t = 12)
    ),
    axis.title.y = element_text(
      size = 20,
      face = "bold",
      margin = margin(r = 12)
    ),
    axis.text = element_text(
      size = 17,
      color = "grey20"
    ),
    
    panel.grid.major = element_line(
      color = "grey86",
      linewidth = 0.45
    ),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey78",
      fill = NA,
      linewidth = 0.45
    ),
    
    plot.margin = margin(16, 18, 16, 18)
  )

sso_partial

ggsave(
  "Figures/sso-partial-effect-poster.png",
  sso_partial,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)


# Partial-effect figures (faceted poster)


partial_effects <- bind_rows(
  sst_pred %>%
    transmute(
      x_value   = temp_surf,
      fit_resp,
      lwr_resp,
      upr_resp,
      effect    = "Sea surface temperature (°C)"
    ),
  sso_pred %>%
    transmute(
      x_value   = o2_surf,
      fit_resp,
      lwr_resp,
      upr_resp,
      effect    = "Sea surface dissolved oxygen (ml/L)"
    )
)%>%
  mutate(
    effect = factor(
      effect,
      levels = c(
        "Sea surface temperature (°C)",
        "Sea surface dissolved oxygen (ml/L)"
      )
    )
  )

partial_faceted <- ggplot(
  partial_effects,
  aes(x = x_value, y = fit_resp, fill = effect, color = effect)
) +
  geom_ribbon(
    aes(ymin = lwr_resp, ymax = upr_resp),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(
    linewidth = 1.5
  ) +
  facet_wrap(~ effect, scales = "free_x", nrow = 1) +
  scale_color_manual(
    values = c(
      "Sea surface temperature (°C)" = "dodgerblue4",
      "Sea surface dissolved oxygen (ml/L)" = "darkolivegreen"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Sea surface temperature (°C)" = "dodgerblue4",
      "Sea surface dissolved oxygen (ml/L)" = "darkolivegreen"
    )
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.85),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    title = "Larval occurrence varies with sea surface temperature and oxygen",
    subtitle = "Partial effects from presence/absence model; other covariates held at their means",
    x = NULL,
    y = "Predicted probability of larval presence"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 26, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 18, color = "grey40", margin = margin(b = 14)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 12)),
    axis.text = element_text(size = 16, color = "grey40"),
    
    strip.text = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey80", linewidth = 0.4),
    
    panel.grid.major = element_line(color = "grey86", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(color = "grey78", fill = NA, linewidth = 0.45),
    plot.margin = margin(16, 18, 16, 18),
    legend.position = "none"
  )

partial_faceted


ggsave(
  "Figures/partial-effects-faceted-poster.png",
  partial_faceted,
  width = 16,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

# Environmental model-comparison table


library(flextable)
library(dplyr)

aic_vals <- AIC(m1_sdm, m2_sdm, m3_sdm, m4_sdm, m5_sdm)

# Build comparison table
sdm_aic_table <- aic_vals %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Model") %>%
  mutate(
    Model = case_when(
      Model == "m1_sdm" ~ "M1",
      Model == "m2_sdm" ~ "M2",
      Model == "m3_sdm" ~ "M3",
      Model == "m4_sdm" ~ "M4",
      Model == "m5_sdm" ~ "M5"
    ),
    Covariates = case_when(
      Model == "M1" ~ "SST + SSS",
      Model == "M2" ~ "SST + SSS + SSO",
      Model == "M3" ~ "SST + SSS + SSCh",
      Model == "M4" ~ "SSS + SSO + SSCh",
      Model == "M5" ~ "SST + SSS + SSO + SSCh"
    ),
    dAIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC))
  ) %>%
  arrange(dAIC) %>%
  rename(K = df) %>%
  dplyr::select(Model, Covariates, K, AIC, dAIC, AIC_weight)

# All models include: spatial random field + mean_cpue_z offset
# noted in table footnote below

ft_env <- flextable(sdm_aic_table) %>%
  set_header_labels(
    Model      = "Model",
    Covariates = "Covariates",
    K          = "K",
    AIC        = "AIC",
    dAIC       = "\u0394AIC",
    AIC_weight = "AIC Weight"
  ) %>%
  colformat_double(j = "AIC",        digits = 2) %>%
  colformat_double(j = "dAIC",       digits = 2) %>%
  colformat_double(j = "AIC_weight", digits = 3) %>%
  bold(i = ~ dAIC == 0) %>%          # bold the preferred model row
  add_footer_lines(
    "All models include mean annual CPUE (mean_cpue_z) and a spatial random field. 
     Response variable: larval presence (binomial). K = number of estimated parameters. 
     Preferred model shown in bold."
  ) %>%
  set_table_properties(layout = "autofit") %>%
  theme_booktabs()

ft_env

save_as_docx(ft_env, path = "Figures/env_model_comparison_table.docx")



# =####################################################################=
# A.  APPENDIX — SUPPLEMENTAL / SENSITIVITY / ALTERNATIVE MODELS    ----
# =####################################################################=



# -==============================================================================-
## ==== A.1  Alternative index model: DLn fixed-year (ss_cab_fit_dln4) ====
# -==============================================================================-

# Parallel diagnostics / STAR comparison / regression for the
# fixed-year-effect delta-lognormal model. Not used by Sections 4-8.


### ---- A.1.a  Diagnostics & index -------------------------------------------

# Calibration summaries

set.seed(my.seed)
ss_dln_fy_pred_obs <- predict(ss_cab_fit_dln4, type = "response")

# Overall mean calibration
ss_dln_fy_pred_obs %>%
  summarise(
    obs_mean      = mean(larvae_100m3, na.rm = TRUE),
    pred_mean     = mean(est, na.rm = TRUE),
    obs_occ       = mean(larvae_100m3 > 0, na.rm = TRUE),
    pred_occ      = mean(est1, na.rm = TRUE),
    obs_pos_mean  = mean(larvae_100m3[larvae_100m3 > 0], na.rm = TRUE),
    pred_pos_mean = mean(est2, na.rm = TRUE)
  )

# Binomial calibration
ss_dln_fy_pred_obs$presence <- ss_dln_fy_pred_obs$larvae_100m3 > 0

ss_dln_fy_calibration_summary <- ss_dln_fy_pred_obs %>%
  mutate(prob_bin = cut(est1, breaks = seq(0, 1, by = 0.1))) %>%
  group_by(prob_bin) %>%
  summarise(
    mean_pred = mean(est1),
    observed  = mean(presence),
    n         = n()
  )

plot(ss_dln_fy_calibration_summary$mean_pred, ss_dln_fy_calibration_summary$observed,
     xlab = "Predicted presence probability",
     ylab = "Observed presence frequency")
abline(0, 1, lty = 2)

# Predicted vs. observed plots

# Raw scale
plot(
  x    = ss_dln_fy_pred_obs$est,
  y    = cabezon_shelf_spawn$larvae_100m3,
  xlab = "Predicted",
  ylab = "Observed"
)
abline(0, 1, lty = 2)

# Log scale
plot(log1p(ss_dln_fy_pred_obs$est), log1p(cabezon_shelf_spawn$larvae_100m3),
     xlab = "log(Predicted + 1)", ylab = "log(Observed + 1)")
abline(0, 1, lty = 2)

# Lognormal (positive density) component only
pos_idx <- cabezon_shelf_spawn$larvae_100m3 > 0
plot(log(ss_dln_fy_pred_obs$est2[pos_idx]), log(cabezon_shelf_spawn$larvae_100m3[pos_idx]),
     xlab = "log(Predicted positive abundance)",
     ylab = "log(Observed positive abundance)")
abline(0, 1, lty = 2)

# Randomized quantile residuals

set.seed(my.seed)
ss_dln_fy_rq_res <- residuals(ss_cab_fit_dln4, type = "mle-mvn")

# QQ plot
qqnorm(ss_dln_fy_rq_res)
qqline(ss_dln_fy_rq_res)

# Spatial residual map by year
ss_dln_fy_pred_obs$rq_resid <- ss_dln_fy_rq_res

ggplot(ss_dln_fy_pred_obs, aes(x = longitude, y = latitude, color = rq_resid)) +
  geom_point() +
  scale_color_gradient2() +
  facet_wrap(~ year)

# Temporal calibration

ss_dln_fy_pred_obs %>%
  group_by(year) %>%
  summarize(
    obs  = mean(larvae_100m3),
    pred = mean(est)
  ) %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = obs),  color = "black") +
  geom_line(aes(y = pred), color = "blue", linetype = "dashed")

# Abundance index: regular spatial grid

ss_dln_fy_pred_grid <- expand.grid(
  X = seq(min(cabezon_shelf_spawn$X), max(cabezon_shelf_spawn$X), by = 5),
  Y = seq(min(cabezon_shelf_spawn$Y), max(cabezon_shelf_spawn$Y), by = 5)
) %>%
  tidyr::crossing(year = as.integer(sort(unique(cabezon_shelf_spawn$year)))) %>%
  mutate(fyear = as.factor(year))

ss_dln_fy_cab_map   <- predict(ss_cab_fit_dln4, newdata = ss_dln_fy_pred_grid, return_tmb_object = TRUE)
ss_dln_fy_cab_index <- get_index(ss_dln_fy_cab_map, area = 1, bias_correct = TRUE)

ggplot(ss_dln_fy_cab_index, aes(x = year, y = est)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  labs(y = "Estimated abundance index", x = "Year",
       title = "Abundance Index: Delta-Lognormal Fixed Year Effect, Spatial + Spatiotemporal Off")

# Abundance index: observed station grid (preferred)

ss_dln_fy_pred_grid_obs <- cabezon_shelf_spawn %>%
  dplyr::select(X, Y, year, fyear) %>%
  distinct()

ss_dln_fy_cab_map_obs   <- predict(ss_cab_fit_dln4,
                                   newdata = ss_dln_fy_pred_grid_obs,
                                   return_tmb_object = TRUE)

ss_dln_fy_cab_index_obs <- get_index(ss_dln_fy_cab_map_obs,
                                     area = 1,
                                     bias_correct = TRUE)

ggplot(ss_dln_fy_cab_index_obs, aes(x = year, y = est)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  labs(y = "Estimated abundance index", x = "Year",
       title = "Abundance Index: Observed Station Grid (Delta-Lognormal Fixed Year)")


### ---- A.1.b  STAR comparisons -------------------------------------------

# ss_cab_fit_dln4 STAR comparisons and CCF
# Create one data frame to house all STAR + larval index data
# Note! For these comparisons, we move forward with the constrained 4.2.f. abundance index, assuming it's more stable.
ss_dln_fy_combined <- ss_dln_fy_cab_index %>%
  dplyr::select(year, est, lwr, upr) %>%
  rename(larvae_index = est) %>%
  left_join(
    STAR_SSB %>% rename(ssb = value, ssb_lo = lo, ssb_hi = hi),
    by = "year"
  ) %>%
  left_join(
    STAR_recdevs %>% rename(rec_dev = value, rec_dev_lo = lo, rec_dev_hi = hi),
    by = "year"
  ) %>%
  left_join(
    STAR_age0 %>% rename(age0 = value, age0_lo = lo, age0_hi = hi),
    by = "year"
  )

# Filter to observed years
ss_dln_fy_combined_obs <- ss_dln_fy_combined %>%
  filter(!is.na(larvae_index))

# CCF: Restricted larval index vs. age-0 recruits
ccf_boot(
  x = ss_dln_fy_combined_obs$larvae_index,
  y = ss_dln_fy_combined_obs$age0,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. Age-0 Recruits (Delta-Lognormal Fixed Year)")

# CCF: Restricted larval index vs. recruitment deviations
ccf_boot(
  x = ss_dln_fy_combined_obs$larvae_index,
  y = ss_dln_fy_combined_obs$rec_dev,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. Recruitment Deviations (Delta-Lognormal Fixed Year)")

# Nothing really :( when we try against SSB, we get some strong signal, but this is probably due to a shared declining trend. We'll try to parse that.

# CCF: Restricted larval index vs. SSB
ccf_boot(
  x = ss_dln_fy_combined_obs$larvae_index,
  y = ss_dln_fy_combined_obs$ssb,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. SSB (Delta-Lognormal Fixed Year)")

# First-differenced
ccf_boot(
  x = diff(ss_dln_fy_combined_obs$larvae_index),
  y = diff(ss_dln_fy_combined_obs$ssb),
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Δ Larval Index vs. Δ SSB (Delta-Lognormal Fixed Year)")

# Assuming ecologically that SSB would have the strongest relationship to larval abundance, we'll plot those two time series visually for a qualitative look
dln_fy_scale_factor_ssb <- max(ss_dln_fy_combined_obs$larvae_index, na.rm = TRUE) /
  max(ss_dln_fy_combined_obs$ssb, na.rm = TRUE)

ggplot(ss_dln_fy_combined_obs %>% filter(year >= 1984, year <= 2015),
       aes(x = year)) +
  # SSB scaled UP to primary axis units
  geom_ribbon(aes(ymin = ssb_lo * dln_fy_scale_factor_ssb,
                  ymax = ssb_hi * dln_fy_scale_factor_ssb),
              fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = ssb * dln_fy_scale_factor_ssb), color = "steelblue", linewidth = 0.8) +
  # Larvae index on primary axis
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "coral", alpha = 0.3) +
  geom_line(aes(y = larvae_index), color = "coral", linewidth = 0.8) +
  scale_y_continuous(
    name = "Estimated Larval Abundance Index; Delta-Lognormal Fixed Year",
    sec.axis = sec_axis(~ . / dln_fy_scale_factor_ssb,
                        name = "SSB (mt)")
  ) +
  scale_x_continuous(breaks = seq(1984, 2015, by = 4)) +
  labs(title = "SoCal Cabezon: Larval Index vs. SSB (1984–2015), Delta-Lognormal Fixed Year",
       x = "Year") +
  theme_bw()

### ---- A.1.c  Regression analysis ----------------------------------------

# RQ residuals for delta-lognormal fixed year effect model
ss_dln_fy_pred_obs %>%
  group_by(year) %>%
  mutate(med_resid = median(rq_resid)) %>%
  ggplot(aes(x = factor(year), y = rq_resid)) +
  geom_boxplot(outlier.size = 0.6, fill = "steelblue", alpha = 0.4) +
  geom_smooth(
    aes(x = as.numeric(factor(year)), y = rq_resid),
    method  = "loess", span = 0.4,
    se      = TRUE, color = "firebrick", linewidth = 0.8
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x     = "Year",
    y     = "RQ Residual",
    title = "Annual RQ Residual Distributions — ss_cab_fit_dln4"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 7))

ss_dln_fy_pred_obs %>%
  group_by(year) %>%
  summarise(med_resid = median(rq_resid), .groups = "drop") %>%
  arrange(year) %>%
  mutate(cusum = cumsum(med_resid)) %>%
  ggplot(aes(x = year, y = cusum)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    x     = "Year",
    y     = "Cumulative Sum of Median RQ Residual",
    title = "CUSUM — ss_cab_fit_dln4"
  ) +
  theme_minimal()

# 1. Is the underprediction concentrated in one delta component?
set.seed(my.seed)
ss_dln_fy_pred_obs <- ss_dln_fy_pred_obs %>%
  mutate(
    rq_resid1 = residuals(ss_cab_fit_dln4, model = 1, type = "mle-mvn"),
    rq_resid2 = residuals(ss_cab_fit_dln4, model = 2, type = "mle-mvn")
  )

ss_dln_fy_pred_obs %>%
  group_by(year) %>%
  summarise(
    med_resid1 = median(rq_resid1, na.rm = TRUE),
    med_resid2 = median(rq_resid2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year) %>%
  mutate(
    cusum1 = cumsum(replace_na(med_resid1, 0)),
    cusum2 = cumsum(replace_na(med_resid2, 0))
  ) %>%
  pivot_longer(cols = c(cusum1, cusum2),
               names_to  = "component",
               values_to = "cusum") %>%
  ggplot(aes(x = year, y = cusum, color = component)) +
  geom_line(linewidth = 0.9) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("cusum1" = "steelblue", "cusum2" = "firebrick"),
                     labels = c("cusum1" = "Delta (binomial)", "cusum2" = "Lognormal")) +
  labs(x = "Year", y = "CUSUM",
       title = "CUSUM by Delta Model Component — ss_cab_fit_dln4") +
  theme_minimal()

# 2. Is the underprediction spatially structured?
ss_dln_fy_pred_obs %>%
  group_by(station) %>%
  summarise(mean_resid = mean(rq_resid), .groups = "drop") %>%
  ggplot(aes(x = station, y = mean_resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = TRUE) +
  labs(x = "Station", y = "Mean RQ Residual",
       title = "Mean Residual by Station — ss_cab_fit_dln4") +
  theme_minimal()

#---- Simple lagged regression at lag -1

ss_dln_fy_combined_lagged <- ss_dln_fy_combined_obs %>%
  arrange(year) %>%
  left_join(
    ss_dln_fy_combined_obs %>%
      dplyr::select(year, age0_lead1 = age0, rec_dev_lead1 = rec_dev) %>%
      mutate(year = year -1),
    by = "year"
  )

ss_dln_fy_lm_age0   <- lm(age0_lead1    ~ larvae_index, data = ss_dln_fy_combined_lagged)
ss_dln_fy_lm_recdev <- lm(rec_dev_lead1 ~ larvae_index, data = ss_dln_fy_combined_lagged)

summary(ss_dln_fy_lm_age0)
summary(ss_dln_fy_lm_recdev)

ss_dln_fy_combined_lagged %>%
  filter(!is.na(age0_lead1)) %>%
  ggplot(aes(x = larvae_index, y = age0_lead1, label = year)) +
  geom_point() +
  geom_text(nudge_y = 8, size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Larval Index (year t)", y = "Age-0 Recruits (year t+1)",
       title = "Lagged Larval Index vs. Age-0 Recruits — ss_cab_fit_dln4") +
  theme_minimal()

# Sensitivity: remove 1983
ss_dln_fy_combined_lagged_no83 <- ss_dln_fy_combined_lagged %>% filter(year != 1983)
ss_dln_fy_lm_age0_no83 <- lm(age0_lead1 ~ larvae_index, data = ss_dln_fy_combined_lagged_no83)
summary(ss_dln_fy_lm_age0_no83)

# Detrended regressions
ss_dln_fy_lm_age0_detrend   <- lm(age0_lead1    ~ larvae_index + year, data = ss_dln_fy_combined_lagged)
ss_dln_fy_lm_recdev_detrend <- lm(rec_dev_lead1 ~ larvae_index + year, data = ss_dln_fy_combined_lagged)
summary(ss_dln_fy_lm_age0_detrend)
summary(ss_dln_fy_lm_recdev_detrend)

# Cook's distance
plot(ss_dln_fy_lm_age0, which = 4)

# Trimmed analysis: remove three most anomalous years
ss_dln_fy_combined_lagged_trimmed <- ss_dln_fy_combined_lagged %>%
  filter(!year %in% c(1983, 1986, 1998))

ss_dln_fy_combined_lagged_trimmed %>%
  filter(!is.na(age0_lead1)) %>%
  ggplot(aes(x = larvae_index, y = age0_lead1, label = year)) +
  geom_point() +
  geom_text(nudge_y = 8, size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Larval Index (year t)", y = "Age-0 Recruits (year t+1)",
       title = "Lagged Larval Index vs. Age-0 Recruits (trimmed) — ss_cab_fit_dln4") +
  theme_minimal()

ss_dln_fy_lm_age0_trimmed <- lm(age0_lead1 ~ larvae_index, data = ss_dln_fy_combined_lagged_trimmed)
summary(ss_dln_fy_lm_age0_trimmed)

ss_dln_fy_lm_age0_trimmed_detrend <- lm(age0_lead1 ~ larvae_index + year, data = ss_dln_fy_combined_lagged_trimmed)
summary(ss_dln_fy_lm_age0_trimmed_detrend)




# -==============================================================================-
## ==== A.2  Observed-station-grid index variants (sensitivity) ====
# -==============================================================================-

# Re-runs the STAR comparisons on the observed-station-grid indices
# instead of the regular grid. Objects are suffixed _sg so they never
# overwrite the canonical ss_*_combined objects from Sections 5-8.


### ---- A.2.a  dln3, observed-station grid --------------------------------

# ss_cab_fit_dln3 STAR comparisons and CCF
# Create one data frame to house all STAR + larval index data
# Note! For these comparisons, we move forward with the constrained 4.1.f. abundance index, assuming it's more stable.
ss_dln3_combined_sg <- ss_dln3_cab_index_obs %>%
  dplyr::select(year, est, lwr, upr) %>%
  rename(larvae_index = est) %>%
  left_join(
    STAR_SSB %>% rename(ssb = value, ssb_lo = lo, ssb_hi = hi),
    by = "year"
  ) %>%
  left_join(
    STAR_recdevs %>% rename(rec_dev = value, rec_dev_lo = lo, rec_dev_hi = hi),
    by = "year"
  ) %>%
  left_join(
    STAR_age0 %>% rename(age0 = value, age0_lo = lo, age0_hi = hi),
    by = "year"
  )
# Filter to observed years
ss_dln3_combined_obs_sg <- ss_dln3_combined_sg %>%
  filter(!is.na(larvae_index))
# CCF: Restricted larval index vs. age-0 recruits
ccf_boot(
  x = ss_dln3_combined_obs_sg$larvae_index,
  y = ss_dln3_combined_obs_sg$age0,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. Age-0 Recruits (Delta-Lognormal)")
# CCF: Restricted larval index vs. recruitment deviations
ccf_boot(
  x = ss_dln3_combined_obs_sg$larvae_index,
  y = ss_dln3_combined_obs_sg$rec_dev,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. Recruitment Deviations (Delta-Lognormal)")
# Nothing really :( when we try against SSB, we get some strong signal, but this is probably due to a shared declining trend. We'll try to parse that.
# CCF: Restricted larval index vs. SSB
ccf_boot(
  x = ss_dln3_combined_obs_sg$larvae_index,
  y = ss_dln3_combined_obs_sg$ssb,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. SSB (Delta-Lognormal)")
# First-differenced
ccf_boot(
  x = diff(ss_dln3_combined_obs_sg$larvae_index),
  y = diff(ss_dln3_combined_obs_sg$ssb),
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Δ Larval Index vs. Δ SSB (Delta-Lognormal)")
# And only differenced SSB on raw larvae
ccf_boot(
  x = ss_dln3_combined_obs_sg$larvae_index,
  y = diff(ss_dln3_combined_obs_sg$ssb),
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Larval Index vs. Δ SSB (Delta-Lognormal)")
# Join differenced SSB back to index
ssb_diff_df_dln3 <- ss_dln3_combined_obs_sg %>%
  arrange(year) %>%
  mutate(ssb_diff = c(NA, diff(ssb))) %>%
  filter(!is.na(ssb_diff))
lm_lag0_dln3 <- lm(ssb_diff ~ larvae_index, data = ssb_diff_df_dln3)
summary(lm_lag0_dln3)
# Visual
ssb_diff_df_dln3 %>%
  ggplot(aes(x = larvae_index, y = ssb_diff, label = year)) +
  geom_point() +
  geom_text(nudge_y = 50, size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Larval Index", y = "ΔSSB (mt)",
       title = "Contemporaneous larval index vs. change in SSB (Delta-Lognormal)") +
  theme_minimal()
# Assuming ecologically that SSB would have the strongest relationship to larval abundance, we'll plot those two time series visually for a qualitative look
scale_factor_ssb_dln3 <- max(ss_dln3_combined_obs_sg$larvae_index, na.rm = TRUE) /
  max(ss_dln3_combined_obs_sg$ssb, na.rm = TRUE)
ggplot(ss_dln3_combined_obs_sg %>% filter(year >= 1984, year <= 2015),
       aes(x = year)) +
  # SSB scaled UP to primary axis units
  geom_ribbon(aes(ymin = ssb_lo * scale_factor_ssb_dln3,
                  ymax = ssb_hi * scale_factor_ssb_dln3),
              fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = ssb * scale_factor_ssb_dln3), color = "steelblue", linewidth = 0.8) +
  # Larvae index on primary axis
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "coral", alpha = 0.3) +
  geom_line(aes(y = larvae_index), color = "coral", linewidth = 0.8) +
  scale_y_continuous(
    name = "Estimated Larval Abundance Index; Delta-Lognormal AR1 (encounter)",
    sec.axis = sec_axis(~ . / scale_factor_ssb_dln3,
                        name = "SSB (mt)")
  ) +
  scale_x_continuous(breaks = seq(1984, 2015, by = 4)) +
  labs(title = "SoCal Cabezon: Larval Index vs. SSB (1984–2015), Delta-Lognormal",
       x = "Year") +
  theme_bw()

### ---- A.2.b  dln4, observed-station grid --------------------------------

ss_dln_fy_combined_sg <- ss_dln_fy_cab_index_obs %>%
  dplyr::select(year, est, lwr, upr) %>%
  rename(larvae_index = est) %>%
  left_join(
    STAR_SSB %>% rename(ssb = value, ssb_lo = lo, ssb_hi = hi),
    by = "year"
  ) %>%
  left_join(
    STAR_recdevs %>% rename(rec_dev = value, rec_dev_lo = lo, rec_dev_hi = hi),
    by = "year"
  ) %>%
  left_join(
    STAR_age0 %>% rename(age0 = value, age0_lo = lo, age0_hi = hi),
    by = "year"
  )

# Filter to observed years
ss_dln_fy_combined_obs_sg <- ss_dln_fy_combined_sg %>%
  filter(!is.na(larvae_index))

# CCF: Restricted larval index vs. age-0 recruits
ccf_boot(
  x = ss_dln_fy_combined_obs_sg$larvae_index,
  y = ss_dln_fy_combined_obs_sg$age0,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. Age-0 Recruits (Delta-Lognormal Fixed Year)")

# CCF: Restricted larval index vs. recruitment deviations
ccf_boot(
  x = ss_dln_fy_combined_obs_sg$larvae_index,
  y = ss_dln_fy_combined_obs_sg$rec_dev,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. Recruitment Deviations (Delta-Lognormal Fixed Year)")

# Nothing really :( when we try against SSB, we get some strong signal, but this is probably due to a shared declining trend. We'll try to parse that.

# CCF: Restricted larval index vs. SSB
ccf_boot(
  x = ss_dln_fy_combined_obs_sg$larvae_index,
  y = ss_dln_fy_combined_obs_sg$ssb,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. SSB (Delta-Lognormal Fixed Year)")

# First-differenced
ccf_boot(
  x = diff(ss_dln_fy_combined_obs_sg$larvae_index),
  y = diff(ss_dln_fy_combined_obs_sg$ssb),
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Δ Larval Index vs. Δ SSB (Delta-Lognormal Fixed Year)")

# Assuming ecologically that SSB would have the strongest relationship to larval abundance, we'll plot those two time series visually for a qualitative look
dln_fy_scale_factor_ssb <- max(ss_dln_fy_combined_obs_sg$larvae_index, na.rm = TRUE) /
  max(ss_dln_fy_combined_obs_sg$ssb, na.rm = TRUE)

ggplot(ss_dln_fy_combined_obs_sg %>% filter(year >= 1984, year <= 2015),
       aes(x = year)) +
  # SSB scaled UP to primary axis units
  geom_ribbon(aes(ymin = ssb_lo * dln_fy_scale_factor_ssb,
                  ymax = ssb_hi * dln_fy_scale_factor_ssb),
              fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = ssb * dln_fy_scale_factor_ssb), color = "steelblue", linewidth = 0.8) +
  # Larvae index on primary axis
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "coral", alpha = 0.3) +
  geom_line(aes(y = larvae_index), color = "coral", linewidth = 0.8) +
  scale_y_continuous(
    name = "Estimated Larval Abundance Index; Delta-Lognormal Fixed Year",
    sec.axis = sec_axis(~ . / dln_fy_scale_factor_ssb,
                        name = "SSB (mt)")
  ) +
  scale_x_continuous(breaks = seq(1984, 2015, by = 4)) +
  labs(title = "SoCal Cabezon: Larval Index vs. SSB (1984–2015), Delta-Lognormal Fixed Year",
       x = "Year") +
  theme_bw()


# -==============================================================================-
## ==== A.3  Alternative family sensitivity check: Tweedie IID index ====
# -==============================================================================-

tw_iid_pred_grid <- expand.grid(
  X = seq(min(cabezon_shelf_spawn$X), max(cabezon_shelf_spawn$X), by = 5),
  Y = seq(min(cabezon_shelf_spawn$Y), max(cabezon_shelf_spawn$Y), by = 5)
) %>%
  tidyr::crossing(year = as.integer(sort(unique(cabezon_shelf_spawn$year)))) %>%
  mutate(fyear = as.factor(year))

tw_iid_pred_map <- predict(tw_iid_fit, newdata = tw_iid_pred_grid, return_tmb_object = TRUE)
tw_iid_pred_index <- get_index(tw_iid_pred_map, area = 1, bias_correct = TRUE)

library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

utm_crs <- 32611

# Predict over full grid x year
tw_iid_pred_map <- predict(
  tw_iid_fit,
  newdata = tw_iid_pred_grid,
  type = "response"
)

# Average predictions across years at each spatial grid cell
tw_iid_avg_pred_plot <- tw_iid_pred_map |>
  mutate(
    X_m = X * 1000,
    Y_m = Y * 1000,
    pred = est
  ) |>
  group_by(X, Y, X_m, Y_m) |>
  summarize(
    avg_pred = mean(pred, na.rm = TRUE),
    .groups = "drop"
  )

# California coastline
ca_coast <- ne_states(
  country = "United States of America",
  returnclass = "sf"
) |>
  dplyr::filter(name == "California") |>
  st_transform(utm_crs)

# Observed-data bounds for study-region crop
obs_plot <- cabezon_shelf_spawn |>
  mutate(
    X_m = X * 1000,
    Y_m = Y * 1000
  )

# Average prediction heatmap across all years
ggplot() +
  geom_raster(
    data = tw_iid_avg_pred_plot,
    aes(x = X_m, y = Y_m, fill = log1p(avg_pred))
  ) +
  geom_sf(
    data = ca_coast,
    inherit.aes = FALSE,
    fill = "grey80",
    color = "grey40"
  ) +
  geom_sf(
    data = calcofi_lines,
    inherit.aes = FALSE,
    color = "grey95",
    linewidth = 0.25
  ) +
  coord_sf(
    crs = st_crs(utm_crs),
    xlim = range(obs_plot$X_m, na.rm = TRUE),
    ylim = range(obs_plot$Y_m, na.rm = TRUE),
    expand = FALSE
  ) +
  scale_fill_viridis_c(
    name = "Mean predicted\nlarvae / 100 m³"
  ) +
  labs(
    title = "Average predicted cabezon larval abundance across years",
    subtitle = "Predictions averaged across years at each spatial grid cell",
    x = NULL,
    y = NULL
  ) +
  theme_minimal()

ggplot(tw_iid_pred_index, aes(x = year, y = est)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  labs(y = "Estimated abundance index", x = "Year",
       title = "Abundance Index: Tweedie on Global Intercept with IID (encounter) Spatiotemporal Field")

df <- dplyr::inner_join(ss_dln3_cab_index, tw_iid_pred_index, by = "year")  # est_dln, est_twd

cor(df$est.x, df$est.y, method = "spearman")        # primary: rank, scale- and outlier-robust
cor(log(df$est.x), log(df$est.y))                   # secondary: Pearson on log scale tames the spike

ccf_boot(
  x = tw_iid_pred_index$est,
  y = ss_dln3_cab_index$est,
  lag.max = 10,
  plot = "Spearman",
  B = 1000
); title(main = "CCF: Restricted Larval Index vs. Age-0 Recruits (Delta-Lognormal)")

funtimes::ccf_boot(df$est.x, df$est.y,
                   lag.max = 3, plot = "Spearman")

# leverage check: the influential years
cor(df$est.x[df$year != 2004], df$est.y[df$year != 2004], method = "spearman")


# -==============================================================================-
## ==== A.4  Bottle data QA/QC ====
# -==============================================================================-

# ── BOTTLE DATA QA/QC ──────────────────────────────────────────────────────────
# Builds on bottle_clean and bottle_surface from data prep script
# Checks: effort, quality flags, covariate completeness, distributions over time, station coverage, and match rate to larval data

library(tidyverse)
library(patchwork)

# 1. PARSE QUALITY FLAGS
# CalCOFI uses 9 = not measured/applicable; treat as NA for completeness checks
# Flag codes: 1 = good, 2 = good (replicate), 3 = questionable, 4 = bad, 9 = missing

bottle_flagged <- bottle_clean %>%
  filter(month %in% c(1, 2, 3),
         Depthm <= 10) %>%
  mutate(
    temp_ok  = T_qual  %in% c(NA, 1, 2),   # NA qual = assumed good in early years
    sal_ok   = S_qual  %in% c(NA, 1, 2),
    o2_ok    = O_qual  %in% c(NA, 1, 2),
    chl_ok   = Chlqua  %in% c(NA, 1, 2),
    temp_measured = !is.na(T_degC)  & T_qual  != 9,
    sal_measured  = !is.na(Salnty)  & S_qual  != 9,
    o2_measured   = !is.na(O2ml_L)  & O_qual  != 9,
    chl_measured  = !is.na(ChlorA)  & Chlqua  != 9
  )

# 2. EFFORT: CAST COUNT PER YEAR
cast_effort <- bottle_clean %>%
  filter(month %in% c(1, 2, 3)) %>%
  distinct(year, line, sta, Cst_Cnt) %>%
  group_by(year) %>%
  summarise(
    n_casts    = n_distinct(Cst_Cnt),
    n_stations = n_distinct(paste(line, sta)),
    .groups = "drop"
  )

p_effort <- ggplot(cast_effort, aes(x = year, y = n_casts)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, color = "firebrick", linewidth = 0.8) +
  labs(x = NULL, y = "N casts",
       title = "Cast effort (Jan–Mar, surface) by year") +
  theme_minimal()

# 3. COVARIATE COMPLETENESS OVER TIME
completeness <- bottle_flagged %>%
  group_by(year) %>%
  summarise(
    pct_temp = mean(temp_measured),
    pct_sal  = mean(sal_measured),
    pct_o2   = mean(o2_measured),
    pct_chl  = mean(chl_measured),
    .groups = "drop"
  ) %>%
  pivot_longer(-year, names_to = "variable", values_to = "pct_complete",
               names_prefix = "pct_")

p_complete <- ggplot(completeness, aes(x = year, y = pct_complete, color = variable)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_color_brewer(palette = "Set2",
                     labels = c("Chl-a", "O2", "Salinity", "Temperature")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey50") +
  labs(x = NULL, y = "% complete", color = NULL,
       title = "Covariate completeness by year (surface, Jan–Mar)") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 4. RAW DISTRIBUTIONS OVER TIME
# Mirrors the manta net tow presence/absence + density check
# Expectation: ENSO years (1983, 1992, 1998) show anomalies; protocol changes show jumps

bottle_surface_long <- bottle_surface %>%
  pivot_longer(cols = c(temp_surf, sal_surf, o2_surf, chla_surf),
               names_to = "variable", values_to = "value") %>%
  mutate(variable = recode(variable,
                           temp_surf  = "Temperature (°C)",
                           sal_surf   = "Salinity",
                           o2_surf    = "O2 (ml/L)",
                           chla_surf  = "Chl-a (µg/L)"
  ))

p_dist <- ggplot(bottle_surface_long, aes(x = factor(year), y = value)) +
  geom_boxplot(outlier.size = 0.4, fill = "steelblue", alpha = 0.4, linewidth = 0.3) +
  geom_smooth(aes(x = as.numeric(factor(year))),
              method = "loess", se = TRUE, color = "firebrick",
              linewidth = 0.7, alpha = 0.15) +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  labs(x = "Year", y = NULL,
       title = "Surface covariate distributions by year (Jan–Mar)",
       subtitle = "LOESS smoother; flag abrupt jumps vs. gradual ENSO-driven anomalies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 5),
        strip.text = element_text(face = "bold"))

# 5. STATION COVERAGE HEATMAP
# Which stations are consistently sampled vs. intermittently dropped?

station_grid <- bottle_clean %>%
  filter(month %in% c(1, 2, 3), Depthm <= 10) %>%
  distinct(year, line, sta) %>%
  mutate(
    sampled   = 1,
    station   = paste0(line, "_", sta)
  ) %>%
  complete(year, station, fill = list(sampled = 0)) %>%
  # Retain only stations sampled in >= 25% of years (filter noise)
  group_by(station) %>%
  filter(mean(sampled) >= 0.25) %>%
  ungroup()

p_stations <- ggplot(station_grid, aes(x = year, y = station, fill = factor(sampled))) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_manual(values = c("0" = "grey92", "1" = "steelblue"),
                    labels = c("Not sampled", "Sampled")) +
  labs(x = "Year", y = "Station (line_sta)", fill = NULL,
       title = "Bottle station coverage by year (Jan–Mar, surface)") +
  theme_minimal() +
  theme(axis.text.y  = element_text(size = 4),
        legend.position = "bottom")

# 6. SALINITY UNIT CHECK
# R_Sal appears to be in a different unit than Salnty — flag if so

sal_check <- bottle_clean %>%
  filter(month %in% c(1, 2, 3), Depthm <= 10, !is.na(Salnty), !is.na(R_Sal)) %>%
  summarise(
    Salnty_range = paste(round(min(Salnty), 2), "–", round(max(Salnty), 2)),
    R_Sal_range  = paste(round(min(R_Sal),  2), "–", round(max(R_Sal),  2)),
    Salnty_mean  = round(mean(Salnty), 3),
    R_Sal_mean   = round(mean(R_Sal),  3)
  )

cat("\n── Salinity column comparison ──\n")
print(sal_check)
cat("Note: if R_Sal >> 35, it is likely in units of 10^-3 (‰ × 10) or a raw conductance value.\n")
cat("Use Salnty (PSU) as the primary salinity column.\n\n")

# 7. MATCH RATE TO LARVAL DATA
# Requires cabezon_env_complete already in environment (from larval pipeline)
# This is the final gate: does bottle data exist where larvae were sampled?

if (exists("cabezon_env_complete")) {
  match_rate <- cabezon_env_complete %>%
    group_by(year) %>%
    summarise(
      n_tows      = n(),
      n_temp_match = sum(!is.na(temp_surf_z)),
      n_o2_match   = sum(!is.na(o2_surf_z)),
      pct_temp     = n_temp_match / n_tows,
      pct_o2       = n_o2_match   / n_tows,
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(pct_temp, pct_o2),
                 names_to = "variable", values_to = "pct_matched",
                 names_prefix = "pct_")
  
  p_match <- ggplot(match_rate, aes(x = year, y = pct_matched, color = variable)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c(temp = "steelblue", o2 = "darkorange"),
                       labels = c("O2", "Temperature")) +
    geom_hline(yintercept = 0.85, linetype = "dashed", color = "grey40") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(x = "Year", y = "Match rate", color = NULL,
         title = "Larval tow → bottle cast match rate by year",
         subtitle = "Dashed line = 85% threshold; years below worth flagging") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p_match)
} else {
  message("cabezon_env_complete not found — skipping match rate check.")
}

# 8. PRINT QC SUMMARY

effort_flag  <- cast_effort %>% filter(n_casts < 0.5 * median(n_casts))
low_complete <- completeness %>% filter(pct_complete < 0.8)

cat("── QC SUMMARY ──────────────────────────────────────────────────────────\n")
cat(sprintf("Total surface records (Jan–Mar, ≤10m): %d\n",
            nrow(filter(bottle_clean, month %in% c(1,2,3), Depthm <= 10))))
cat(sprintf("Year range: %d – %d\n",
            min(cast_effort$year), max(cast_effort$year)))
cat("\nYears with < 50% of median cast effort:\n")
print(effort_flag)
cat("\nCovariate-years with < 80% completeness:\n")
print(low_complete, n = Inf)
cat("────────────────────────────────────────────────────────────────────────\n")

# 9. COMBINED OUTPUT FIGURE
(p_effort / p_complete) | p_dist

p_effort
p_complete
p_dist


# -==============================================================================-
## ==== A.5  CCF power analysis (sieve bootstrap) ====
# -==============================================================================-

# Power analysis for the sieve-bootstrap CCF (funtimes::ccf_boot).
# Plants a known lagged correlation into surrogate series matched to the AR structure of the observed series, then measures how often ccf_boot detects it. Implementation notes:
#   * `level` is the confidence-level argument (`cl` sets the number of cores).
#   * Result columns are r_S / lower_S / upper_S (Spearman, matching the main
#     analysis).
#   * ARest() defaults to ar.method = "HVK", which can return non-stationary
#     coefficients on short surrogate draws and crash arima.sim; ar.method =
#     "yw" (Yule-Walker) is always stationary and is used here.
#   * cl = 1L forces a sequential bootstrap, avoiding cluster spawn/teardown
#     on every call.
library(funtimes)

index_ts  <- ss_dln3_combined_obs$larvae_index
assess_ts <- ss_dln3_combined_obs$rec_dev

# 1. Fit AR models so surrogates inherit each series' autocorrelation
fit_ar <- function(x) {
  x <- as.numeric(scale(x))
  m <- ar(x, order.max = 3, method = "yule-walker", aic = TRUE)
  list(phi = as.numeric(m$ar), order = m$order)
}
sim_ar <- function(phi, n, burn = 200) {
  if (length(phi) == 0) return(rnorm(n))
  as.numeric(arima.sim(list(ar = phi), n = n, n.start = burn))
}

# 2. One simulated dataset with a planted lag-L correlation r
make_pair <- function(r, n, lag, fitS, fitI) {
  S <- as.numeric(scale(sim_ar(fitS$phi, n)))
  e <- as.numeric(scale(sim_ar(fitI$phi, n)))
  I <- numeric(n)
  if (lag >= 0) {
    idx <- (lag + 1):n
    I[idx] <- r * S[1:(n - lag)] + sqrt(1 - r^2) * e[idx]
    if (lag > 0) I[1:lag] <- e[1:lag]
  } else {
    L <- -lag; idx <- 1:(n - L)
    I[idx] <- r * S[(L + 1):n] + sqrt(1 - r^2) * e[idx]
    I[(n - L + 1):n] <- e[(n - L + 1):n]
  }
  data.frame(I = I, S = S)
}

# 3. Detection via YOUR procedure: ccf_boot, check the lag-L band
detect <- function(pair, lag, B = 300, level = 0.95, ar.method = "yw") {
  cc <- tryCatch(
    ccf_boot(pair$I, pair$S,
             lag.max   = max(5, abs(lag)),
             plot      = "none",
             level     = level,
             B         = B,
             cl        = 1L,            # sequential: no cluster spawned per call
             ar.method = ar.method),    # forwarded to ARest(); "yw" = stationary
    error = function(e) NULL)
  if (is.null(cc)) return(NA)
  row <- which(cc$Lag == lag)
  obs <- cc$r_S[row]                    # Spearman, matching the main analysis
  (obs < cc$lower_S[row]) || (obs > cc$upper_S[row])
}

# 4. Power at a given r (redraws on the rare residual failure)
power_at_r <- function(r, n, lag, fitS, fitI, M = 200, B = 300,
                       level = 0.95, ar.method = "yw", max_tries = 5) {
  hits <- 0L; done <- 0L; fails <- 0L
  while (done < M) {
    d <- NA; tries <- 0L
    while (is.na(d) && tries < max_tries) {
      d <- detect(make_pair(r, n, lag, fitS, fitI), lag,
                  B = B, level = level, ar.method = ar.method)
      tries <- tries + 1L
      if (is.na(d)) fails <- fails + 1L
    }
    if (is.na(d)) next
    hits <- hits + d; done <- done + 1L
  }
  if (fails > 0) message(sprintf("  (r=%.2f: %d bootstrap refit failures, redrawn)", r, fails))
  hits / M
}

# 5. Sweep r, report power curve and minimum detectable r
run_sweep <- function(index_ts, assess_ts, lag = 0,
                      r_grid = seq(0.3, 0.8, 0.05),
                      M = 200, B = 300, level = 0.95, ar.method = "yw") {
  n    <- min(length(index_ts), length(assess_ts))
  fitI <- fit_ar(index_ts)
  fitS <- fit_ar(assess_ts)
  cat(sprintf("n = %d | index AR(%d) | assessment AR(%d), phi1 = %.2f | ar.method = %s\n",
              n, fitI$order, fitS$order,
              if (length(fitS$phi)) fitS$phi[1] else 0, ar.method))
  pw  <- sapply(r_grid, power_at_r, n = n, lag = lag,
                fitS = fitS, fitI = fitI, M = M, B = B,
                level = level, ar.method = ar.method)
  mdr <- if (any(pw >= 0.8) && any(pw < 0.8)) approx(pw, r_grid, xout = 0.8)$y else NA
  print(data.frame(r = r_grid, power = round(pw, 2)))
  cat(sprintf("Minimum detectable |r| at 80%% power: %.2f\n", mdr))
  invisible(list(grid = r_grid, power = pw, mdr = mdr))
}

# 6. Run it
set.seed(my.seed)
run_sweep(index_ts, assess_ts, lag = 0)    # repeat per assessment quantity and lag as needed
