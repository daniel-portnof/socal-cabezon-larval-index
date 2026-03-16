#Load packages
library(r4ss)
library(here)
library(dplyr)
library(ggplot2)

wd <- "~/Library/CloudStorage/SynologyDrive-Basilica/Dan's Files/SIO/0_capstone/models/cabezon/Cab_SCS_BC_STAR" # Mac

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
y <- data.frame(year=devyrs, value=devmed,lo=devlower,hi=devupper)

ggplot(data=y, aes(x=year, ymin=devlower, ymax=devupper)) +
  geom_ribbon(fill="skyblue", alpha=0.5) +
  geom_line(aes(y=devmed)) +
  labs(title="Cabezon SoCal RecDevs") +
  theme_bw()

# Plot age-0 recruits
age0yrs <- 1970:2018
startrowindex <- which(rownames(derived_quants)=="Recr_1970")
endrowindex <- which(rownames(derived_quants)=="Recr_2018")
age0 <- derived_quants[174:222, "Value"]
age0sd <- derived_quants[174:222, "StdDev"]
age0lower <- (age0-1.96*age0sd)
age0upper <- (age0+1.96*age0sd)

x <- data.frame(year=age0yrs, value=age0,lo=age0lower,hi=age0upper)

ggplot(data=x, aes(x=year, ymin=age0lower, ymax=age0upper)) +
  geom_ribbon(fill="skyblue", alpha=0.5) +
  geom_line(aes(y=age0)) +
  labs(title="Cabezon SoCal Age 0 Recruits") +
  theme_bw()

# Plot SSB
ssbyrs <- 1916:2018
startrowindex <- which(rownames(derived_quants)=="SSB_1916")
endrowindex <- which(rownames(derived_quants)=="SSB_2018")
ssb <- derived_quants[3:105, "Value"]
ssbsd <- derived_quants[3:105, "StdDev"]
ssblower <- (ssb-1.96*ssbsd)
ssbupper <- (ssb+1.96*ssbsd)

z <- data.frame(year=ssbyrs, value=ssb,lo=ssblower,hi=ssbupper)

ggplot(data=z, aes(x=year, ymin=ssblower, ymax=ssbupper)) +
  geom_ribbon(fill="skyblue", alpha=0.5) +
  geom_line(aes(y=ssb)) +
  labs(title="Cabezon SoCal SSB") +
  theme_bw()
