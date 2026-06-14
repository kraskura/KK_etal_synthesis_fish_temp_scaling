# Author: Krista Kraskura
# Description: 
# File to compare scaling slopes with Z-scores 
# ********************************************
# ********************************************

library(here)
library(tidyverse)

# colors -------
cols.amr<-c(colorBlindness::SteppedSequential5Steps[c(21, 22, 23, 24,  25)]) 
cols.rmr<-c(colorBlindness::SteppedSequential5Steps[c(11, 12, 13,14,  15)]) # "#6B990F" "#A3CC51" "#E5FFB2"
cols.fas<-c(colorBlindness::SteppedSequential5Steps[c(1,2, 3,4, 5)])  # "#990F0F" "#CC5151" "#FFB2B2"
cols.as<-c(colorBlindness::SteppedSequential5Steps[c(16,17,  18, 19, 20)]) 


cols<-c(cols.amr[1],cols.rmr[1],cols.amr[3],cols.rmr[3], cols.amr[5],cols.rmr[5], cols.as[2], cols.fas[2])
        # AMR -rmr- AMR dark - rmr dark - light - as - fas

# import data -------
scaling.params<-read.csv(here("./Data_exports/models/scaling_parameters.csv")) # main model
scaling.params.j<-read.csv(here("./Data_exports/juvenilemodels/scaling_parameters.csv")) # main
scaling.params.a<-read.csv(here("./Data_exports/adultmodels/scaling_parameters.csv")) # main model

# Written with help from Cloude Code. Sonnet 4.6 version.
# Slope + SE by term name (safer than row index) ---

z_test_slopes <- function(name1, name2, b1, se1, b2, se2) {
  z   <- (b1 - b2) / sqrt(se1^2 + se2^2)
  p   <- 2 * pnorm(abs(z), lower.tail = FALSE)
  sig <- ifelse(p < 0.001, "***",
         ifelse(p < 0.01,  "**",
         ifelse(p < 0.05,  "*",
         ifelse(p < 0.1,   ".", "ns"))))

  data.frame(
    model_1    = name1,
    model_2    = name2,
    slope_1    = b1,
    se_1       = se1,
    slope_2    = b2,
    se_2       = se2,
    difference = b1 - b2,
    z_stat     = z,
    p_value    = p,
    sig        = sig
  )
}

# multiple tests: stack into one table ---
#  Inter-specific global relationships have no comparisons to make (warm vs optimal, or MMR vs RMR) because they all have interactions in at least in temperature condition. 
results <- rbind(
  # Adults
  z_test_slopes("warmMMR", "warmRMR",
                b1 = scaling.params.a[scaling.params.a$temp_categ == "warm" &
                                        scaling.params.a$performance == "MMR", "Slope"][1],
                se1 = scaling.params.a[scaling.params.a$temp_categ == "warm" &
                                        scaling.params.a$performance == "MMR", "SE_slope"][1],
                b2 = scaling.params.a[scaling.params.a$temp_categ == "warm" &
                                        scaling.params.a$performance == "RMR", "Slope"][1],
                se2 = scaling.params.a[scaling.params.a$temp_categ == "warm" &
                                        scaling.params.a$performance == "RMR", "SE_slope"][1]), 
  # Adults  
  z_test_slopes("optimalRMR", "warmRMR",
                b1 = scaling.params.a[scaling.params.a$temp_categ == "ecol_relev" &
                                        scaling.params.a$performance == "RMR", "Slope"][1],
                se1 = scaling.params.a[scaling.params.a$temp_categ == "ecol_relev" &
                                        scaling.params.a$performance == "RMR", "SE_slope"][1],
                b2 = scaling.params.a[scaling.params.a$temp_categ == "warm" &
                                        scaling.params.a$performance == "RMR", "Slope"][1],
                se2 = scaling.params.a[scaling.params.a$temp_categ == "warm" &
                                        scaling.params.a$performance == "RMR", "SE_slope"][1]), 
  # SIGNIFICANT  
  # Adults 
  z_test_slopes("optimalFAS", "warmFAS", 
                b1 = scaling.params.a[scaling.params.a$temp_categ == "ecol_relev" &
                                        scaling.params.a$performance == "FAS", "Slope"][1],
                se1 = scaling.params.a[scaling.params.a$temp_categ == "ecol_relev" &
                                        scaling.params.a$performance == "FAS", "SE_slope"][1],
                b2 = scaling.params.a[scaling.params.a$temp_categ == "warm" &
                                        scaling.params.a$performance == "FAS", "Slope"][1],
                se2 = scaling.params.a[scaling.params.a$temp_categ == "warm" &
                                        scaling.params.a$performance == "FAS", "SE_slope"][1]), 
  # Juveniles
  z_test_slopes("warmMMR", "warmRMR",
                b1 = scaling.params.j[scaling.params.j$temp_categ == "warm" &
                                        scaling.params.j$performance == "MMR", "Slope"][1],
                se1 = scaling.params.j[scaling.params.j$temp_categ == "warm" &
                                        scaling.params.j$performance == "MMR", "SE_slope"][1],
                b2 = scaling.params.j[scaling.params.j$temp_categ == "warm" &
                                        scaling.params.j$performance == "RMR", "Slope"][1],
                se2 = scaling.params.j[scaling.params.j$temp_categ == "warm" &
                                        scaling.params.j$performance == "RMR", "SE_slope"][1]), 
  # Juveniles
  z_test_slopes("optimalRMR", "warmRMR",
                b1 = scaling.params.j[scaling.params.j$temp_categ == "ecol_relev" &
                                        scaling.params.j$performance == "RMR", "Slope"][1],
                se1 = scaling.params.j[scaling.params.j$temp_categ == "ecol_relev" &
                                        scaling.params.j$performance == "RMR", "SE_slope"][1],
                b2 = scaling.params.j[scaling.params.j$temp_categ == "warm" &
                                        scaling.params.j$performance == "RMR", "Slope"][1],
                se2 = scaling.params.j[scaling.params.j$temp_categ == "warm" &
                                        scaling.params.j$performance == "RMR", "SE_slope"][1])
  
)

print(results)

# --- optional: pretty table ---
library(knitr)
kable(results, digits = 4)

