# Purpose: Power analysis for juvenile survival studies.
# Author: Ryan N. Kinzer
# Created: 4 March 2025



# Load packages
library(tidyverse)
library(marked)

# Load functions
source('./R/simulate_cjs_data.R')
source('./R/cjs_power_analysis.R')

# Define parameters
mark_size <- seq(100, 1500, by=100)  # Varying mark-release sizes
phi <- c(0.75, 0.75, 0.50, 0.65) # Survival from marking to 1st recapture (time = 1, occ = 1)
p <- c(0.20, 0.1, 0.4, 0.3) # Detection probability at 1st recapture, 2, 3, ....)
# each juvenile facility is maybe 0.10 detection; 6 JBS sites after LGR (1 - .1)^6 = P(no detection)
# LGO, LMO, MCN and JDA are about 0.05 (P(no detection) = .81). IHD has very low detection. BON close to .15 b/c of BCC

# Run power analysis
power_results <- cjs_power_analysis(mark_size, phi, p, iterations = 100)

dat <- power_results %>%
  mutate(
    actual = case_when(
      param == 'Phi' & time == 1 ~ phi[1],
      param == 'Phi' & time == 2 ~ phi[2],
      param == 'Phi' & time == 3 ~ phi[3],
      param == 'Phi' & time == 4 ~ phi[4],      
      param == 'p' & time == 2 ~ p[1],
      param == 'p' & time == 3 ~ p[2],
      param == 'p' & time == 4 ~ p[3],
      param == 'p' & time == 5 ~ p[4]
    ),
    period = case_when(
      param == 'Phi' & time == 1 ~ '1 - Marking to ESS/ZEN',
      param == 'Phi' & time == 2 ~ '2 - ESS/ZEN to SFG',
      param == 'Phi' & time == 3 ~ '3 - SFG to Lower Granite Dam',
      param == 'Phi' & time == 4 ~ '4 - Lower Granite Dam to Downstream',
      param == 'p' & time == 2 ~ '2 - ESS/ZEN',
      param == 'p' & time == 3 ~ '3 - SFG',
      param == 'p' & time == 4 ~ '4 - Lower Granite Dam',
      param == 'p' & time == 5 ~ '5 - Downstream',
    ))

saveRDS(dat, file = './data/power_analysis/juv_survival_power.rds')






