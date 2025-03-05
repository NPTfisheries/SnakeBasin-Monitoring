# Purpose: Power analysis for juvenile survival studies.
# Author: Ryan N. Kinzer
# Created: 4 March 2025



# Load packages
library(tidyverse)
library(marked)

# simulate CJS data
simulate_cjs_data <- function(n_marks, phi1, phi2, phi3, p2, p3, p4) {
  # n_marks: Number of marked individuals
  # phi1, phi2, phi3: survival probabilities between time steps
  # p2, p3, p4: detection probabilities at recapture events
  
  # unknown state and detection observations
  z1 <- rbinom(n_marks, 1, phi1)
  z2 <- rbinom(n_marks, 1, phi2) * z1
  z3 <- rbinom(n_marks, 1, phi3) * z2
  
  y2 <- rbinom(n_marks, 1, p2) * z1
  y3 <- rbinom(n_marks, 1, p3) * z2
  y4 <- rbinom(n_marks, 1, p4) * z3
  
  # capture histories
  ch <- paste0(1, y2, y3, y4)
  
  data <- data.frame(ch = ch, freq = 1)

  return(data)
}

power_analysis <- function(mark_size, phi1, phi2, phi3, p2, p3, p4, iterations = 100){
  df <- NULL
  df_mark <- NULL
  for(n in 1:length(mark_size)){
    for(i in 1:iterations){
      dat <- simulate_cjs_data(mark_size[n], phi1, phi2, phi3, p2, p3, p4)
      dat.proc <- process.data(dat, model = "cjs", begin.time = 1)
      dat.ddl <- make.design.data(dat.proc)
      tryCatch({
        cjs.model <- crm(dat.proc, dat.ddl, hessian = TRUE, model.parameters = list(Phi = list(formula = ~time), p = list(formula = ~time)))
        tmp <- bind_rows(cjs.model$results$reals, .id = 'param') %>%
          mutate(n = mark_size[n],
                 iter = i)
        df_mark <- bind_rows(df_mark, tmp)
      },
      error = function(e){
        return(NULL)
        }
      )
    }
    df <- bind_rows(df, df_mark)
  }
  return(df)
}

# Define parameters
mark_size <- seq(100, 1500, by=100)  # Varying mark-release sizes
phi1 <- 0.65  # Survival from marking to 1st recapture (time = 1, occ = 1)
phi2 <- 0.35 # Survival from 1st to 2nd recapture
phi3 <- .60
p2 <- 0.10    # Detection probability at 1st recapture
p3 <- 0.20    # Detection probability at 2nd recapture
p4 <- 0.50

# Run power analysis
power_results <- power_analysis(mark_size, phi1, phi2, phi3, p2, p3, p4)

dat <- power_results %>%
  mutate(
    actual = case_when(
      param == 'Phi' & time == 1 ~ phi1,
      param == 'Phi' & time == 2 ~ phi2,
      param == 'Phi' & time == 3 ~ phi3,
      param == 'p' & time == 2 ~ p2,
      param == 'p' & time == 3 ~ p3,
      param == 'p' & time == 4 ~ p4,
    ),
    period = case_when(
      param == 'Phi' & time == 1 ~ 'Marking to IPTDS',
      param == 'Phi' & time == 2 ~ 'IPTDS to Lower Granite Dam',
      param == 'Phi' & time == 3 ~ 'Lower Granite Dam to Downstream',
      param == 'p' & time == 2 ~ 'IPTDS',
      param == 'p' & time == 3 ~ 'Lower Granite Dam',
      param == 'p' & time == 4 ~ 'Downstream',
    ),
    period = factor(period, levels = c('Marking to IPTDS','IPTDS to Lower Granite Dam','Lower Granite Dam to Downstream','IPTDS','Lower Granite Dam','Downstream'))
  )


# Results
subtitle_text <- bquote("Simulated survival estimates of juveniles to IPTDS and Lower Granite Dam; " ~ 
                          phi[1] == .(phi1) ~ ", " ~ 
                          phi[2] == .(phi2) ~ ", " ~ 
                          rho[2] == .(p2) ~ ", and " ~ 
                          rho[3] == .(p3))

## survival
dat %>%
  filter(param == 'Phi',
         time %in% c(1,2)) %>%
  ggplot(aes(x = as.factor(n), y = estimate)) +
    geom_boxplot() +
    geom_hline(aes(yintercept = actual), colour = 'firebrick') +
  facet_wrap(~period) + 
  labs(subtitle = subtitle_text,
       x = 'Tags Released',
       y = 'Survival Probability')


subtitle_text <- bquote("Simulated detection estimates of juveniles at IPTDS and Lower Granite Dam; " ~ 
                          rho[2] == .(phi1) ~ ", and " ~ 
                          rho[3] == .(phi2))

## detection
dat %>%
  filter(param == 'p',
         time %in% c(2,3)) %>%
  ggplot(aes(x = as.factor(n), y = estimate)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = actual), colour = 'firebrick') +
  facet_wrap(~period) + 
  labs(subtitle = subtitle_text,
       x = 'Tags Released',
       y = 'Detection Probability')

df_sum <- dat %>%
  group_by(param, time, occ, period, n) %>%
  summarise(valid_sims = n_distinct(iter),
            actual = mean(actual),
            mu_estimate = mean(estimate),
            sd_estimate = sd(estimate),
            cv = sd_estimate/mu_estimate,
            bias = mu_estimate - actual) %>%
  ungroup()

## precision
df_sum %>%
  filter(param == 'Phi',
         time %in% c(1,2)) %>%
  ggplot(aes(x = as.factor(n), y = cv, group = period)) +
  geom_line() +
  geom_hline(yintercept = 0.15, colour = 'dodgerblue') +
  facet_wrap(~period) + 
  labs(subtitle = 'Coefficient of variation around survival estimates.',
       x = 'Tags Released',
       y = 'Coefficient of Variation')


saveRDS(dat, file = './data/power_analysis/juv_survival_power.rds')






