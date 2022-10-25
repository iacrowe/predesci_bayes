library(tidyverse)
library(cmprsk)

#### Functions ####
## survival for a pair of trial participants: one treated with placebo, one with beta-blocker
lifespan_estimate_mortality <- 
  function(x, days, overallRiskIncrease, lrm_riskReduction, lrmRisk, risk_decomp, risk_hcc) {
    Pooled <- 
      x %>% 
      mutate(
        DeathRisk = OverallMortRate / 365 * overallRiskIncrease,
        decomp_risk = DeathRisk + risk_decomp / 365,
        hcc_risk = decomp_risk + risk_hcc / 365
      ) %>%
      
      mutate(
        DeathRiskRx = DeathRisk,
        decomp_risk_rx = DeathRiskRx + (risk_decomp / 365 * (1 - lrm_riskReduction)),
        hcc_risk_rx = decomp_risk_rx + (risk_hcc / 365 * (1 - lrm_riskReduction))
      )
    
    Pooled$Throw <- runif(days)
    Pooled$Throw_rx <- runif(days)
    
    Pooled$Morbid <- 
      case_when(Pooled$Throw < Pooled$DeathRisk ~ 2, #Other
                Pooled$Throw > Pooled$DeathRisk & Pooled$Throw < Pooled$decomp_risk ~ 3, #Decomp
                Pooled$Throw > Pooled$decomp_risk & Pooled$Throw < Pooled$hcc_risk ~ 4, #HCC
                TRUE ~ 0
      )
    
    Pooled$MorbidRx <- 
      case_when(Pooled$Throw_rx < Pooled$DeathRiskRx ~ 2, #Other
                Pooled$Throw_rx > Pooled$DeathRiskRx & Pooled$Throw_rx < Pooled$decomp_risk_rx ~ 3, #Decomp
                Pooled$Throw_rx > Pooled$decomp_risk_rx & Pooled$Throw_rx < Pooled$hcc_risk_rx ~ 4, #HCC
                TRUE ~ 0
      )
    
    
    OutputMorbid <- 
      Pooled %>% 
      group_by(Morbid) %>% 
      filter((Morbid > 0)) %>% 
      summarise(SurvMorbid = min(UI), .groups = 'drop') %>% 
      select(SurvMorbid, Morbid) %>%
      slice(which.min(SurvMorbid))
    
    OutputMorbidRx <- 
      Pooled %>% 
      group_by(MorbidRx) %>% 
      filter((MorbidRx > 0)) %>% 
      summarise(SurvMorbidRx = min(UI), .groups = 'drop') %>% 
      select(SurvMorbidRx, MorbidRx) %>%
      slice(which.min(SurvMorbidRx))
    
    FinalDeathOutput <- 
      bind_cols(
        OutputMorbid, OutputMorbidRx) %>%
      mutate(nsbb = lrm_riskReduction)
    
    return(FinalDeathOutput)
  }

## determine cumulative incidence using cmprsk::cuminc()
cuminc_curve_morb_plot <- function(x) {
  
  a <- cuminc(x$SurvMorbid, x$Morbid)
  
  b <- timepoints(a, times)
  
  c <- 
    as_tibble(b$est) %>%
    mutate(Cause = c("Other", "Decompensation", "HCC"))
  
  untreated <- 
    c %>% 
    pivot_longer(!Cause, names_to = "time", values_to = "est") %>% 
    mutate(
      time = as.numeric(time))
  
  untreated <- 
    untreated %>%
    pivot_wider(names_from = "Cause", values_from = "est")
  
  a <- cuminc(x$SurvMorbidRx, x$MorbidRx)
  
  b <- timepoints(a, times)
  
  c <- 
    as_tibble(b$est) %>%
    mutate(Cause = c("Other_rx", "Decompensation_rx", "HCC_rx"))
  
  treated <- 
    c %>% 
    pivot_longer(!Cause, names_to = "time", values_to = "est") %>% 
    mutate(
      time = as.numeric(time))
  
  treated <-
    treated %>%
    pivot_wider(names_from = "Cause", values_from = "est")
  
  out <- 
    left_join(untreated, treated, by = "time")
  
  return(out)
  
}


# Prepare for simulation
## read in life tables
RawMortality1990 <- read_csv("RawMortality1990.csv")
## convert to daily mortality risk
RawMortality1990Days <- RawMortality1990 %>% slice(rep(1:n(), each = 365))
## start lifetable at age 50
mortalityStart50y <- 
  RawMortality1990Days %>%
  filter(Years >=50) %>%
  rowid_to_column("UI") %>%
  select(UI, Years, OverallMortRate)
## compute number of days of life possible
n_days50 = 365 * 51


# Run Simuations
## parameters can be adjusted as required

## Paramters
set.seed(44)
n_sim = 1000 # equivalent to 1000 patients per arm in trial, reduce to speed up simulation
times <- seq(0, 3660, by = 10) # extracts data points for cumulative incidence every 10 days
treatment_effect <- 0.3 # 1 - posterior HR from moderate neutral prior
prob_decompensation <- 0.05 # 5% annual probaility of decompensation

## Single trial
## pr_decomp = 0.05
single_trial <- 
  replicate(
    n_sim,
    lifespan_estimate_mortality(
      mortalityStart50y, 
      days = n_days50,
      overallRiskIncrease = 6, # to reflect increased global mortality in cirrhosis
      lrm_riskReduction = treatment_effect,          
      lrmRisk = 0.5, # overall probability of liver-related death, Jepsen et al, Shearer et al
      risk_decomp = prob_decompensation,
      risk_hcc = 0.01),
    simplify = FALSE
  ) %>%
  bind_rows()

## extract decompensation-free lifespan gain
gain <- 
  single_trial %>% 
  mutate(
    SurvMorbid = if_else(SurvMorbid > 3650, 3650, as.numeric(SurvMorbid)),
    SurvMorbidRx = if_else(SurvMorbidRx > 3650, 3650, as.numeric(SurvMorbidRx))
  ) %>%
  summarise(untreated = mean(SurvMorbid), treated = mean(SurvMorbidRx)) %>% 
  mutate(diff = (treated - untreated) / 30)

## extract cum inc data points
single_trial_cum_inc <- 
  cuminc_curve_morb_plot(single_trial)

# Plot single trial data
## text annotations
gain_text <-
  tibble(
    gain_label = glue::glue("Mean gain in lifespan\n{round(gain[[1,3]], 1)} months per patient"),
    pr_decomp = "Annual Pr(decompensation) 0.05"
  )

arm_text <-
  tibble(
    pr_decomp = rep(c("Annual Pr(decompensation) 0.05"), each = 2),
    arm_label = rep(c("Standard of\ncare", "Beta\nblocker"), 1),
    x_pos = rep(9, 2),
    y_pos = c(max(single_trial_cum_inc$Decompensation) + 0.02, max(single_trial_cum_inc$Decompensation_rx) - 0.05),
    label_colour = rep(c("#e76254", "#1e466e"), 1)
  )


single_trial_plot <-
  ggplot(single_trial_cum_inc) +
  geom_step(
    aes(x = time / 365.25, y = Decompensation_rx),
    colour = "#1e466e"
  ) +
  geom_step(
    aes(x = time / 365.25, y = Decompensation),
    colour = "#e76254"
  ) +
  geom_text(data = gain_text, aes(label = gain_label), x = 7.5, y = 0.02) +
  geom_text(data = arm_text, aes(label = arm_label, x = x_pos, y = y_pos, colour = arm_label)) +
  scale_colour_manual(
    values = c("#1e466e", "#e76254"),
    labels = c("Beta\nblocker", "Standard of\ncare")) +
  scale_y_continuous(
    limits = c(0, max(single_trial_cum_inc$Decompensation) + 0.05),
    name = "Cumulative incidence\n(decompensation / liver-related mortality)"
  ) +
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    limits = c(0, 10),
    name = "Time (years)"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle(arm_text$pr_decomp)

single_trial_plot



