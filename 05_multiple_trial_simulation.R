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
      mutate(
        nsbb = lrm_riskReduction)
    
    return(FinalDeathOutput)
  }

## determine relevant outcomes for plotting from tibble of input parameters
### outputs including:
#### mean gain in decompensation-free lifespan
#### significance test for cumulative incidence comparison
probab_sim_purr <- function(mort_increase, risk_reduction, pr_decomp){
  
  te_25 <- 
    replicate(
      n_sim,
      lifespan_estimate_mortality(
        mortalityStart50y, 
        days = n_days50,
        overallRiskIncrease = mort_increase, 
        lrm_riskReduction = risk_reduction,
        lrmRisk = 0.5,
        risk_decomp = pr_decomp,
        risk_hcc = 0.01),
      simplify = FALSE
    ) %>%
    bind_rows() 
  
  decomp <-
    te_25 %>%
    filter(Morbid == 2)
  
  n_decomp <- nrow(decomp)  
  
  decomp_rx <-
    te_25 %>%
    filter(MorbidRx == 2)
  
  n_decomp_rx <- nrow(decomp_rx)
  
  morb_gain <-
    te_25 %>%
    mutate(
      SurvMorbidRx = if_else(SurvMorbidRx >3650, 3650, as.double(SurvMorbidRx)),
      SurvMorbid = if_else(SurvMorbid >3650, 3650, as.double(SurvMorbid)),
      morbid_gain = SurvMorbidRx - SurvMorbid)
  
  te_25_morbid <-
    te_25 %>%
    select(SurvMorbid, Morbid) %>%
    rename(
      time = SurvMorbid,
      event  = Morbid
    ) %>%
    mutate(treatment = "none")
  
  te_25_morbid_rx <-
    te_25 %>%
    select(SurvMorbidRx, MorbidRx) %>%
    rename(
      time = SurvMorbidRx,
      event  = MorbidRx
    ) %>%
    mutate(treatment = "nsbb")
  
  te_25_comp <-
    bind_rows(te_25_morbid, te_25_morbid_rx)
  
  p_val <- cuminc(te_25_comp$time, te_25_comp$event, te_25_comp$treatment)
  
  
  output_table <-
    tibble(
      lifespan_gain_days = mean(morb_gain$morbid_gain),
      p_value = p_val$Tests[[2, 2]],
      hr_reduction = risk_reduction,
      decomp_prev = n_decomp - n_decomp_rx
    )
  
  output_table
  
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
## Paramters
set.seed(44)
n_trials = 20 # number of trials
n_sim = 200 # equivalent to 200 patients per arm in trial
times <- seq(0, 3660, by = 10)

## set up tibble for simulation
### risk reduction is sampled from posterior HR distribution per relevant prior 
new_data_skeptical <-
  tibble(
    mort_increase = rep(6, n_trials),
    risk_reduction = 1 - exp(rnorm(n_trials, mean = -0.355, sd = 0.238)),
    pr_decomp = rep(0.05, n_trials)
  )


new_data_optimistic <-
  tibble(
    mort_increase = rep(6, n_trials),
    risk_reduction = 1 - exp(rnorm(n_trials, mean = -0.916, sd = 0.313)),
    pr_decomp = rep(0.05, n_trials)
  )

## run simulations
bayes_sim_skeptical <- pmap_dfr(new_data_skeptical, probab_sim_purr)

bayes_sim_optimistic <- pmap_dfr(new_data_optimistic, probab_sim_purr)


# Plot multiple trial outcomes
## prepare data
skeptical_data <-
  bayes_sim_skeptical %>%
  mutate(prior = "skeptical")

optimistic_data <-
  bayes_sim_optimistic %>%
  mutate(prior = "optimistic")

plot_data <- 
  bind_rows(skeptical_data, optimistic_data)

## plot data - dot plot, as per manuscript
dot_overall <-
  ggplot(plot_data) +
  geom_jitter(
    aes(x = 1 - hr_reduction, y = lifespan_gain_days / 30, colour = prior)
  ) +
  scale_colour_manual(
    values = c("#376795", "#ef8a47"), 
    labels = c("optimistic", "skeptical"),
    guide = guide_legend(reverse = FALSE)) +
  annotate(
    "text", x = 0.75, y = 22, 
    label = "Optimistic prior",
    colour = "#376795") +
  annotate(
    "text", x = 1.2, y = 9, 
    label = "Neutral ('skeptical')\nprior",
    colour = "#ef8a47") +
  scale_x_continuous(
    name = "Hazard ratio\n(decompensation / liver-related mortality)",
    limits = c(0, 1.5)) +
  ylab("Mean lifespan gain at 10 years\n(months per participant)") +
  theme_classic() +
  theme(legend.position = "none") 

dot_overall


