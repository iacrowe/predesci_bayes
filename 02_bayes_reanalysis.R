library(tidyverse)
library(survival)
library(brms)
library(rstan)
library(splines2)
library(bayesmeta)


## Read in data
predesci <- 
  read_csv("predesci_estimated_data.csv") %>%
  mutate(
    binary = if_else(outcome == 1, 1, 0),
    group_num = if_else(group == "bb", 1, 0)
  )

cox_fit <- coxph(Surv(time, binary) ~ group_num, data = predesci) # frequentist cause specific

brm_fit <- brm(time | cens(1 - binary) ~ group_num, data = predesci, family = "cox") # Bayesian survival model

summary(cox_fit)


## Priors
### neutral moderate
neutral_prior <- prior(normal(0, 0.355), class = "b")

neutral_model <- 
  brm(
    time | cens(1 - binary) ~ group_num,
    data = predesci,
    prior = neutral_prior,
    family = "cox",
    seed = 796,
    iter = 50000
  )

neutral_posteriors <- insight::get_parameters(neutral_model)

summary(neutral_model)


### optimistic moderate
#### parameters: OR estimate from PREDESCI power calc & fixed for Pr(OR >1) = 0.15
optimistic_prior <- prior(normal(-0.916, 0.88), class = "b") 

optimistic_model <- 
  brm(
    time | cens(1 - binary) ~ group_num,
    data = predesci,
    prior = optimistic_prior,
    family = "cox",
    seed = 796,
    iter = 50000
  )

optimistic_posteriors <- insight::get_parameters(optimistic_model)

summary(optimistic_model)


### pessimistic moderate
#### parameters: OR estimate from PREDESCI power calc & fixed for Pr(OR <1) = 0.3
pesimistic_prior <- prior(normal(0.916, 1.76), class = "b") 

pesimistic_model <- 
  brm(
    time | cens(1 - binary) ~ group_num,
    data = predesci,
    prior = pesimistic_prior,
    family = "cox",
    seed = 796,
    iter = 50000
  )

pesimistic_posteriors <- insight::get_parameters(pesimistic_model)

summary(pesimistic_model)


## Meta-analysis
y <- c(
  summary(neutral_model)$fixed[2,1],
  summary(optimistic_model)$fixed[2,1],
  summary(pesimistic_model)$fixed[2,1]
)

sigma = c(
  summary(neutral_model)$fixed[2,2],
  summary(optimistic_model)$fixed[2,2],
  summary(pesimistic_model)$fixed[2,2]
)

ll <- c("Neutral", "Optimistic", "Pesimistic")

meta_priors <- 
  bayesmeta(
    y = y, 
    sigma = sigma, 
    mu.prior = c(0, 0.355), 
    label = ll, 
    tau.porior = "DuMouchel") 

meta_rec <- forest(meta_priors, trans = exp, refline = 1, xlab = "Hazard ratio")
meta_rec <- as.data.frame(meta_rec)
meta_rec$name <- rownames(meta_rec)
meta_rec <- meta_rec[-5,]
colnames(meta_rec) <- c("Low", "Est", "High", "name")

ggplot(meta_rec, aes(x = Est, y = name)) +
  geom_point() +
  geom_errorbarh(aes(xmin = Low, xmax = High), height = 0) +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_minimal() +
  labs(
    x = "log(HR); 95% credible interval",
    y = ""
  )

summary(meta_priors)
meta_priors$I2(tau=0.316) # I2 = 0.545

### Data outputs
## Posteriors
posteriors <-
  tibble(
    neutral = neutral_posteriors[,2],
    optimistic = optimistic_posteriors[,2],
    pessimistic = pesimistic_posteriors[,2]
  )

write_csv(posteriors, "posteriors.csv")

## Summary posterior distributions
posterior_distributions <-
  tibble(
    prior = c(
      "Neutral", "Optimistic", "Pessimistic"
    ),
    mean = c(
      summary(neutral_model)$fixed[2,1],
      summary(optimistic_model)$fixed[2,1],
      summary(pesimistic_model)$fixed[2,1]
    ),
    sd = c(
      summary(neutral_model)$fixed[2,2],
      summary(optimistic_model)$fixed[2,2],
      summary(pesimistic_model)$fixed[2,2]
    )
  )

write_csv(posterior_distributions, "posterior_distributions.csv")

## Bayes meta output
bayes_meta_out <-
  meta_rec

write_csv(bayes_meta_out, "bayes_meta_out.csv")