library(tidyverse)
library(patchwork)

# Plot outputs from Bayesian re-analysis
## Read in data

posterior_dist <- read_csv("posterior_distributions.csv")
posteriors <- read_csv("posteriors.csv")


## Posteriors
post_prob <-
  posteriors %>%
  pivot_longer(
    cols = everything(),
    names_to = "prior",
    values_to = "log_hr"
  ) %>% 
  mutate(
    outcome = if_else(log_hr < 0, "benefit", "harm"),
    rope = if_else(log_hr > log(1 / 1.1) & log_hr < log(1.1), "rope", "not"),
    better_rope = if_else(log_hr < log(1 / 1.1), "clin_sig_benefit", "none"),
    outstanding = 
      case_when(
        log_hr < log(0.7) ~ "benefit",
        log_hr > log(1.3) ~ "harm",
        TRUE ~ "none")
  )

## probability summaries
prob_out <-
  post_prob %>%
  group_by(prior) %>%
  count(outcome) %>%
  mutate(prop_outcome = n / 100000) %>%
  select(-n)

rope_out <-
  post_prob %>%
  group_by(prior) %>%
  count(rope) %>%
  mutate(rope_outcome = n / 100000) %>%
  select(-n) %>% 
  filter(rope == "rope")

outstanding <-
  post_prob %>%
  group_by(prior) %>%
  count(outstanding) %>%
  mutate(outstanding_outcome = n / 100000) %>%
  select(-n)

above_rope <-
  post_prob %>%
  group_by(prior) %>%
  count(better_rope) %>%
  mutate(better_rope_outcome = n / 100000) %>%
  select(-n)

## Posterior plots
### Neutral
neutral_mean <- posterior_dist %>% filter(prior == "Neutral") %>% select(mean)
neutral_sd <- posterior_dist %>% filter(prior == "Neutral") %>% select(sd)
neutral_post_prob <- prob_out %>% filter(prior == "neutral")

neutral_plot <-
  ggplot(data.frame(x = c(-1.8, 0.7)), aes(x)) +
  stat_function(
    fun = dnorm, 
    n = 100, 
    args = list(
      mean = neutral_mean[[1,1]],
      sd = neutral_sd[[1,1]]),
    linetype = 1
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = neutral_mean[[1,1]],
      sd = neutral_sd[[1,1]]),
    fill = "#528fad",
    xlim = c(-1.8, log(1 / 1.1))
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = neutral_mean[[1,1]],
      sd = neutral_sd[[1,1]]),
    fill = "#528fad", alpha = 0.5,
    xlim = c(log(1 / 1.1), 0)
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = neutral_mean[[1,1]],
      sd = neutral_sd[[1,1]]),
    fill = "#ef8a47", alpha = 0.5,
    xlim = c(0, log(1.1))
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = neutral_mean[[1,1]],
      sd = neutral_sd[[1,1]]),
    fill = "#ef8a47",
    xlim = c(log(1.1),0.7)
  ) +
  geom_vline(
    xintercept = 0
  ) +
  # annotate(
  # "label", x = 0, y = 0.1, 
  # label = "ROPE") +
  annotate(
    "text", x = 0.4, y = 1.8, 
    label = glue::glue("Pr(harm) = {round(neutral_post_prob[[2, 3]], 2)}"),
    size = 3,
    colour = "#ef8a47") + 
  annotate(
    "text", x = -0.45, y = 1.8, 
    label = glue::glue("Pr(benefit) = {round(neutral_post_prob[[1, 3]], 2)}"),
    size = 3,
    colour = "#528fad") +
  scale_y_continuous(
    limits = c(0, 1.95),
    expand = c(0, 0)) +
  theme_classic() +
  labs(
    x = "log(HR)",
    y = "Density",
    subtitle = "Neutral ('skeptical') prior: mean = 0, sd = 0.355"
  )


### Optimistic
optimistic_mean <- posterior_dist %>% filter(prior == "Optimistic") %>% select(mean)
optimistic_sd <- posterior_dist %>% filter(prior == "Optimistic") %>% select(sd)
optimistic_post_prob <- prob_out %>% filter(prior == "optimistic")

optimistic_plot <-
  ggplot(data.frame(x = c(-1.8, 0.7)), aes(x)) +
  stat_function(
    fun = dnorm, 
    n = 100, 
    args = list(
      mean = optimistic_mean[[1,1]],
      sd = optimistic_sd[[1,1]]),
    linetype = 1
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = optimistic_mean[[1,1]],
      sd = optimistic_sd[[1,1]]),
    fill = "#528fad",
    xlim = c(-1.8, log(1 / 1.1))
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = optimistic_mean[[1,1]],
      sd = optimistic_sd[[1,1]]),
    fill = "#528fad", alpha = 0.5,
    xlim = c(log(1 / 1.1), 0)
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = optimistic_mean[[1,1]],
      sd = optimistic_sd[[1,1]]),
    fill = "#ef8a47", alpha = 0.5,
    xlim = c(0, log(1.1))
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = optimistic_mean[[1,1]],
      sd = optimistic_sd[[1,1]]),
    fill = "#ef8a47",
    xlim = c(log(1.1),0.7)
  ) +
  geom_vline(
    xintercept = 0
  ) +
  # annotate(
  # "label", x = 0, y = 0.1, 
  # label = "ROPE") +
  annotate(
    "text", x = 0.4, y = 1.8, 
    label = glue::glue("Pr(harm) = {round(optimistic_post_prob[[2, 3]], 2)}"),
    size = 3,
    colour = "#ef8a47") + 
  annotate(
    "text", x = -0.45, y = 1.8, 
    label = glue::glue("Pr(benefit) = {round(optimistic_post_prob[[1, 3]], 2)}"),
    size = 3,
    colour = "#528fad") +
  scale_y_continuous(
    limits = c(0, 1.95),
    expand = c(0, 0)) +
  theme_classic() +
  labs(
    x = "log(HR)",
    y = "Density",
    subtitle = "Optimistic prior: mean = -0.916, sd = 0.88"
  )

### Pessimistic
pessimistic_mean <- posterior_dist %>% filter(prior == "Pessimistic") %>% select(mean)
pessimistic_sd <- posterior_dist %>% filter(prior == "Pessimistic") %>% select(sd)
pessimistic_post_prob <- prob_out %>% filter(prior == "pessimistic")

pessimistic_plot <-
  ggplot(data.frame(x = c(-1.8, 0.7)), aes(x)) +
  stat_function(
    fun = dnorm, 
    n = 100, 
    args = list(
      mean = pessimistic_mean[[1,1]],
      sd = pessimistic_sd[[1,1]]),
    linetype = 1
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = pessimistic_mean[[1,1]],
      sd = pessimistic_sd[[1,1]]),
    fill = "#528fad",
    xlim = c(-1.8, log(1 / 1.1))
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = pessimistic_mean[[1,1]],
      sd = pessimistic_sd[[1,1]]),
    fill = "#528fad", alpha = 0.5,
    xlim = c(log(1 / 1.1), 0)
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = pessimistic_mean[[1,1]],
      sd = pessimistic_sd[[1,1]]),
    fill = "#ef8a47", alpha = 0.5,
    xlim = c(0, log(1.1))
  ) +
  geom_area(
    stat = "function", fun = dnorm,
    args = list(
      mean = pessimistic_mean[[1,1]],
      sd = pessimistic_sd[[1,1]]),
    fill = "#ef8a47",
    xlim = c(log(1.1),0.7)
  ) +
  geom_vline(
    xintercept = 0
  ) +
  # annotate(
  # "label", x = 0, y = 0.1, 
  # label = "ROPE") +
  annotate(
    "text", x = 0.4, y = 1.8, 
    label = glue::glue("Pr(harm) = {round(pessimistic_post_prob[[2, 3]], 2)}"),
    size = 3,
    colour = "#ef8a47") + 
  annotate(
    "text", x = -0.45, y = 1.8, 
    label = glue::glue("Pr(benefit) = {round(pessimistic_post_prob[[1, 3]], 2)}"),
    size = 3,
    colour = "#528fad") +
  scale_y_continuous(
    limits = c(0, 1.95),
    expand = c(0, 0)) +
  theme_classic() +
  labs(
    x = "log(HR)",
    y = "Density",
    subtitle = "Pessimistic prior: mean = 0.916, sd = 1.76"
  )

jpeg("poterior_plot.jpeg", width = 4, height = 10, units = "in", res = 300)
neutral_plot / optimistic_plot / pessimistic_plot +
  plot_annotation(tag_levels = 'A')
dev.off()
