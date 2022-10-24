library(tidyverse)
library(patchwork)

# Plot priors
prior_neutral_plot <-
  ggplot(data.frame(x = c(-4, 4)), aes(x)) +
  stat_function(
    fun = dnorm, 
    n = 100000, 
    args = list(
      mean = 0,
      sd = 0.355),
    linetype = 1
  ) +
  geom_vline(
    xintercept = 0
  ) +
  scale_y_continuous(
    limits = c(0, 1.3),
    expand = c(0, 0)) +
  theme_classic() +
  labs(
    x = "log(HR)",
    y = "Density",
    subtitle = "Neutral ('skeptical') prior: mean = 0, sd = 0.355"
  )

prior_optimistic_plot <-
  ggplot(data.frame(x = c(-4, 4)), aes(x)) +
  stat_function(
    fun = dnorm, 
    n = 100000, 
    args = list(
      mean = -0.916,
      sd = 0.88),
    linetype = 1
  ) +
  geom_vline(
    xintercept = 0
  ) +
  scale_y_continuous(
    limits = c(0, 1.3),
    expand = c(0, 0)) +
  theme_classic() +
  labs(
    x = "log(HR)",
    y = "Density",
    subtitle = "Optimistic prior: mean = -0.916, sd = 0.88"
  )

prior_pessimistic_plot <-
  ggplot(data.frame(x = c(-4, 4)), aes(x)) +
  stat_function(
    fun = dnorm, 
    n = 100000, 
    args = list(
      mean = 0.916,
      sd = 1.76),
    linetype = 1
  ) +
  geom_vline(
    xintercept = 0
  ) +
  scale_y_continuous(
    limits = c(0, 1.3),
    expand = c(0, 0)) +
  theme_classic() +
  labs(
    x = "log(HR)",
    y = "Density",
    subtitle = "Pessimistic prior: mean = 0.916, sd = 1.76"
  )

jpeg("prior_plot.jpeg", width = 4, height = 10, units = "in", res = 300)
prior_neutral_plot / prior_optimistic_plot / prior_pessimistic_plot +
  plot_annotation(tag_levels = 'A')
dev.off()
