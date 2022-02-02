library(tidyverse)
library(survival)
library(cmprsk)


## Data load
predesci <- 
  read_csv("predesci_data.csv") %>%
  mutate(
    binary = if_else(outcome == 1, 1, 0),
    group_num = if_else(group == "bb", 0, 1)
  )

## Survival comparisons
### Cause specific
cox_fit <- coxph(Surv(time, binary) ~ group_num, data = predesci)
cox_pval <- round(summary(cox_fit)$coefficients[1,5],3)

### Competing risk
cum_inc_out <- cuminc(predesci$time, predesci$outcome, predesci$group, cencode = 0)
cmprsk_pval <-round(cum_inc_out$Tests[[1,2]], 3)

## Plot cumulative incidence decompensation / liver-related mortality
### Cumulative incidence curve data prep
cum_inc_data <-
  cum_inc_out %>%
  list_modify("Tests" = NULL) %>%
  map_df(`[`, c("time", "est", "var"), .id = "id") %>%
  mutate(id = recode(
    id,
    "bb 1" = "bb:Decomp/dead",
    "pbo 1" = "pbo:Decomp/dead",
    "bb 2" = "bb:Transplant",
    "pbo 2" = "pbo:Transplant",
    "bb 3" = "bb:Non-liver",
    "pbo 3" = "pbo:Non-liver")
  ) %>%
  separate(id, c("Group", "Cause"), ":")

### Plot to match PREDESCI manuscript figure 2
predesci_plot <-
  ggplot(
    cum_inc_data %>% filter(Cause == "Decomp/dead")) +
  geom_step(
    aes(x = time, y = est, colour = Group)
  ) +
  scale_colour_manual(values = c("#1e466e", "#e76254")) +
  scale_x_continuous(
    breaks = c(12, 24, 36, 48, 60),
    limits = c(0, 60),
    name = "Time (months)"
  ) +
  ylab("Cumulative incidence\n(decompensation  or liver-related death)") +
  ylim(0, 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  # annotate("text", x = 8, y = 0.40, label = glue::glue("p = {p_val}"))  +
  annotate("text", x = 55, y = 0.44, label = "Placebo", colour = "#e76254") +
  annotate("text", x = 55, y = 0.24, label = expression(beta ~ " blocker"), colour = "#1e466e") +
  annotate("text", x = 0, y = 0.48, 
           label = glue::glue("Cox cause specific p value = {cox_pval}\nFine & Gray p value = {cmprsk_pval}"),
           hjust = 0)


jpeg("PREDESCI_data_plot.jpeg", width = 4, height = 4, units = "in", res = 300)
predesci_plot
dev.off()
