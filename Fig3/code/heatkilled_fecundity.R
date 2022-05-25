# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(ggridges)
library(ZamanianLabThemes)
library(wesanderson)

# misc
library(conflicted)
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")


# import data and organize ------------------------------------------------

data <- read_rds(here("Fig3/data/fecundity_data.rds")) %>%
  separate(other, c("time_point", "exp_date"), "_") %>%
  # filter for the heat killed data
  filter(assay_date == "20210825") %>%
  arrange(plate, col, well) %>%
  group_by(time_point, species, treatment, conc, exp_date) %>%
  mutate(worm_number = row_number()) %>%
  ungroup() %>%
  mutate(time_point = factor(time_point, levels = c("0hr", "24hr", "48hr"), ordered = TRUE))

(hk_plot <- data %>%
  mutate(time_point = factor(time_point,
    levels = c("0hr", "24hr", "48hr"),
    labels = c("0 hr", "24 hr", "48 hr")
  )) %>%
  pivot_wider(id_cols = c(exp_date, treatment, worm_number), names_from = time_point, values_from = scaled) %>%
  filter(`0 hr` > 200) %>%
  mutate(
    prop = (`24 hr` + `48 hr`) / `0 hr`,
    treatment = case_when(
      treatment == "Alive" ~ "Alive",
      TRUE ~ "Heat-<br>killed"
    )
  ) %>%
  ggplot(aes(x = treatment, y = prop)) +
  geom_errorbar(
    position = position_dodge(width = 1),
    stat = "summary", width = 0.5, alpha = 0.75,
    color = "black"
  ) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean", fill = "black") +
  geom_quasirandom(
    shape = 21, fill = "white", size = 2,
    alpha = 0.75, dodge.width = 0.9, width = 0.1,
    show.legend = FALSE, color = "black"
  ) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "", y = "Inferred progeny<br>(percent 0 hr)", color = "Treatment") +
  theme_nw2() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_markdown(angle = 0, hjust = 0.5)
  ) +
  NULL)

save_plot(here("Fig3/subplots/Fig3e.pdf"), hk_plot, base_width = 2)
write_rds(hk_plot, here("Fig3/subplots/heatkilled_fecundity.rds"))
