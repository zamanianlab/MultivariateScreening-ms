# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(wesanderson)
library(ZamanianLabThemes)

# misc
library(conflicted)
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")


# analysis ----------------------------------------------------------------

primary_hits <-
  read_csv(here("Fig2/tables/Table1.csv")) %>%
  select(number, treatment) %>% 
  mutate(number = as.character(number))

tidy_data <-
  read_rds(here("Fig4/data/lactate_data.rds")) %>%
  left_join(primary_hits) %>%
  mutate(
    assay_date = str_extract(plate, "^20[0-9]{6}"),
    number = case_when(
      treatment %in% c("Ivermectin", "DMSO") ~ treatment,
      !is.na(number) ~ as.character(number),
      is.na(number) ~ NA_character_
    )
  ) %>%
  drop_na(number) %>%
  group_by(stages, treatment) %>%
  mutate(worm_number = row_number())

control <- tidy_data %>%
  filter(treatment == "DMSO") %>%
  group_by(assay_date, stages) %>%
  summarise(control_mean = mean(rfu))

# remove the worms that weren't moving at time_point 0
dead <- tribble(
  ~treatment, ~stages, ~well,
  "Ivermectin", "Adult male", "F02",
  "PJ 34 hydrochloride", "Adult male", "E06",
  "PJ 34 hydrochloride", "Adult male", "G06",
  "SSR 146977 hydrochloride", "Adult male", "E07",
  "SSR 146977 hydrochloride", "Adult male", "H07",
  "Sildenafil citrate", "Adult male", "D01",
  "GSK J4", "Adult male", "A06",
  "N106", "Adult male", "F01",
  "ML 418", "Adult male", "G12",
  "Ciclopirox", "Adult male", "G09",
  "Sal 003", "Adult female", "B11",
  "Sal 003", "Adult female", "C11"
) %>%
  mutate(dead = "dead")

lactate_data <- tidy_data %>%
  left_join(dead) %>%
  filter(is.na(dead)) %>%
  select(-dead) %>%
  filter(str_detect(plate, "20211017")) %>%
  left_join(control) %>%
  mutate(rfu_norm = rfu / control_mean) %>%
  ungroup() %>%
  # hackish way of reordering by group
  mutate(treatment = str_remove(treatment, " hydrochloride| maleate| citrate")) %>%
  mutate(reorder = str_c(treatment, stages, sep = "_")) %>%
  mutate(reorder = fct_reorder(reorder, rfu_norm, mean, .desc = TRUE)) %>%
  mutate(
    reorder = fct_relevel(reorder, "DMSO_Adult female", "DMSO_Adult male", "Ivermectin_Adult female", "Ivermectin_Adult male")
  )

metabolism_hits <- lactate_data %>%
  mutate(individual_hit = case_when(
    rfu_norm < 0.5 ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  ungroup() %>%
  group_by(stages, treatment) %>%
  summarise(
    treatment_mean = mean(rfu_norm),
    n_hit = sum(individual_hit),
    n = n()
  ) %>%
  mutate(
    hit = case_when(
      treatment %in% c("Ivermectin", "DMSO") ~ treatment,
      (n_hit > 2) ~ "Hit",
      treatment_mean < 0.75 ~ "Hit",
      TRUE ~ "Non-Hit"
    )
  )

(lactate_plot <- lactate_data %>%
  left_join(metabolism_hits) %>%
  ggplot(aes(x = reorder, y = rfu_norm, group = stages, shape = stages, color = hit, fill = hit)) +
  geom_errorbar(
    position = position_dodge(width = 1),
    stat = "summary", width = 0.5, alpha = 0.75
  ) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean") +
  geom_quasirandom(
    shape = 21, fill = "white", size = 1,
    alpha = 0.75, dodge.width = 0.9, width = 0.1,
    show.legend = FALSE
  ) +
  scale_color_manual(
    limits = c("Non-Hit", "DMSO", "Ivermectin", "Hit"),
    values = c("#273046", "#FAD510", "#CB2314", "#00AFB5")
  ) +
  scale_fill_manual(
    limits = c("Non-Hit", "DMSO", "Ivermectin", "Hit"),
    values = c("#273046", "#FAD510", "#CB2314", "#00AFB5")
  ) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = ~ str_remove(.x, "_.*")) +
  facet_wrap(facets = vars(stages), nrow = 2, scales = "free") +
  labs(x = "Treatment", y = "Normalized RFU<br>(Percent control)", color = "", fill = "") +
  theme_nw2() +
  theme(
    panel.grid.major.x = element_blank(),
    strip.text.x = element_markdown(face = "italic")
  ) +
  NULL)

save_plot(here('Fig4/subplots/Fig4d.pdf'), lactate_plot, base_width = 7)
write_rds(lactate_plot, here("Fig4/subplots/lactate.rds"))

write_rds(metabolism_hits, here("Fig4/data/metabolism_hits.rds"))
