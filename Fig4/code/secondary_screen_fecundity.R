# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(ZamanianLabThemes)

# misc
library(conflicted)
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")


# import and fit data -----------------------------------------------------

pixel_data <- read_rds(here("Fig4/data/raw_fecundity_data.rds"))

model <- read_rds(here("Fig4/data/fecundity_model.rds"))

fitted_data <- pixel_data %>%
  group_by(plate) %>%
  group_nest() %>%
  mutate(predic = map(data, ~ broom::augment(model, new_data = .x))) %>%
  unnest(cols = predic) %>%
  select(-data)

rescaled_data <- ungroup(fitted_data) %>%
  group_by(assay_date) %>%
  mutate(scaled = scales::rescale(.pred, to = c(0, max(.pred)))) %>%
  ungroup()


# analysis ----------------------------------------------------------------

primary_hits <- read_csv(here("Fig2/tables/Table1.csv")) %>%
  select(number, treatment)

day_data <- rescaled_data %>%
  left_join(primary_hits) %>%
  mutate(number = case_when(
    treatment %in% c("Ivermectin", "DMSO") ~ treatment,
    !is.na(number) ~ as.character(number),
    is.na(number) ~ NA_character_
  )) %>%
  drop_na(number) %>%
  group_by(assay_date, time_point, species, number, treatment, conc) %>%
  mutate(worm_number = row_number()) %>%
  mutate(sex = "Female") %>%
  ungroup() %>%
  mutate(time_point = factor(time_point, levels = c("0hr", "48hr"), ordered = TRUE)) %>%
  select(-plate) %>%
  pivot_wider(id_cols = c(assay_date, species, treatment, conc, worm_number), names_from = "time_point", values_from = "scaled") %>%
  mutate(treatment = str_remove(treatment, " hydrochloride| maleate| citrate"))

proportions <- day_data %>%
  filter(
    !is.na(`0hr`),
    `0hr` > 200
  ) %>%
  mutate(
    prop = `48hr` / `0hr`,
    conc = factor(conc,
      levels = c("0.1p", "1uM"),
      labels = c("0.1%", "1 µM")
    )
  ) %>%
  select(-`0hr`, -`48hr`)

fecundity_hits <- proportions %>%
  mutate(individual_hit = case_when(
    prop < 0.5 ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  ungroup() %>%
  group_by(treatment) %>%
  summarise(
    treatment_mean = mean(prop),
    n_hit = sum(individual_hit),
    n = n()
  ) %>%
  mutate(hit = case_when(
    treatment %in% c("Ivermectin", "DMSO") ~ treatment,
    (n_hit > 2) ~ "Hit",
    treatment_mean < 0.5 ~ "Hit",
    TRUE ~ "Non-Hit"
  )) %>%
  # remove treatments that had too few worms because of the 200 cutoff
  filter(n > 1)

# show the distribution of fecundity pre-treatment
n <- length(day_data$`0hr`)

(baseline <- day_data %>%
  left_join(fecundity_hits, by = "treatment") %>%
  # mutate(hit_type = case_when(
  #   treatment == 'DMSO' ~ 'DMSO',
  #   is.na(hit_type) ~ 'Non-hit',
  #   TRUE ~ hit_type
  # )) %>%
  # mutate(hit_type = factor(hit_type, levels = c('DMSO', 'Non-hit', 'Ivermectin', 'Male', 'Female', 'Both'))) %>%
  ggplot(aes(x = species, y = `0hr`, color = hit)) +
  geom_beeswarm(cex = 3.5, alpha = 0.5) +
  geom_hline(yintercept = 200, linetype = "dashed", color = "black") +
  annotate("text", label = str_glue("N = {n}"), x = 1, y = -500, size = 3) +
  annotate("text", label = "200", x = 1.24, y = 500, size = 2, fontface = "italic") +
  scale_color_manual(
    limits = c("Non-Hit", "DMSO", "Ivermectin", "Hit"),
    values = c("#273046", "#FAD510", "#CB2314", "#00AFB5")
  ) +
  labs(x = "", y = "Pre-treatment") +
  coord_cartesian(
    ylim = c(0, 3000),
    clip = "off"
  ) +
  theme_nw2() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "empty"
  ) +
  NULL)

write_rds(baseline, here("Fig4/subplots/baseline.rds"))
save_plot(here("Fig4/subplots/Fig4c_inset.pdf"), baseline, base_width = 3)

# plot the percent 0hr proportions
(proportions_plot <- proportions %>%
  left_join(fecundity_hits, by = "treatment") %>%
  drop_na() %>%
  mutate(
    # reorder so the bar chart is a waterfall plot
    treatment = fct_reorder(treatment, prop, mean, .desc = TRUE),
    hit = factor(hit, levels = c("DMSO", "Non-Hit", "Ivermectin", "Hit")),
    treatment = fct_relevel(treatment, "DMSO", "Ivermectin"),
    stages = "Adult female"
  ) %>%
  ggplot(aes(x = treatment, y = prop, color = hit, fill = hit)) +
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
  facet_grid(cols = vars(stages), scales = "free_x") +
  labs(x = "", y = "Inferred progeny, 48 hr<br>(percent 0 hr)", fill = "", color = "") +
  theme_nw2() +
  theme(
    panel.grid.major.x = element_blank(),
    strip.text.x = element_markdown(face = "italic")
  ) +
  guides(col = guide_legend(ncol = 6)) +
  NULL
)

save_plot(here("Fig4/subplots/Fig4c.pdf"), proportions_plot, base_width = 7)
write_rds(proportions_plot, here("Fig4/subplots/proportions.rds"))

# export data for merge ---------------------------------------------------

prop_summary <- proportions %>%
  group_by(treatment) %>%
  summarize(mean_prop = mean(prop))

write_rds(prop_summary, here("Fig4/data/prop_summary.rds"))
write_rds(proportions, here("Fig4/data/proportions.rds"))

write_rds(fecundity_hits, here("Fig4/data/fecundity_hits.rds"))
