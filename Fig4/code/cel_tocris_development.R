# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(ZamanianLabThemes)
library(wesanderson)
library(ggrepel)

# stats
library(broom)
library(tidymodels)

# misc
library(here)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# snake order of the IX recording
record_order <- tibble(
  row = LETTERS[1:8], `01` = seq(1, 8, 1), `02` = seq(16, 9, -1),
  `03` = `01` + 16, `04` = `02` + 16, `05` = `03` + 16, `06` = `04` + 16,
  `07` = `05` + 16, `08` = `06` + 16, `09` = `07` + 16, `10` = `08` + 16,
  `11` = `09` + 16, `12` = `10` + 16
) %>%
  pivot_longer(`01`:`12`, names_to = "col", values_to = "well_number")

pruned_data <- read_rds(here("Fig4/data/cel_tocris_development_data.rds"))

# check for drift
(drift <- pruned_data %>%
  left_join(record_order) %>%
  select(plate:conc, well_number, type, area = area_shape_area, major_axis_length = area_shape_major_axis_length) %>%
  arrange(well_number) %>%
  mutate(smooth_row = as.numeric(factor(row))) %>%
  filter(!col %in% c("01", "12")) %>%
  ggplot(aes(x = smooth_row, y = area)) +
  geom_point(size = 1, alpha = 0.75) +
  geom_smooth(method = "lm", color = "indianred") +
  scale_x_continuous(labels = c("A", "C", "E", "G"), breaks = seq(1, 8, 2)) +
  facet_grid(cols = vars(col)) +
  labs(x = "Row", y = "Raw motility") +
  theme_nw2() +
  theme(legend.position = "empty") +
  NULL)

# z-score -----------------------------------------------------------------

well_summary <- pruned_data %>%
  select(plate:conc, type, area = area_shape_area, major_axis_length = area_shape_major_axis_length) %>%
  group_by(across(!where(is.numeric))) %>%
  summarise(mean_major_axis_length = mean(major_axis_length, na.rm = TRUE))

plate_summary <- well_summary %>%
  ungroup() %>%
  group_by(plate) %>%
  summarise(
    mean_plate_major_axis_length = mean(mean_major_axis_length, na.rm = TRUE),
    sd_plate_major_axis_length = sd(mean_major_axis_length, na.rm = TRUE)
  )

z <- well_summary %>%
  left_join(plate_summary) %>%
  mutate(z_major_axis_length = (mean_major_axis_length - mean_plate_major_axis_length) / sd_plate_major_axis_length)

primary_hits <- read_csv(here("Fig2/tables/Table1.csv")) %>%
  select(number, treatment)

cutoff <- -1

(z_plot <- z %>%
  filter(!plate %in% c("p03", "p04")) %>%
  mutate(
    type = case_when(
      treatment %in% primary_hits$treatment ~ "Primary Hit",
      !treatment %in% c("Negative Control", "Positive Control") ~ "Treated",
      TRUE ~ treatment
    )
  ) %>%
  mutate(type = factor(type, levels = c("Primary Hit", "Negative Control", "Positive Control", "Treated"))) %>%
  ggplot(aes(x = species, y = z_major_axis_length, color = type)) +
  annotate("rect",
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = cutoff,
    alpha = .2, color = NA, fill = "#354823"
  ) +
  geom_beeswarm() +
  # label overlap
  geom_label_repel(
    data = . %>% filter(z_major_axis_length < cutoff, type == "Primary Hit"),
    aes(label = treatment, fill = type),
    size = 3,
    x = 1.2, direction = "y", point.size = NA, min.segment.length = 1000,
    color = "white", fontface = "italic", show.legend = FALSE, max.overlaps = 1000
  ) +
  # label non-overlap
  geom_label_repel(
    data = . %>% filter(z_major_axis_length < cutoff, type == "Treated"),
    aes(label = treatment, fill = type),
    size = 3,
    x = 0.8, direction = "y", point.size = NA, min.segment.length = 1000,
    color = "white", fontface = "italic", show.legend = FALSE
  ) +
  scale_color_manual(values = c("#00AFB5", "#FAD510", "#CB2314", "#273046")) +
  scale_fill_manual(values = c("#00AFB5", "#273046")) +
  labs(x = "", y = "Development Z-score", color = "", fill = "") +
  theme_nw2() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  NULL)

save_plot(here("Fig4/subplots/Fig4a.pdf"), z_plot, base_width = 6, base_height = 6)
write_rds(z_plot, here("Fig4/subplots/cel_screen_dev.rds"))

one_hits <- z %>%
  ungroup() %>%
  filter(z_major_axis_length < cutoff, treatment != "Positive Control", !plate %in% c("p03", "p04")) %>%
  select(treatment, conc)


# 10 µM follow up ---------------------------------------------------------

followup <- pruned_data %>%
  filter(plate %in% c("p03", "p04")) %>%
  mutate(
    type = case_when(
      str_detect(type, "Control") ~ type,
      TRUE ~ "Primary Hit"
    ),
    treatment = str_remove(treatment, " hydrochloride| dihydrochloride| maleate| citrate"),
    treatment = fct_reorder(treatment, area_shape_major_axis_length, median, .desc = TRUE),
    treatment = fct_relevel(treatment, "Negative Control", "Positive Control", after = Inf)
  )

followup_stat_layer <- followup %>%
  mutate(treatment = fct_relevel(treatment, "Negative Control")) %>%
  group_nest(strains) %>%
  mutate(otm = map(data, ~ PMCMRplus::kwManyOneDunnTest(.x$area_shape_major_axis_length ~ as.factor(.x$treatment), p.adjust.method = "bonferroni"))) %>%
  mutate(summary = map(otm, ~ summary(.x))) %>%
  mutate(pval = map(summary, ~ as_tibble(pluck(summary, 1, 3), rownames = "treatment"))) %>%
  mutate(estimate = map(summary, ~ as_tibble(pluck(summary, 1, 4), rownames = "treatment"))) %>%
  mutate(tidy = map2(pval, estimate, ~ left_join(.x, .y, by = "treatment"))) %>%
  unnest(tidy) %>%
  select(-data, -otm, -summary, -pval, -estimate, pval = `Negative Control.x`, estimate = `Negative Control.y`) %>%
  mutate(sig = case_when(
    pval <= 0.0001 ~ "****",
    pval <= 0.001 ~ "***",
    pval <= 0.01 ~ "**",
    pval <= 0.05 ~ "*",
    pval > 0.05 ~ ""
  )) %>%
  left_join(select(followup, treatment, type) %>% distinct(treatment, type)) %>%
  select(treatment, type, pval, sig) %>%
  mutate(
    area_shape_major_axis_length = 400
  ) %>%
  filter(sig != "")

(followup_plot <- followup %>%
    ggplot(aes(x = treatment, y = area_shape_major_axis_length)) +
    ggdist::stat_halfeye(aes(fill = type),
                         point_size = 1, interval_alpha = 0.5) +
    geom_text(
      data = followup_stat_layer,
      aes(label = sig, color = type),
      show.legend = FALSE, fontface = "plain"
    ) +
    # scale_y_continuous(limits = c(0, 7000)) +
    scale_fill_manual(values = c("#FAD510", "#CB2314", "#00AFB5")) +
    scale_color_manual(values = c("#CB2314", "#00AFB5")) +
    labs(y = "Worm length (pixels)") +
    coord_flip() +
    theme_nw2() +
    theme(
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      legend.position = "empty"
    ) +
    NULL)

save_plot(here("Fig4/subplots/Fig4b.pdf"), followup_plot, base_width = 6, base_height = 8)
write_rds(followup_plot, here("Fig4/subplots/cel_followup.rds"))

ten_hits <- followup_stat_layer %>%
  filter(treatment != "Positive Control") %>%
  select(treatment) %>%
  mutate(conc = "10uM")

cel_hits <- bind_rows(one_hits, ten_hits) %>%
  write_rds(here("Fig4/data/cel_hits.rds"))
