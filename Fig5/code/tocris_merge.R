# data wrangling/plotting
library(tidyverse)

# other plotting
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(ZamanianLabThemes)
library(ggdendro)
library(patchwork)

# misc
library(conflicted)
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")

# import data -------------------------------------------------------------

### mf data

numbers <- read_rds(here("Fig2/data/z_scores.rds")) %>%
  select(number, treatment)

primary_hits <- read_rds(here("Fig2/data/z_scores.rds")) %>%
  pivot_longer(contains("z"), names_to = "assay", values_to = "value") %>%
  mutate(
    stages = "Mf",
    assay = str_remove(assay, "z_")
  ) %>%
  select(-number)

dr <- read_rds(here("Fig2/data/ec50.rds")) %>%
  select(treatment, value = motility_36hr) %>%
  mutate(
    assay = "dose-response",
    stages = "Mf"
  )

### adult data
motility <- read_rds(here("Fig4/data/motility_end.rds")) %>%
  rename(value = total_motility) %>%
  mutate(assay = "motility")

fecundity <- read_rds(here("Fig4/data/proportions.rds")) %>%
  select(treatment, value = prop) %>%
  mutate(
    stages = "Adult female",
    assay = "fecundity"
  )

lactate <- read_rds(here("Fig4/data/lactate_data.rds")) %>%
  ungroup() %>%
  select(treatment, stages, value = rfu_norm) %>%
  mutate(assay = "metabolism")

# join it together!
all_data <- bind_rows(motility, fecundity, lactate, primary_hits, dr) %>%
  left_join(numbers) %>%
  mutate(
    assay = str_c(stages, assay, sep = " ") %>% str_remove(., "Adult ") %>% str_to_sentence(),
    assay = factor(assay, levels = c("Mf viability", "Mf motility", "Mf dose-response", "Female motility", "Female fecundity", "Female metabolism", "Male motility", "Male metabolism"))
  ) %>%
  group_by(assay, stages) %>%
  mutate(rescaled_value = scales::rescale(value, to = c(0, 1))) %>%
  arrange(-rescaled_value) %>%
  select(number, treatment, stages, assay, value, rescaled_value) %>%
  mutate(number = case_when(
    treatment == "DMSO" ~ "DMSO",
    treatment == "Ivermectin" ~ "Ivermectin",
    is.numeric(number) == TRUE ~ as.character(number),
    TRUE ~ NA_character_
  )) %>%
  drop_na(number)

means <- all_data %>%
  group_by(number, treatment, stages, assay) %>%
  summarise(
    mean_rescaled = mean(rescaled_value, na.rm = TRUE),
    mean_value = mean(value, na.rm = TRUE)
  ) %>%
  group_by(stages, assay) %>%
  arrange(-mean_rescaled) %>%
  mutate(order = row_number())

# heatmap -----------------------------------------------------------------

wide <- means %>%
  ungroup() %>%
  select(number, treatment, assay, mean_value) %>%
  # filter(assay != 'Mf dose-response') %>%
  pivot_wider(names_from = assay, values_from = mean_value) %>%
  replace_na(list("Mf dose-response" = 2)) %>%
  # mutate(
  #   `Mf viability` = case_when(
  #     treatment == 'DMSO' ~ 0,
  #     treatment == 'Ivermectin' ~ max(.$`Mf viability`, na.rm = TRUE),
  #     TRUE ~ `Mf viability`),
  #   `Mf motility` = case_when(
  #     treatment == 'DMSO' ~ 0,
  #     treatment == 'Ivermectin' ~ max(.$`Mf motility`, na.rm = TRUE),
  #     TRUE ~ `Mf motility`)) %>%
  drop_na(`Female motility`) %>%
  # remove controls to filter in an unbiased way
  # filter(!treatment %in% c('DMSO', 'Ivermectin')) %>%
  mutate(
    # across(where(is.numeric), ~scale(log2(.x + 1), scale = T, center = T)),
    across(where(is.numeric), scales::rescale),
    treatment = str_replace(treatment, "Ivermectin", "Positive Control")
  )

wide %>%
  pivot_longer(where(is.numeric), names_to = "phenotype", values_to = "val") %>%
  ggplot(aes(x = phenotype, y = val)) +
  geom_quasirandom() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  NULL

dists <- wide %>%
  select(where(is.numeric), treatment) %>%
  column_to_rownames(var = "treatment") %>%
  dist(method = "euclidean")

clust <- hclust(dists, method = "ward.D2")

ord <- clust$order

hm <- as.data.frame(as.matrix(dists)) %>%
  rownames_to_column(var = "treatment") %>%
  left_join(numbers)

hm$treatment <- factor(hm$treatment,
                       levels = c(hm$treatment[ord])
)

hm <- hm %>%
  pivot_longer(c(-treatment, -number), names_to = "treatment2", values_to = "dist")

hm$treatment2 <- factor(hm$treatment2,
                        levels = c(hm$treatment2[ord])
)

(control_dists <- hm %>%
    mutate(
      log_dist = log10(dist),
      rescale = scales::rescale(dist, to = c(-1, 1))
    ) %>%
    filter(treatment %in% c("DMSO", "Positive Control")) %>%
    ggplot() +
    geom_density(aes(x = dist, color = treatment)) +
    labs(x = "euclidean distance"))

labels <- hm %>%
  distinct(treatment, number) %>%
  mutate(number = case_when(
    treatment == "DMSO" ~ "-",
    treatment == "Positive Control" ~ "+",
    TRUE ~ as.character(number)
  )) %>%
  tibble::deframe()

cluster_labels <- tribble(
  ~x, ~y, ~label,
  14, 10, "**↑ Mf EC<sub>50</sub><br>↓ Adult effects**",
  10, 19, "**↓ Mf EC<sub>50</sub><br>↓ Adult effects**",
  23, 20.5, "**↑ Mf EC<sub>50</sub><br>↑ Adult effects**",
  19, 30.5, "**↓ Mf EC<sub>50</sub><br>↑ Adult effects**"
)

# heatmap showing euclidean distances between treatments
(dist_map <- hm %>%
    mutate(
      dist = case_when(
        dist == 0 ~ NA_real_,
        TRUE ~ dist
      )
    ) %>%
    mutate(
      log_dist = log10(dist),
      rescale = scales::rescale(dist, to = c(-1, 1))
    ) %>%
    ggplot() +
    geom_tile(aes(x = treatment, y = treatment2, fill = rescale)) +
    # highlight positive control
    annotate(
      geom = "rect",
      xmin = 0.5, ymin = 19.5,
      xmax = 32.5, ymax = 18.5,
      color = NA, fill = "#CB2314", alpha = 0.5
    ) +
    annotate(
      geom = "rect",
      xmin = 18.5, ymin = 0.5,
      xmax = 19.5, ymax = 32.5,
      color = NA, fill = "#CB2314", alpha = 0.5
    ) +
    # highlight DMSO
    annotate(
      geom = "rect",
      xmin = 0.5, ymin = 5.5,
      xmax = 32.5, ymax = 6.5,
      color = NA, fill = "#FAD510", alpha = 0.5
    ) +
    annotate(
      geom = "rect",
      xmin = 5.5, ymin = 0.5,
      xmax = 6.5, ymax = 32.5,
      color = NA, fill = "#FAD510", alpha = 0.5
    ) +
    # boxes around hits
    annotate(
      geom = "rect",
      xmin = 14.6, ymin = 14.6,
      xmax = 18.4, ymax = 18.4,
      color = "#00AFB5", fill = NA, linetype = "solid", size = 1
    ) +
    annotate(
      geom = "rect",
      xmin = 18.6, ymin = 18.6,
      xmax = 23.4, ymax = 23.4,
      color = "#00AFB5", fill = NA, linetype = "solid", size = 1
    ) +
    annotate(
      geom = "rect",
      xmin = 23.6, ymin = 23.6,
      xmax = 32.4, ymax = 32.4,
      color = "#00AFB5", fill = NA, linetype = "solid", size = 1
    ) +
    # box around DMSO-like
    annotate(
      geom = "rect",
      xmin = 0.625, ymin = 0.625,
      xmax = 14.4, ymax = 14.4,
      color = "black", fill = NA, alpha = 0.5, linetype = "solid", size = 1
    ) +
    geom_richtext(
      data = cluster_labels,
      aes(x = x, y = y, label = label),
      color = "white", fill = alpha("grey30", 0.75),
      hjust = 0, size = 8 / 14 * 5
    ) +
    scale_x_discrete(position = "top", labels = labels) +
    scale_fill_gradient2(
      low = "#ffbf46", high = "grey20", mid = "#575761",
      na.value = "black"
    ) +
    coord_fixed() +
    theme_void() +
    theme(
      axis.text.x = element_text(size = 8, margin = margin(b = 3)),
      axis.text.y = element_text(size = 8, hjust = 1, margin = margin(r = 3)),
      axis.ticks = element_line(size = 0.25, linetype = "solid"),
      axis.ticks.length = unit(3, "pt"),
      plot.margin = unit(c(7, 0, 2, 7), "pt"),
      legend.position = "empty"
    ) +
    NULL)

(pheno_map <- wide %>%
    mutate(across(.cols = c("Mf viability", "Mf motility"), ~ scales::rescale(.x, to = c(1, 0)))) %>%
    pivot_longer(where(is.numeric), names_to = "phenotype", values_to = "value") %>%
    ggplot() +
    geom_tile(aes(x = treatment, y = phenotype, fill = value)) +
    annotate(
      geom = "rect",
      xmin = 5.5, ymin = 0.5,
      xmax = 6.5, ymax = 8.5,
      color = NA, fill = "#FAD510", alpha = 0.5
    ) +
    annotate(
      geom = "rect",
      xmin = 18.5, ymin = 0.5,
      xmax = 19.5, ymax = 8.5,
      color = NA, fill = "#CB2314", alpha = 0.5
    ) +
    geom_vline(
      data = tibble(xintercept = c(14.63, 18.5, 23.5)),
      aes(xintercept = xintercept),
      color = "#00AFB5", linetype = "solid", size = 1
    ) +
    geom_vline(
      xintercept = 14.37,
      color = "#273046", linetype = "solid", size = 1
    ) +
    scale_fill_gradient(
      low = "#ffbf46", high = "grey20",
      na.value = "black"
    ) +
    coord_fixed() +
    theme_void() +
    theme(
      axis.text.y = element_markdown(size = 8, hjust = 1, margin = margin(r = 3)),
      axis.ticks.y = element_line(size = 0.25, linetype = "solid"),
      axis.ticks.length = unit(3, "pt"),
      plot.margin = unit(c(0, 7, 7, 7), "pt"),
      legend.position = "empty"
    ) +
    NULL)

(dendro <- ggdendro::ggdendrogram(clust, labels = FALSE) +
    annotate("segment",
             x = 1, xend = 14,
             y = -0.1, yend = -0.1,
             colour = "black", size = 2
    ) +
    annotate("segment",
             x = 15, xend = 18,
             y = -0.1, yend = -0.1,
             colour = "#00AFB5", size = 2
    ) +
    annotate("segment",
             x = 19, xend = 23,
             y = -0.1, yend = -0.1,
             colour = "#00AFB5", size = 2
    ) +
    annotate("segment",
             x = 24, xend = 32,
             y = -0.1, yend = -0.1,
             colour = "#00AFB5", size = 2
    ) +
    scale_x_continuous(limits = c(1, 32), expand = c(0.05, 0.05)) +
    theme_void() +
    NULL)

heatmaps <- dist_map / pheno_map

join <- ggdraw(heatmaps, xlim = c(0, 1), ylim = c(0, 1.3)) +
  draw_plot(dendro, 0.182, 0.98, 0.8, 0.3)

save_plot(here("Fig5/Fig5.pdf"), join, base_height = 8.8, base_width = 6.5)
ggsave(here("Fig5/Fig5.png"), join, height = 8.8, width = 6.5, units = "in", bg = "white")
