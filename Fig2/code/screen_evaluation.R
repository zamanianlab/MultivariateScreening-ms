# data wrangling/plotting
library(tidyverse)
library(janitor)

# stats
library(broom)
library(tidymodels)
library(e1071)

# other plotting
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(wesanderson)
library(ZamanianLabThemes)
library(gt)
# library(ggrepel)

# misc
library(conflicted)
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

rescaled_data <- readRDS(here("Fig2/data/rescaled_data.rds"))

control <- rescaled_data %>%
  filter(treatment %in% c("Negative Control", "Positive Control"))


# z-factor ----------------------------------------------------------------

control_summary <- control %>%
  select(assay_date, treatment, motility, viability) %>%
  pivot_longer(cols = c(motility, viability), values_to = "value", names_to = "assay_type") %>%
  group_by(assay_date, treatment, assay_type) %>%
  summarise(mean = mean(value), sd = sd(value), n = n())

z_factor <- control_summary %>%
  pivot_wider(names_from = treatment, values_from = c(mean, sd, n)) %>%
  ungroup() %>%
  mutate(
    z_factor = 1 - (((3 * `sd_Negative Control`) + (3 * `sd_Positive Control`)) /
      abs(`mean_Negative Control` - `mean_Positive Control`)),
    x = 1.25,
    y = 40
  )

(z_factor_p <- control %>%
    select(assay_date, treatment, motility, viability) %>%
    pivot_longer(cols = c(motility, viability), values_to = "value", names_to = "assay_type") %>%
    mutate(assay_type = str_to_sentence(assay_type)) %>%
    ggplot() +
    geom_rect(
      data = control_summary %>% mutate(assay_type = str_to_sentence(assay_type)),
      aes(
        xmin = mean - (3 * sd), xmax = mean + (3 * sd),
        ymin = -Inf, ymax = Inf, fill = treatment
      ), alpha = 0.25
    ) +
    geom_histogram(aes(x = value, fill = treatment), alpha = 0.75, bins = 50) +
    geom_vline(
      data = control_summary %>% mutate(assay_type = str_to_sentence(assay_type)),
      aes(xintercept = mean, color = treatment),
      linetype = "dotted", size = 1
    ) +
    geom_text(
      data = z_factor %>% mutate(assay_type = str_to_sentence(assay_type)),
      aes(x = x, y = y, label = str_c("Z' = ", round(z_factor, digits = 2))),
      fontface = "italic", color = "black"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = wes_palette("BottleRocket2")) +
    scale_color_manual(values = wes_palette("BottleRocket2")) +
    facet_grid(rows = vars(assay_type), scales = "free_x") +
    labs(x = "Scaled phenotypic value", y = "Count", fill = "", color = "") +
    theme_nw2() +
    theme(axis.text.x = element_markdown(angle = 0, hjust = 0.5)) +
    NULL)

save_plot(here("Fig2/subplots/Fig2e.pdf"), z_factor_p, base_height = 4, base_width = 6)
saveRDS(z_factor_p, here("Fig2/subplots/primary_screen_z_factor.rds"))
