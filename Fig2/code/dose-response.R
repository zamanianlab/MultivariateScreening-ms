# data wrangling/plotting
library(tidyverse)
library(janitor)

# stats
library(broom)
library(drc)

# other plotting
library(ggbeeswarm)
library(cowplot)
library(ggtext)
library(ZamanianLabThemes)
library(gt)
library(ggrepel)

library(here)

# misc
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# import ------------------------------------------------------------------

m_data <- read_rds(here("Fig2/data/motility_data.rds")) %>%
  mutate(assay_date = str_extract(plate, "20[0-9]{6}"), .before = plate) %>%
  mutate(
    norm_flow = optical_flow / log2(worm_area), .after = worm_area,
    plate = str_extract(plate, "p[0-9]{2}")
  ) %>%
  filter(assay_date %in% c("20210916", "20210917")) %>%
  # GNF 5837 stained the worms yellow and autofluoresced
  filter(treatment != "GNF 5837") %>%
  rename(time_point = other) %>%
  filter(case_when(!(assay_date == "20210917" & plate == "p08") ~ TRUE)) %>%
  rename(motility = norm_flow)

v_data <- read_rds(here("Fig2/data/celltox_data.rds")) %>%
  mutate(assay_date = as.character(as.numeric(str_extract(plate, "20[0-9]{6}"))), .before = plate) %>%
  filter(assay_date == "20210917") %>%
  mutate(plate = str_extract(plate, "p[0-9]{2}")) %>%
  mutate(plate = str_c("p", str_pad(as.numeric(str_extract(plate, "[0-9]{2}")) - 14, width = 2, side = "left", pad = 0))) %>%
  filter(plate != "p08") %>%
  # GNF 5837 stained the worms yellow and autofluoresced
  filter(treatment != "GNF 5837") %>%
  rename(time_point = other) %>%
  group_by(across(assay_date:time_point)) %>%
  summarize(viability = mean(area_occupied_area_occupied_green_worms)) %>%
  select(assay_date, plate:time_point, viability)

dr_data <- left_join(m_data, v_data) %>%
  select(-optical_flow, -worm_area) %>%
  mutate(
    treatment = case_when(
      treatment == "HeatKilled" ~ "Positive Control",
      treatment == "DMSO" ~ "Negative Control",
      TRUE ~ treatment
    ),
    type = case_when(
      str_detect(treatment, "Control") ~ treatment,
      TRUE ~ "Treated"
    )
  ) %>%
  pivot_longer(cols = c(viability, motility), names_to = "assay_type", values_to = "value") %>%
  filter(!is.na(value))


# trim outliers -----------------------------------------------------------

distribution <- dr_data %>%
  ggplot(aes(x = 1, y = value)) +
  geom_quasirandom(aes(color = type), alpha = 0.25, size = 3) +
  geom_boxplot(
    data = . %>% filter(str_detect(type, "Control")),
    aes(color = type),
    x = 0, position = position_dodge(0)
  ) +
  scale_color_manual(values = c("#FAD510", "#CB2314", "#273046")) +
  facet_wrap(facets = vars(assay_type, time_point), scales = "free_y") +
  theme_nw2() +
  NULL

coef <- 3

outliers <- dr_data %>%
  filter(str_detect(type, "Control")) %>%
  group_by(assay_type, time_point, type) %>%
  group_nest() %>%
  mutate(quant = map(data, ~ as.numeric(quantile(.x$value, c(0, 0.25, 0.5, 0.75, 1))))) %>%
  mutate(iqr = unlist(map(quant, ~ diff(.x[c(2, 4)])))) %>%
  mutate(cutoff = case_when(
    assay_type == "viability" & type == "Negative Control" ~ unlist(map(quant, ~ pluck(.x, 4))),
    assay_type == "viability" & type == "Positive Control" ~ unlist(map(quant, ~ pluck(.x, 2))),
    assay_type == "motility" & type == "Negative Control" ~ unlist(map(quant, ~ pluck(.x, 2))),
    assay_type == "motility" & type == "Positive Control" ~ unlist(map(quant, ~ pluck(.x, 4)))
  )) %>%
  mutate(cutoff = case_when(
    assay_type == "viability" & type == "Negative Control" ~ cutoff + coef * iqr,
    assay_type == "viability" & type == "Positive Control" ~ cutoff - coef * iqr,
    assay_type == "motility" & type == "Negative Control" ~ cutoff - coef * iqr,
    assay_type == "motility" & type == "Positive Control" ~ cutoff + coef * iqr,
  )) %>%
  select(assay_type:type, cutoff)

distribution <- dr_data %>%
  ggplot(aes(x = 1, y = value)) +
  geom_quasirandom(aes(color = type), alpha = 0.25, size = 3) +
  geom_boxplot(
    data = . %>% filter(str_detect(type, "Control")),
    aes(color = type),
    x = 0, position = position_dodge(0)
  ) +
  geom_hline(
    data = outliers, aes(yintercept = cutoff, color = type),
    linetype = "dashed"
  ) +
  scale_color_manual(values = c("#FAD510", "#CB2314", "#273046")) +
  facet_wrap(facets = vars(assay_type, time_point), scales = "free_y") +
  theme_nw2() +
  NULL

# remove 5 outliers in the controls
trimmed <- dr_data %>%
  left_join(outliers) %>%
  mutate(outlier = case_when(
    assay_type == "viability" & type == "Negative Control" & value > cutoff ~ TRUE,
    assay_type == "viability" & type == "Positive Control" & value < cutoff ~ TRUE,
    assay_type == "motility" & type == "Negative Control" & value < cutoff ~ TRUE,
    assay_type == "motility" & type == "Positive Control" & value > cutoff ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  filter(outlier == FALSE)

# correct for drift and normalize -----------------------------------------

# show the drift
drift <- trimmed %>%
  select(assay_date:value) %>%
  mutate(smooth_row = as.numeric(factor(row))) %>%
  ggplot(aes(x = smooth_row, y = value)) +
  geom_point(aes(color = plate), size = 2) +
  facet_grid(cols = vars(col), rows = vars(assay_type), scales = "free_y") +
  labs(x = "Row", y = "Value") +
  theme_nw2() +
  NULL

c_model <- trimmed %>%
  mutate(smooth_row = as.numeric(factor(row))) %>%
  mutate(smooth_col = as.numeric(col)) %>%
  filter(str_detect(treatment, "Control")) %>%
  group_by(time_point, assay_type) %>%
  group_nest() %>%
  mutate(model = map(data, ~ lm(value ~ smooth_row + smooth_col, data = .x))) %>%
  mutate(tidy = map(model, tidy)) %>%
  unnest(tidy) %>%
  select(time_point, assay_type, term, estimate) %>%
  filter(term != "(Intercept)") %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(row_coef = smooth_row, col_coef = smooth_col)

corrected <- trimmed %>%
  mutate(smooth_row = as.numeric(factor(row))) %>%
  mutate(smooth_col = as.numeric(col)) %>%
  left_join(c_model) %>%
  mutate(corrected_value = value - smooth_row * row_coef - smooth_col * col_coef)

corrected %>%
  ggplot(aes(x = smooth_row, y = value)) +
  geom_point(aes(color = plate), size = 2) +
  facet_grid(cols = vars(col), rows = vars(assay_type), scales = "free_y") +
  labs(x = "Row", y = "Value") +
  theme_nw2() +
  NULL


# rescale -----------------------------------------------------------------

control_summary <- corrected %>%
  filter(str_detect(treatment, "Control") == TRUE) %>%
  group_by(assay_date, plate, assay_type, treatment) %>%
  summarise(control_mean = mean(corrected_value)) %>%
  pivot_wider(names_from = treatment, values_from = control_mean) %>%
  rename(n_mean = 4, p_mean = 5)

rescaled_data <- corrected %>%
  left_join(., control_summary) %>%
  mutate(rescaled = (corrected_value - n_mean) / (p_mean - n_mean)) %>%
  select(-contains("mean"))


# 1 uM repeated -----------------------------------------------------------

# Of the 32 original hits, how many were "hits" (>25% increase in viability/decrease in movement compared to control) again in the DR?
primary_hits <- read_csv(here("Fig2/tables/Table1.csv"))

repeated <- rescaled_data %>%
  left_join(select(primary_hits, treatment, hit_type)) %>%
  filter(
    conc %in% c("1p", "HK", "1uM"),
    time_point == "36hr",
    !str_detect(treatment, "Control"),
    !is.na(hit_type)
  )

(repeated_plot <- repeated %>%
    mutate(Phenotype = str_to_sentence(assay_type),
           treatment = str_remove(treatment, " hydrochloride| maleate| citrate")) %>%
    ggplot(aes(x = treatment, y = rescaled, color = Phenotype, fill = Phenotype)) +
    geom_errorbar(
      position = position_dodge(width = 1),
      stat = "summary", width = 0.5, alpha = 0.75
    ) +
    geom_bar(position = "dodge", stat = "summary", fun = "mean", alpha = 0.95) +
    geom_quasirandom(
      shape = 21, fill = "white", size = 1,
      alpha = 0.75, dodge.width = 0.9, width = 0.1,
      show.legend = FALSE
    ) +
    geom_hline(yintercept = 0.25, color = 'grey50', linetype = 'dashed', size = 1) +
    scale_y_continuous(limits = c(-0.35, 1.32), breaks = seq(-0.2, 1.2, 0.2), labels = scales::percent) +
    scale_fill_manual(values = c("#A62C5F", "#35872E")) +
    scale_color_manual(values = c("#A62C5F", "#35872E")) +
    facet_grid(cols = vars(hit_type), scales = "free_x", space = "free_x") +
    labs( y = "Inhibition<br>(percent of control)") +
    theme_nw2() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = 'right'
    ) + 
    NULL
)

write_rds(repeated_plot, here('Fig2', 'supplementary', 'SupplementaryFig2a.rds'))
save_plot(here('Fig2', 'supplementary', 'SupplementaryFig2a.pdf'), 
          repeated_plot, base_width = 8, base_height = 4)

repeated_summary <- repeated %>%
  group_by(treatment, assay_type) %>%
  summarise(mean = mean(rescaled)) %>%
  filter(mean > 0.25) %>%
  ungroup()

tally(distinct(repeated_summary, treatment))

filter(primary_hits, !treatment %in% repeated_summary$treatment)

no_dr <- primary_hits %>%
  filter(!str_detect(treatment, "Control")) %>%
  filter(!treatment %in% dr_data$treatment)

# dose-response -----------------------------------------------------------

dr_plot <- rescaled_data %>%
  filter(!str_detect(treatment, "Control")) %>%
  mutate(conc = str_remove(conc, "uM")) %>%
  ggplot(aes(x = as.numeric(conc), y = rescaled, interaction(time_point, assay_type))) +
  geom_point(aes(shape = assay_type), alpha = 0.75) +
  scale_x_log10(limits = c(0.001, 1)) +
  scale_color_manual(
    limits = c("Non-hit", "DMSO", "Ivermectin", "Both", "Male", "Female"),
    values = c("#273046", "#FAD510", "#CB2314", "#00AFB5", "#EE92C2", "#9D75CB")
  ) +
  facet_wrap(facets = vars(treatment)) +
  labs(x = "Concentration (µM)", y = "Motility (Percent of control)") +
  theme_nw() +
  NULL

# fit a 4-parameter log-logistic model with the asymptotes set to 0 and 1
fit <- rescaled_data %>%
  filter(!str_detect(treatment, "Control")) %>%
  mutate(
    conc = as.numeric(str_remove(conc, "uM")),
    keep = case_when(
      treatment == "Bay 11-7085" ~ TRUE,
      # treatment == 'LDN 193189 dihydrochloride' ~ TRUE,
      treatment == "Polygodial" & assay_type == "motility" ~ TRUE,
      treatment == "SMER 3" ~ TRUE,
      treatment == "SP 100030" & assay_type == "motility" ~ TRUE,
      treatment == "SU 3327" ~ TRUE,
      treatment == "VAS 2870" & assay_type == "motility" ~ TRUE,
      treatment == "Nexinhib20" & time_point == "36hr" & assay_type == "motility" ~ TRUE,
      treatment == "N106" & time_point == "36hr" & assay_type == "motility" ~ TRUE,
      treatment == "MLS 1547" ~ TRUE,
      treatment == "Ciclopirox" & time_point == "36hr" & assay_type == "motility" ~ TRUE,
      treatment == "JIB 04" & time_point == "36hr" ~ TRUE,
      treatment == "NAV 2729" ~ TRUE,
      treatment == "NSC 319726" & assay_type == "motility" ~ TRUE,
      treatment == "SC 144 hydrochloride" & time_point == "36hr" & assay_type == "motility" ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(keep == TRUE) %>%
  group_by(treatment, time_point, assay_type) %>%
  group_nest() %>%
  mutate(drc = map(data, ~ drc::drm(.x$rescaled ~ .x$conc, fct = LL.4(fixed = c(NA, 0, 1, NA))))) %>%
  mutate(glance = map(drc, glance)) %>%
  mutate(tidy = map(drc, tidy)) %>%
  # get the EC50
  mutate(ec50 = map(drc, ~ ED(.x, 50)))

ec50 <- fit %>%
  unnest_longer(ec50) %>%
  select(-drc, -data)

temp <- as_tibble(as.vector(ec50$ec50) %>% matrix(nrow(ec50), 2))

fit_e <- bind_cols(fit, temp) %>%
  rename(estimate = V1, std_error = V2) %>%
  select(-ec50) %>%
  # left_join(select(adult_hits, treatment, adult_hit_type)) %>%
  mutate(estimate = case_when(
    estimate > 1 ~ ">1 µM",
    TRUE ~ str_c(round(estimate, 3), " µM")
  )
  )

newdata <- tibble(conc = exp(seq(log(0.001), log(1), length = 100)))

# add prediction for line/ribbon drawing
predict <- fit %>%
  mutate(pred = map(drc, ~ augment(.x, newdata = newdata, conf.int = TRUE))) %>%
  select(treatment, assay_type, time_point, pred) %>%
  unnest(pred) %>%
  mutate(
    assay_type = str_to_sentence(assay_type)
  )

curve_plot <- rescaled_data %>%
  filter(!str_detect(treatment, 'Control')) %>% 
  mutate(conc = str_remove(conc, "uM")) %>%
  ggplot(aes(x = as.numeric(conc), y = rescaled)) +
  geom_point(aes(shape = assay_type), alpha = 0.75) +
  geom_line(data = predict, aes(
    x = conc, group = interaction(assay_type, time_point),
    y = .fitted, linetype = assay_type
  )) +
  geom_text_repel(
    data = fit_e,
    aes(label = estimate),
    direction = "y", x = log(0.001), y = 1, size = 2.5, min.segment.length = 10000
  ) +
  scale_x_log10(limits = c(0.001, 1)) +
  scale_color_manual(
    limits = c("Non-hit", "Both", "Male", "Female"),
    values = c("#273046", "#00AFB5", "#EE92C2", "#9D75CB")
  ) +
  scale_fill_manual(
    limits = c("Non-hit", "Both", "Male", "Female"),
    values = c("#273046", "#00AFB5", "#EE92C2", "#9D75CB")
  ) +
  facet_wrap(facets = vars(treatment)) +
  scale_y_continuous(limits = c(-.4, 1.3), breaks = seq(0, 1, 0.2), labels = scales::percent) +
  labs(
    x = "Concentration (µM)", y = "Phenotype (Percent of control)",
    shape = "Phenotype", linetype = "Phenotype",
    fill = "Hit type", color = "Hit type"
  ) +
  theme_nw2() +
  theme(legend.position = "bottom") +
  NULL


# ec50 bar ----------------------------------------------------------------

# base <- hits %>%
#   select(treatment, viability, motility) %>%
#   pivot_longer(c(viability, motility), names_to = 'assay_type', values_to = 'value') %>%
#   select(-value) %>%
#   mutate(time_point = '36hr')

ec50_bar <- ec50 %>%
  rename(estimate = ec50) %>%
  select(-glance, -tidy) %>%
  filter(time_point == "36hr") %>%
  mutate(estimate = estimate[, 1])

(ec50_bar_plot <- ec50_bar %>%
    group_by(assay_type, group = cut(estimate,
                                     breaks = c(0, 0.2, 0.5, 1, Inf),
                                     labels = c("0-200 nM", "200-500 nM", "500 nM - 1 µM", ">1 µM")
    )) %>%
    summarise(n = n()) %>%
    filter(!is.na(group)) %>%
    mutate(assay_type = str_to_sentence(assay_type)) %>%
    ggplot(aes(x = group, y = n, fill = assay_type)) +
    geom_label(aes(label = str_extract(assay_type, "^.")),
               size = 2, color = "white"
    ) +
    scale_fill_manual(values = c("#A62C5F", "#35872E")) +
    scale_y_continuous(limits = c(0, 9), breaks = seq(0, 8, 2)) +
    # scale_shape_manual(values = c('M', 'V')) +
    labs(x = "EC<sub>50</sub>", y = "N", color = "Assay type") +
    coord_flip() +
    theme_minimal(base_size = 16, base_family = "Helvetica") +
    theme(
      axis.text.x = element_text(size = 8, angle = 0, vjust = 1, hjust = 0.5),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_markdown(angle = 90, size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.5),
      axis.line = element_line(size = 0.25, color = "black"),
      axis.ticks = element_line(size = 0.25, color = "black"),
      legend.position = "empty"
    ) +
    NULL)

write_rds(ec50_bar_plot, here("Fig2/subplots/ec50_bar.rds"))
save_plot(here('Fig2/subplots/Fig2i.pdf'), ec50_bar_plot,
          base_height = 2, base_width = 3)


# time points -------------------------------------------------------------

comp_data <- rescaled_data %>% 
  filter(assay_type == 'motility',
         treatment %in% ec50$treatment) %>% 
  mutate(conc = as.numeric(str_remove(conc, 'uM')),
         treatment = str_remove(treatment, " hydrochloride| maleate| citrate"),
         time_point = str_remove(time_point, 'hr')) 

comp_predict <- predict %>% 
  filter(assay_type == 'Motility',
         treatment %in% ec50$treatment) %>% 
  select(treatment:conc, -assay_type, -.lower, -.upper) %>% 
  distinct() %>% 
  mutate(treatment = as.factor(str_remove(treatment, " hydrochloride| maleate| citrate")))

comp_ribbon <- comp_predict %>% 
  select(treatment, conc, time_point, .fitted) %>% 
  pivot_wider(names_from = time_point, values_from = .fitted) %>% 
  replace_na(list(`12hr` = 0))


auc <- function(x, y) {
  
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2 
  
}

area_between_curves <- comp_predict %>% 
  group_nest(treatment, time_point) %>% 
  mutate(auc = map(data, ~auc(.x$conc, .x$.fitted))) %>%  
  select(-data) %>% 
  mutate(auc = unlist(auc)) %>% 
  pivot_wider(names_from = time_point, values_from = auc, values_fill = 0) %>% 
  mutate(abc = `36hr` - `12hr`) %>% 
  mutate(treatment = fct_reorder(treatment, abc, .desc = TRUE))


(comp_plot <- comp_predict %>% 
    filter(!treatment %in% c('N106', 'Nexinhib20', 'JIB 04', 'SC 144', 'Ciclopirox')) %>% 
    mutate(treatment = factor(treatment, levels = levels(area_between_curves$treatment))) %>%
    ggplot(aes(x = conc)) +
    geom_ribbon(data = comp_ribbon %>% filter(!treatment %in% c('N106', 'Nexinhib20', 'JIB 04', 'SC 144', 'Ciclopirox')),
                aes(x = conc, ymin = `12hr`, ymax = `36hr`), 
                fill = 'grey', alpha = 0.6) +
    geom_line(aes(y = .fitted, color = time_point),
              size = 1) +
    geom_text(data = area_between_curves %>% filter(!treatment %in% c('N106', 'Nexinhib20', 'JIB 04', 'SC 144', 'Ciclopirox')),
              aes(x = 0.8, y = 0.15, label = round(abc, digits = 3))) +
    labs(x = "Concentration (µM)", y = "Motility inhibition (percent of control)", color = 'Hours post-treatment') +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_color_manual(values = c("#F1A340", "#998EC3")) +
    facet_wrap(~treatment) +
    theme_nw2() +
    theme(
      axis.text.x = element_markdown(angle = 0, hjust = 0.5)
    ) +
    NULL)

write_rds(comp_plot, here('Fig2', 'supplementary', 'SupplementaryFig2b.rds'))
save_plot(here('Fig2/supplementary/SupplementaryFig2b.pdf'), comp_plot, base_width = 8, base_height = 6)


# export data -------------------------------------------------------------

bind_cols(fit, temp) %>%
  rename(estimate = V1, std_error = V2) %>%
  mutate(name = str_c(assay_type, time_point, sep = '_')) %>% 
  select(treatment, name, estimate) %>%
  mutate(treatment = str_remove(treatment, " hydrochloride| maleate| citrate")) %>%
  pivot_wider(names_from = name, values_from = estimate) %>% 
  write_rds(here("Fig2/data/ec50.rds"))
