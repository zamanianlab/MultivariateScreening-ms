# data wrangling/plotting
library(tidyverse)

# other plotting
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(ZamanianLabThemes)

# stats
library(parsnip)
library(tidymodels)
library(drc)
library(broom)

# misc
library(conflicted)
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# import ------------------------------------------------------------------

pruned_data <- read_rds(here('Fig5/data/hypoxia_data.rds'))

# analysis ----------------------------------------------------------------

### trim outliers

coef <- 2.5

outliers <- pruned_data %>%
  group_by(treatment, other, strains, conc) %>%
  group_nest() %>%
  mutate(quant = map(data, ~ as.numeric(quantile(.x$area_shape_major_axis_length, c(0, 0.25, 0.5, 0.75, 1)))),
         iqr = unlist(map(quant, ~ diff(.x[c(2, 4)]))),
         top_cutoff = unlist(map(quant, ~ pluck(.x, 4))) + coef * iqr,
         bottom_cutoff = unlist(map(quant, ~ pluck(.x, 2))) - coef * iqr) %>% 
  select(treatment, other, strains, conc, contains('cutoff'))

trimmed_data <- pruned_data %>% 
  left_join(outliers) %>% 
  mutate(outlier = case_when(
    area_shape_major_axis_length > top_cutoff | area_shape_major_axis_length < bottom_cutoff ~ TRUE,
    TRUE ~ FALSE
  ), .before = image_number) %>% 
  # filter(outlier == FALSE) %>%
  select(-contains('cutoff'), -outlier)

### QC
trimmed_data %>%
  filter(treatment %in% c('DMSO', 'Untreated')) %>%
  ggplot(aes(x = treatment, y = area_shape_major_axis_length, color = other)) +
  geom_beeswarm(dodge.width = 1, alpha = 0.2) +
  stat_summary(geom = 'point', fun = median, size = 3, shape = 'circle',
               position = position_dodge(1)) +
  stat_summary(geom = 'point', fun = mean, size = 3, shape = 'square',
               position = position_dodge(1)) +
  theme_nw2() +
  NULL

# sqrt transform of poisson distribution
trans_data <- trimmed_data %>%
  mutate(trans_length = sqrt(area_shape_major_axis_length))

control_summary <- trans_data %>%
  filter(treatment == 'DMSO') %>%
  group_by(plate) %>%
  summarise(control_median = median(trans_length, na.rm = FALSE),
            control_mean = mean(trans_length, na.rm = FALSE),
            untrans_control_median = median(area_shape_major_axis_length, na.rm = FALSE))

# normalize by dividing by the median of the control (DMSO)
norm_data <- trans_data %>% 
  filter(treatment != 'Untreated') %>% 
  mutate(conc = case_when(
    treatment == 'DMSO' ~ 0.01,
    treatment != 'DMSO' ~ as.numeric(str_remove(conc, 'uM'))
  )) %>% 
  left_join(control_summary) %>% 
  mutate(norm_length = trans_length / control_mean,
         untrans_norm_length = area_shape_major_axis_length / untrans_control_median)

well_summary <- norm_data %>% 
  group_by(metadata_date, plate, other, strains, well, treatment, conc) %>% 
  summarise(mean_norm_length = mean(norm_length),
            median_norm_length = median(norm_length),
            mean_length = mean(untrans_norm_length),
            median_length = median(untrans_norm_length),
            median_raw = median(area_shape_major_axis_length))

norm_data %>%
  filter(
    treatment == 'ML 228'
  ) %>%
  ggplot(aes(x = conc, y = untrans_norm_length)) +
  geom_quasirandom(alpha = 0.5, color = 'indianred') +
  geom_quasirandom(data = well_summary %>% filter(treatment == 'ML 228'),
             aes(x = conc, y = mean_length),
             color = 'steelblue', size = 3) +
  stat_summary(geom = 'point', fun = median, size = 3, shape = 'circle', color = 'black') +
  stat_summary(geom = 'point', fun = mean, size = 3, shape = 'square', color = 'black') +
  scale_x_log10() +
  facet_grid(rows = vars(other)) +
  theme_nw2() +
  NULL


### dose-responses
curves <- norm_data %>% 
  filter(!treatment %in% c('DMSO', 'Untreated')) %>%
  # uncomment below to fit to the summarized well data
  group_by(treatment, strains, other, well, conc) %>%
  summarise(mean_norm_length = mean(norm_length),
            median_norm_length = median(norm_length),
            mean_length = mean(untrans_norm_length),
            median_length = median(untrans_norm_length)) %>%
  ungroup() %>%
  group_nest(treatment, strains, other) %>% 
  mutate(drc = map(data, ~ drm(.x$mean_length ~ .x$conc, fct = drc::LL.4()))) %>%
  mutate(glance = map(drc, glance)) %>%
  mutate(tidy = map(drc, tidy)) %>%
  # get the EC50
  mutate(ec50 = map(drc, ~ drc::ED(.x, 50)))

# get the ec50
ec50 <- curves %>%
  select(ec50) %>%
  unnest(c(ec50)) 

temp <- as_tibble(as.vector(ec50$ec50) %>% matrix(nrow(ec50), 2))

fit_e <- bind_cols(curves, temp) %>%
  rename(estimate = V1, std_error = V2) %>%
  select(treatment, strains, other, estimate, std_error)

# fit the model to draw a line
newdata <- expand.grid(
  plate = unique(trimmed_data$plate),
  strains = unique(trimmed_data$strains),
  other = unique(trimmed_data$other),
  conc = exp(seq(log(0.39), 
                 log(max(norm_data$conc, na.rm = TRUE)), 
                 length = 100))) %>% 
  group_nest(plate, strains, other) %>% 
  rename(newdata = data)

predict <- curves %>%
  left_join(newdata) %>% 
  mutate(pred = map2(drc, newdata, ~augment(.x, newdata = .y, conf.int = TRUE))) %>% 
  unnest(pred)

(curve_plot <- norm_data %>%
    filter(treatment == 'ML 228') %>% 
    ggplot() +
    geom_line(data = predict, 
              aes(x = conc, group = other, y = .fitted),
              size = 1) +
    geom_quasirandom(data = well_summary %>% filter(treatment == 'ML 228'),
                     aes(x = conc, y = mean_length),
                     size = 0.6, width = 0.1, alpha = 0.6) +
    # geom_ribbon(data = predict,
    #             aes(x = conc, ymin = .lower, ymax = .upper),
    #             alpha = 0.25) +
    geom_rect(data = fit_e,
              aes(xmin = estimate - std_error, xmax = estimate + std_error, ymin = -Inf, ymax = Inf),
              fill = 'black', alpha = 0.25) +
    geom_vline(data = fit_e, aes(xintercept = estimate), linetype = 'dashed') +
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    # scale_color_manual(values = c('indianred', '#275D7C', '#B39C4D')) +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(from = 0.5, to = 1.5, by = 0.5)) +
    coord_cartesian(clip = "off") +
    labs(x = 'Concentration (µM)', y = "Normalized length",
         color = 'Strain', fill = 'Strain', shape = 'Date') +
    facet_grid(rows = vars(other)) +
    theme_nw2() +
    theme(
      legend.position = 'right',
      axis.text.x = element_markdown(angle = 0, hjust = 0.5),
      axis.title.x = element_markdown(face = 'plain'),
      axis.title.y = element_markdown(face = 'plain')) +
    remove_legend() +
    NULL)

save_plot(here('Fig5', 'subplots', 'Fig5e.pdf'),
          curve_plot,
          base_height = 3, base_width = 4)

write_rds(curve_plot, here('Fig5', 'subplots', 'k_s_media.rds'), compress = 'gz')





