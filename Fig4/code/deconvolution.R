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

# misc
library(conflicted)
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")


# import ------------------------------------------------------------------
pruned_data <- read_rds(here('Fig4/data/deconvolution_data.rds'))


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
  left_join(control_summary) %>% 
  mutate(norm_length = trans_length / control_mean,
         untrans_norm_length = area_shape_major_axis_length / untrans_control_median)

plot_data <- norm_data %>% 
    mutate(
      conc = factor(conc, 
                    levels = c('1p', '0.5uM', '1uM', '2.5uM', '5uM', '10uM'),
                    labels = c('1% DMSO', '500 nM', '1 µM', '2.5 µM', '5 µM', '10 µM')),
      genotype = factor(genotype,
                        levels = c('N2', '*dat-1<br>(ok157)*', '*ser-1<br>(ok345)*',
                                   '*ador-1<br>(tm3971)*', '*unc-103<br>(e1597n1213)*', '*hif-1<br>(mr20)*',
                                   '*arf-6<br>(tm1447)*'))) %>% 
  filter(!plate %in% c('20220120-p04-NJW_1105', '20220120-p05-NJW_1106', '20220121-p04-NJW_1112', '20220121-p05-NJW_1113'))# %>%
    
(plot <- plot_data %>% 
    filter(other == '48hr',
           panel != 'Polygodial') %>% 
    ggplot(aes(x = genotype, y = untrans_norm_length)) + 
    geom_violin(aes(fill = conc), draw_quantiles = c(0.25, 0.5, 0.75), color = 'white', size = 0.2,
                position = "dodge", scale = 'width', trim = TRUE) +
    # geom_quasirandom(aes(color = conc),
    #                  size = 1, shape = 21, fill = 'white',
    #                  alpha = 0.75, dodge.width = 0.9, width = 0.2) +
    # geom_text(
    #   data = stat_layer, aes(label = sig, y = y, color = conc),
    #   size = 8, position = position_dodge(1), show.legend = FALSE
    # ) +
    scale_y_continuous(breaks = c(0, 0.5, 1, 1.5),
                       limits = c(0, 1.75)) +
    scale_color_manual(values = c('#333333', '#ADD9F4', '#78290F', '#FF7D00')) +
    scale_fill_manual(values = c('#333333', '#ADD9F4', '#78290F', '#FF7D00')) +
    facet_grid(cols = vars(panel), scales = 'free', space = 'free') +
    labs(y = 'Normalized<br>length', x = '', color = 'Concentration', fill = 'Concentration') +
    theme_nw2() +
    theme(
      axis.text.x = element_markdown(angle = 0, hjust = 0.5)
    ) +
    NULL)

save_plot(here('Fig4/subplots/Fig4c.pdf'), plot, base_width = 8, base_height = 6)
write_rds(plot, here('Fig4/subplots/deconvolution.rds'))





                     