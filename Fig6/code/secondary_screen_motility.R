# data wrangling/plotting
library(tidyverse)

# other plotting
library(cowplot)
library(ggtext)
library(ZamanianLabThemes)
library(ggbeeswarm)
library(ggrepel)

# stats
library(broom)

# misc
library(conflicted)
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")


# prune outliers ----------------------------------------------------------

tidy_data <- read_rds(here("Fig6/data/motility_data.rds")) %>% 
  mutate(total_motility = log10(total_motility))

scores <- read_csv(here("Fig6/data/motility_scores.csv"),
                   col_names = c("well", "plate", "score")
)

# subtract the background by getting the mean of wormless
background <- tidy_data %>%
  filter(is.na(species)) %>%
  pull(total_motility) %>%
  mean()

tidy_data <- tidy_data %>%
  mutate(total_motility = total_motility - background) %>%
  mutate(total_motility = case_when(
    total_motility < 0 ~ 0,
    TRUE ~ total_motility
  )) %>%
  drop_na(treatment)  

init <- tidy_data %>%
  filter(time_point == "0hr") %>%
  left_join(scores)

ggplot(
  init,
  aes(
    y = total_motility,
    fill = as.factor(score),
    color = as.factor(score)
  )
) +
  geom_beeswarm(aes(x = stages))

dead <- init %>%
  filter(score == 1) %>%
  select(species, stages, well, treatment) %>%
  mutate(dead = TRUE)

filtered_data <- tidy_data %>%
  left_join(dead) %>%
  filter(is.na(dead))

# analysis ----------------------------------------------------------------

# filter to keep controls for each experiment
self_norm <- filtered_data %>%
  filter(time_point == "0hr") %>%
  select(well, species, stages, treatment, conc, init_motility = total_motility)

# self-normalize
normalized <- filtered_data %>%
  left_join(self_norm) %>%
  mutate(self_norm_motility = total_motility / init_motility) %>%
  filter(self_norm_motility < 2.5)

(normalized %>%
  ggplot(aes(x = time_point, y = total_motility)) +
  geom_quasirandom() +
  facet_grid(cols = vars(stages)) +
  theme_nw2() +
  NULL)


# plot --------------------------------------------------------------------

primary_hits <- read_csv(here("Fig2/tables/Table1.csv")) %>%
  select(number, treatment)

plot_data <- normalized %>%
  left_join(primary_hits) %>%
  mutate(
    species = factor(species,
      levels = c("Bpa", "Bma"),
      labels = c("*B. pahangi*", "*B. malayi*")
    ),
    stages = factor(stages,
      levels = c("AF", "AM"),
      labels = c("Adult female", "Adult male")
    ),
    smooth_time_point = case_when(
      time_point == "0hr" ~ 0,
      time_point == "0.1hr" ~ 1,
      time_point == "1hr" ~ 2,
      time_point == "24hr" ~ 3,
      time_point == "48hr" ~ 4
    ),
    conc = factor(conc,
      levels = c("0.1p", "1uM"),
      labels = c("0.1%", "1 uM")
    ),
    treatment = str_remove(treatment, " hydrochloride| maleate| citrate"),
    number = case_when(
      treatment %in% c("Ivermectin", "DMSO") ~ treatment,
      !is.na(number) ~ as.character(number),
      is.na(number) ~ NA_character_
    )
  ) %>%
  drop_na(number)

### initial plot

# get the DMSO means
dmso <- plot_data %>%
  filter(treatment == "DMSO", time_point == "0hr") %>%
  group_by(stages) %>%
  summarise(control_mean = mean(total_motility))

# get the DMSO mean at 48 hr
dmso_end <- plot_data %>%
  filter(treatment == "DMSO", time_point == "48hr") %>%
  group_by(stages) %>%
  summarise(control_mean = mean(total_motility))

# 48 hr bar plot
end_data <- plot_data %>%
  ungroup() %>%
  # left_join(adult_hits) %>%
  left_join(dmso_end) %>%
  filter(time_point == "48hr") %>%
  mutate(mean_norm = total_motility / control_mean)

end_hits <- end_data %>%
  mutate(individual_hit = case_when(
    stages == "Adult male" & mean_norm < 0.5 ~ TRUE,
    stages == "Adult female" & mean_norm < 0.5 ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  ungroup() %>%
  group_by(stages, treatment) %>%
  summarise(
    treatment_mean = mean(mean_norm),
    n_hit = sum(individual_hit),
    n = n()
  ) %>%
  mutate(hit = case_when(
    treatment %in% c("Ivermectin", "DMSO") ~ treatment,
    stages == "Adult male" & (n_hit / n > 0.5) & treatment_mean < 0.6 ~ "Hit",
    stages == "Adult female" & (n_hit / n > 0.5) & treatment_mean < 0.4 ~ "Hit",
    TRUE ~ "Non-Hit"
  ))

(end_bar <- end_data %>%
  left_join(end_hits) %>%
  # hackish way of reordering by group
  mutate(reorder = factor(str_c(treatment, stages, sep = "_"))) %>%
  mutate(reorder = fct_reorder(reorder, total_motility, mean, .desc = TRUE)) %>%
  mutate(reorder = fct_relevel(reorder, "DMSO_Adult female", "DMSO_Adult male", "Ivermectin_Adult female", "Ivermectin_Adult male")) %>%
  ggplot(aes(x = reorder, y = mean_norm, color = hit, fill = hit)) +
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
  scale_x_discrete(labels = ~ str_remove(.x, "_.*")) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(
    limits = c("Non-Hit", "DMSO", "Ivermectin", "Hit"),
    values = c("#273046", "#FAD510", "#CB2314", "#00AFB5")
  ) +
  scale_fill_manual(
    limits = c("Non-Hit", "DMSO", "Ivermectin", "Hit"),
    values = c("#273046", "#FAD510", "#CB2314", "#00AFB5")
  ) +
  labs(x = "Treatment", y = "48 hr. motility<br>(Percent control)", color = "", fill = "") +
  facet_wrap(facets = vars(stages), nrow = 2, scales = "free") +
  theme_nw2() +
  theme(
    panel.grid.major.x = element_blank(),
    strip.text.x = element_markdown(face = "italic")
  ) +
  NULL)

save_plot(here("Fig6/subplots/Fig6a.pdf"), end_bar, base_width = 8, base_height = 6)
write_rds(end_bar, here("Fig6/subplots/motility.rds"))

# mean motility at each time point
treatment_summary <- plot_data %>%
  group_by(smooth_time_point, time_point, stages, treatment) %>%
  summarize(mean = mean(total_motility)) %>%
  left_join(end_hits)

# background ribbons showing the SD for each time point
# ribbon <- plot_data %>%
#   left_join(., adult_hits) %>%
#   filter(is.na(hit_type)) %>%
#   group_by(smooth_time_point, species, stages) %>%
#   summarize(
#     sd = sd(self_norm_rescaled, na.rm = TRUE),
#     n = n(),
#     se = sd / sqrt(n))

# add if text labels are going to be added to the right of the plot
# blank <- tibble(smooth_time_point = 5)

### overlapping lines
# line_plot <- plot_data %>%
#   left_join(., adult_hits) %>%
#   ggplot(aes(x = smooth_time_point)) +
#   # geom_blank(data = blank) +
#   geom_ribbon(data = ribbon, aes(ymin = 1 - sd, ymax = 1 + sd),
#               fill = '#354823', alpha = 0.4) +
#   # geom_hline(yintercept = 0.5, linetype = 'dashed') +
#   geom_line(data = filter(treatment_summary, hit_type == 'Non-hit'),
#             aes(y = mean, group = interaction(treatment, stages)),
#             alpha = 0.2, color = '#273046') +
#   # white shadow
#   geom_line(data = filter(treatment_summary, hit_type != 'Non-hit'),
#             aes(y = mean, group = interaction(treatment, stages)),
#             color = 'white', size = 0.5) +
#   # real data
#   geom_line(data = filter(treatment_summary, hit_type != 'Non-hit'),
#             aes(y = mean, group = interaction(treatment, stages),
#                 color = hit_type), alpha = 0.9, size = 0.3) +
#   geom_point(data = filter(treatment_summary, mean < 0.5, hit_type != 'Non-hit') %>% ungroup() %>% arrange(smooth_time_point) %>% distinct(stages, treatment, .keep_all = TRUE),
#              aes(y = mean, group = interaction(treatment, stages), fill = hit_type),
#              color = 'white', shape = 21, alpha = 0.9) +
#   # geom_text_repel(data = filter(treatment_summary, hit_type %in% c('Both', 'Male', 'Female'), time_point == '48hr'),
#   #                 aes(label = treatment, color = hit_type, x = 4.5, y = mean),
#   #                 size = 1, direction = 'y', seed = 1234, show.legend = FALSE,
#   #                 hjust = 0, min.segment.length = 1000) +
#   scale_x_continuous(labels = c('\u00d8', '0', '1', '24', '48')) +
#   scale_color_manual(limits = c('Non-hit', 'DMSO', 'Ivermectin', 'Both', 'Male', 'Female'),
#                      values = c('#273046', '#FAD510', '#CB2314', '#00AFB5', '#EE92C2', '#9D75CB')) +
#   scale_fill_manual(limits = c('Non-hit', 'DMSO', 'Ivermectin', 'Both', 'Male', 'Female'),
#                     values = c('#273046', '#FAD510', '#CB2314', '#00AFB5', '#EE92C2', '#9D75CB')) +
#   # facet_wrap(facets = vars(treatment), nrow = 6) +
#   facet_grid(rows = vars(stages)) +
#   labs(x = "Hours post-treatment", y = "Self-normalized motility",
#        color = 'Hit type', fill = 'Hit type',) +
#   theme_nc +
#   NULL
#
# save_plot('plots/overlapping_lines.pdf', line_plot, base_width = 3.935, base_height = 4)
# save_plot('plots/overlapping_lines.png', line_plot, base_width = 3.935, base_height = 4)
# write_rds(line_plot, 'plots/overlapping_lines.rds')

### distinct lines
(distinct_plot <- plot_data %>%
  left_join(end_hits) %>%
  filter(treatment %in% filter(end_hits, hit == "Hit")$treatment) %>%
  ggplot(aes(x = smooth_time_point, color = hit)) +
  geom_point(aes(y = self_norm_motility, shape = stages), size = 1) +
  geom_line(aes(y = self_norm_motility, group = interaction(treatment, stages, well)),
    alpha = 0.2
  ) +
  stat_summary(
    geom = "line", aes(
      y = self_norm_motility,
      group = interaction(treatment, stages)
    ),
    fun.data = "mean_se", size = 0.5
  ) +
  # geom_hline(yintercept = 0.5) +
  geom_richtext(
    data = . %>% ungroup() %>% distinct(treatment, hit),
    aes(label = str_c("*", treatment, "*")),
    x = 2, y = 2, color = "black", size = 3, show.legend = FALSE
  ) +
  scale_x_continuous(labels = c("\u00d8", "0", "1", "24", "48")) +
  scale_y_continuous(limits = c(-0.25, 2.25)) +
  scale_color_manual(
    limits = c("Hit", "Non-Hit"),
    values = c("#00AFB5", "#273046"), guide = "none"
  ) +
  facet_wrap(facets = vars(treatment), scales = "free_y", ncol = 2) +
  labs(x = "Hours post-treatment", y = "Self-normalized  motility", shape = "", color = "") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(size = 0.5), 
    axis.text.x = ggtext::element_markdown(face = "plain", 
                                           size = 8, angle = 45, vjust = 1, hjust = 1), 
    axis.text.y = ggtext::element_markdown(face = "plain", 
                                           size = 8), 
    axis.title.x = ggtext::element_markdown(face = "plain", 
                                                                                              size = 10),
    axis.title.y = ggtext::element_markdown(face = "plain", 
                                                                                                                                                  angle = 90, size = 10), 
    strip.text.x = element_blank(),
    strip.text.y = ggtext::element_markdown(size = 10), 
    axis.line = element_line(size = 0.25, color = "black"), 
    axis.ticks = element_line(size = 0.25, color = "black"), 
    legend.title = ggtext::element_markdown(face = "plain", 
                                            size = 10), legend.text = ggtext::element_markdown(size = 8), 
    legend.position = "right"
  ) +
  NULL)

distinct_plot <- lemon::reposition_legend(distinct_plot, "bottom right", panel = "panel-2-7")

save_plot(here("Fig6/subplots/Fig6b.pdf"), distinct_plot, base_width = 5, base_height = 6)
save_plot(here("Fig6/subplots/Fig6b.png"), distinct_plot, base_width = 5, base_height = 6)
write_rds(distinct_plot, here("Fig6/subplots/distinct_lines.rds"))


# export data -------------------------------------------------------------

out <- treatment_summary %>%
  pivot_wider(id_cols = c(stages, treatment), names_from = time_point, values_from = mean) %>%
  relocate(treatment, .before = stages) %>%
  write_rds(here("Fig6/data/adult_motility_data.rds"))

write_rds(end_hits, here("Fig6/data/motility_hits.rds"))

# write_csv(end_hits, "tables/adult_motility_hits.csv")

end_data <- plot_data %>%
  ungroup() %>%
  filter(time_point == "48hr") %>%
  select(treatment, stages, total_motility)

write_rds(end_data, here("Fig6/data/motility_end.rds"))
