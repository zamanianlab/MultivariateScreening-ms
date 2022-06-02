# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(cowplot)
library(ggbeeswarm)
library(ggtext)
library(ggridges)
library(ZamanianLabThemes)

# stats
library(broom)

# misc
library(conflicted)
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")

# import data -------------------------------------------------------------

data_dir <- here("Fig3/data/")

data_files <- tibble(base = data_dir, stem = list.files(path = data_dir,
                                                        pattern = ".*.csv$",
                                                        recursive = TRUE)) %>%
  filter(!str_detect(stem, "202[0-9]{5}")) %>%
  mutate(path = str_c(base, stem, sep = "")) %>%
  separate(stem, sep = "/", into = c("species", "exp_design", "stages", "time_point")) %>%
  mutate(time_point = str_remove(time_point, "hr.csv")) %>%
  select(path, base, species, stages, exp_design, time_point)

get_data <- function(...) {
  df <- tibble(...)

  data <- read_csv(df$path,
    col_names = c("well", "raw_flow", "worm_area", "norm_flow"),
    col_types = "cddd", skip = 1
  ) %>%
    mutate(
      stages = df$stages,
      time_point = df$time_point,
      species = df$species,
      exp_design = df$exp_design
    )
}

data <- data_files %>%
  pmap_dfr(get_data) %>%
  mutate(exp_date = case_when(
    species == "Bpa" & exp_design == "incomplete" ~ "20210202",
    species == "Bpa" & exp_design == "complete" ~ "20201110",
    species == "Bma" & exp_design == "incomplete" ~ "20210120",
    species == "Bma" & exp_design == "complete" ~ "20210223"
  )) %>%
  filter(case_when(
    exp_date == "20210120" & stages == "AF" & species == "Bma" ~ well != "C0003",
    exp_date == "20210202" & stages == "AM" & species == "Bpa" ~ !well %in% c("A03", "D03"),
    TRUE ~ TRUE
  )) %>%
  # some videos were recorded at 30 FPS instead of 16
  mutate(raw_flow = case_when(
    time_point %in% c(48, 72, 96, 120) & species == "Bma" & exp_design == "incomplete" ~ raw_flow * 16 / 30,
    time_point %in% c(72, 96, 120) & species == "Bma" & exp_design == "complete" ~ raw_flow * 16 / 30,
    TRUE ~ raw_flow
  )) %>%
  rename(total_motility = raw_flow, norm_motility = norm_flow)


# analysis ----------------------------------------------------------------

# get the max and min flow for each stage
range <- data %>%
  group_by(time_point, stages) %>%
  summarize(ceiling = max(total_motility), floor = min(total_motility)) %>%
  ungroup() %>%
  pivot_longer(ceiling:floor, names_to = "measure", values_to = "total_motility") %>%
  filter(
    (time_point == "0" & measure == "ceiling") |
      (time_point == "96" & measure == "floor")
  ) %>%
  select(-time_point) %>%
  pivot_wider(names_from = measure, values_from = total_motility)

# filter anything that is less than the floor + 5%
female_cutoff <- as.numeric(range[1, 3]) * 1.05
male_cutoff <- as.numeric(range[2, 3]) * 1.05

dead <- data %>%
  filter(case_when(
    time_point == '0' & stages == 'AF' & total_motility < female_cutoff ~ TRUE,
    time_point == '0' & stages == 'AM' & total_motility < male_cutoff ~ TRUE,
  )) %>%
  select(exp_date, species, stages, well) %>%
  mutate(dead = TRUE)

filtered_data <- data %>%
  left_join(dead) %>%
  filter(is.na(dead))

# rescale to the floor/ceiling
scaled <- filtered_data %>%
  left_join(range) %>%
  mutate(rescaled = (total_motility - floor) / (ceiling - floor)) %>%
  select(-ceiling, -floor)

# filter to keep controls for each experiment
self_norm <- scaled %>%
  filter(time_point == "0") %>%
  select(exp_design, well, species, stages, init_motility = total_motility, init_rescaled = rescaled)

# self-normalize
normalized <- scaled %>%
  left_join(self_norm) %>%
  mutate(
    self_norm_rescaled = rescaled / init_rescaled,
    self_norm_motility = total_motility / init_motility
  )

(normalized %>%
  ggplot(aes(x = time_point, y = self_norm_rescaled)) +
  geom_quasirandom() +
  facet_grid(cols = vars(stages)) +
  theme_nw2() +
  NULL)

# prune a few outliers (any worm that ever has >4x its initial motility)
pruned <- normalized %>%
  group_by(species, stages, well) %>%
  mutate(outlier = case_when(
    any(self_norm_rescaled > 4) ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  filter(outlier == FALSE) %>%
  select(exp_design, well, species, stages, time_point, rescaled, self_norm_rescaled)


# plot --------------------------------------------------------------------

(plot <- pruned %>%
  mutate(
    stages = case_when(
      stages == "AF" ~ "Female",
      stages == "AM" ~ "Male"
    ),
    exp_design = factor(exp_design,
      levels = c("complete", "incomplete"),
      labels = c("Complete", "Incomplete")
    ),
    species = factor(species,
      levels = c("Bpa", "Bma"),
      labels = c("*B. pahangi*", "*B. malayi*")
    )
  ) %>%
  ggplot(aes(x = as.numeric(time_point), y = self_norm_rescaled, color = exp_design, fill = exp_design)) +
  geom_line(aes(group = interaction(stages, well)), alpha = 0.3, size = 0.25) +
  geom_hline(yintercept = 1, color = "black", alpha = 0.5, linetype = "dashed") +
  geom_smooth(aes(group = interaction(stages, exp_design)), size = 1, se = FALSE) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96, 120)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  facet_grid(rows = vars(stages), cols = vars(species)) +
  labs(x = "Time point (hr)", y = "Self-normalized motility", color = "Media", fill = "Media") +
  theme_nw2() +
  theme(
    # panel.grid.major.x = element_blank(),
    axis.text.x = element_markdown(angle = 0, hjust = 0.5)
  ) +
  NULL)

save_plot(here("Fig3/subplots/Fig3a.pdf"), plot, base_width = 7, base_height = 4)
write_rds(plot, here("Fig3/subplots/optimized_motility.rds"))
