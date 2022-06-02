# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(ggridges)
library(ZamanianLabThemes)
library(wesanderson)

# misc
library(conflicted)
library(here)

# conflict resolution
conflict_prefer("filter", "dplyr")


# import data and organize ------------------------------------------------

data <- read_rds(here("Fig3/data/fecundity_data.rds")) %>%
  separate(other, c("time_point", "exp_date"), "_")

# species differences -------------------------------------------------------

zero_hour_median <- data %>%
  filter(time_point == "0hr") %>%
  mutate(species = factor(species,
    levels = c("Bpa", "Bma"),
    labels = c("*B. pahangi*", "*B. malayi*")
  )) %>%
  group_by(species) %>%
  summarize(median = median(scaled), n = n()) %>%
  mutate(label = paste0(species, "<br>N = ", n))

stat_layer <- broom::tidy(kruskal.test(
  x = filter(data, time_point == "0hr")$scaled,
  g = filter(data, time_point == "0hr")$species
))

(zero_hour <- data %>%
  filter(time_point == "0hr") %>%
  mutate(species = factor(species,
    levels = c("Bpa", "Bma"),
    labels = c("*B. pahangi*", "*B. malayi*")
  )) %>%
  left_join(., zero_hour_median) %>%
  ggplot(aes(x = scaled, fill = label)) +
  geom_density(position = "identity", alpha = 0.9) +
  geom_vline(data = zero_hour_median, aes(xintercept = median, color = label), linetype = "dashed") +
  geom_richtext(
    data = stat_layer, aes(label = str_c("*p = ", formatC(p.value, format = "e", digits = 2), "*")),
    x = 1500, y = 0.00145, inherit.aes = FALSE, fill = NA, label.color = NA,
    size = 8 * 0.352777778
  ) + # the magic number for conversion
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c('black', 'grey')) +
  scale_color_manual(values = c('black', 'grey')) +
  labs(x = "Inferred progeny", y = "Density", fill = "Species", color = "Species") +
  theme_nw2() +
  NULL)

save_plot(here("Fig3/subplots/Fig3d.pdf"), zero_hour, base_width = 4, base_height = 4)
write_rds(zero_hour, here("Fig3/subplots/baseline_fecundity.rds"))

# optimization experiments ------------------------------------------------

optimization_plates <- tibble(
  plate = c(
    "20201012-p02-NJW_132", "20201019-p01-NJW_136", "20201019-p02-NJW_137",
    "20201203-p01-NJW_275", "20201203-p02-NJW_276", "20210115-p01-NJW_379",
    "20210115-p02-NJW_380", "20210201-p02-KTR_394", "20210201-p01-KTR_392",
    "20210208-p01-KTR_400", "20210208-p02-KTR_401", "20210125-p01-KTR_386",
    "20210201-p01-KTR_392", "20210303-p01-NJW_415", "20210303-p02-NJW_416",
    "20210303-p03-NJW_417", "20210322-p01-KTR_459", "20210322-p02-KTR_460",
    "20210830-p05-KTR_813", "20210830-p04-KTR_812", "20210830-p03-KTR_811",
    "20210830-p02-KTR_810", "20210830-p01-KTR_809"
  ),
  exp_design = c(
    "Drug<br>assay", "Drug<br>assay", "Drug<br>assay",
    "Complete<br>media", "Complete<br>media", "RNAi<br>timeline",
    "RNAi<br>timeline", "Drug<br>assay", "Drug<br>assay",
    "Drug<br>assay", "Drug<br>assay", "Drug<br>assay",
    "Drug<br>assay", "RNAi<br>timeline", "Complete<br>media",
    "", "Complete<br>media", "RNAi<br>timeline",
    "Drug<br>assay", "Drug<br>assay", "Drug<br>assay",
    "Drug<br>assay", "Drug<br>assay"
  )
)

optimization_data <- data %>%
  filter(plate %in% optimization_plates$plate) %>%
  left_join(optimization_plates) %>%
  mutate(
    exp_design = case_when(
      exp_design == "" & col %in% c("01", "02", "03", "04", "05", "06") ~ "Complete<br>media",
      exp_design == "" & col %in% c("07", "08", "09", "10", "11", "12") ~ "RNAi<br>timeline",
      TRUE ~ exp_design
    ),
    treatment = case_when(
      treatment %in% c("Water", "Pooled", "Orig") ~ "Daily replacement",
      treatment == "New" ~ "QOD replacement"
    )
  ) %>%
  arrange(plate, col, well) %>%
  group_by(exp_date, exp_design, time_point, species, treatment, conc) %>%
  mutate(worm_number = row_number())

media_timeline <- tibble(
  exp_design = c(
    "Drug<br>assay",
    "Complete<br>media",
    "RNAi<br>timeline"
  ),
  complete_start = -Inf,
  complete_stop = c(0, 120, 48),
  incomplete_start = c(0, NA, 48),
  incomplete_stop = c(120, NA, 120),
) %>%
  filter(exp_design != "RNAi<br>timeline")

(optimization_plot <- optimization_data %>%
  mutate(
    species = factor(species,
      levels = c("Bpa", "Bma"),
      labels = c("*B. pahangi*", "*B. malayi*")
    ),
    scaled = case_when(
      !exp_date %in% c("20201201", "20210120", "20210126", "20210202", "20210316") & time_point == "120hr" ~ scaled / 2,
      TRUE ~ scaled
    ),
    time_point = case_when(
      time_point == "0hr" ~ 0,
      time_point == "24hr" ~ 24,
      time_point == "48hr" ~ 48,
      time_point == "72hr" ~ 72,
      time_point == "120hr" ~ 120
    )
  ) %>%
  # filter out those that aren't releasing at 0hr
  group_by(exp_date, exp_design, species, treatment, conc, worm_number) %>%
  mutate(
    low_fecundity = case_when(
      any(time_point == 0 & scaled < 200) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(
    # low_fecundity == FALSE,
    time_point != 96,
    exp_design != "RNAi<br>timeline",
    treatment != "QOD replacement"
  ) %>%
  ggplot() +
  geom_rect(
    data = media_timeline, aes(xmin = complete_start, xmax = complete_stop),
    ymin = 2750, ymax = 3000, fill = "#E69F00", alpha = 0.9
  ) +
  geom_rect(
    data = media_timeline, aes(xmin = incomplete_start, xmax = incomplete_stop),
    ymin = 2750, ymax = 3000, fill = "#56B4E9", alpha = 0.9
  ) +
  geom_line(aes(x = time_point, y = scaled, group = interaction(exp_date, worm_number)), alpha = 0.3, size = 0.25) +
  # geom_point(aes(x = time_point, y = scaled), alpha = 0.25, size = 0.3) +
  geom_smooth(aes(x = as.numeric(time_point), y = scaled), size = 1, se = FALSE, color = "black") +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 1) +
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 1) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96, 120)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(500, 1500, 2500), limits = c(0, 3000)) +
  facet_grid(rows = vars(factor(exp_design, levels = c("Complete<br>media", "RNAi<br>timeline", "Drug<br>assay"))), cols = vars(species)) +
  labs(x = "Time point (hr)", y = "Inferred progeny") +
  theme_nw2() +
  theme(
    strip.text.y = element_blank(),
    # panel.grid.major.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_markdown(angle = 0, hjust = 0.5)
  ) +
  NULL
)

save_plot(here("Fig3/subplots/Fig3c.pdf"), optimization_plot, base_width = 4, base_height = 4)
saveRDS(optimization_plot, here("Fig3/subplots/optimized_fecundity.rds"))

# batch differences -------------------------------------------------------

(batch <- data %>%
   filter(time_point == "0hr", species %in% c("Bpa", "Bma"), as.numeric(exp_date) > 20210101) %>%
   mutate(species = factor(species,
                           levels = c("Bpa", "Bma"),
                           labels = c("*B. pahangi*", "*B. malayi*")
   )) %>%
   ggplot(aes(x = scaled, y = exp_date)) +
   geom_density_ridges(aes(color = species, fill = species),
                       alpha = 0.7,
                       jittered_points = TRUE,
                       position = position_points_jitter(width = 0.5, height = 0),
                       point_shape = "|", point_size = 2
   ) +
   scale_x_continuous(expand = c(0, 0)) +
   # scale_y_continuous(expand = expansion(add = c(0, .00001))) +
   scale_color_manual(values = c('grey', 'black')) +
   scale_fill_manual(values = c('grey', 'black')) +
   facet_grid(rows = vars(species), scales = "free") +
   labs(x = "Inferred progeny", y = "Experiment date", fill = "Species", color = "Species") +
   theme_nw2() +
   theme(
     legend.position = "empty"
   ) +
   NULL
)

save_plot(here("Fig3/supplementary/SupplementaryFig3.pdf"), batch, base_width = 6, base_height = 6)
