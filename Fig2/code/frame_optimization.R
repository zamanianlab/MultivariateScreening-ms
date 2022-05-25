library(tidyverse)
library(cowplot)
library(ggtext)
library(corrr)
library(wesanderson)
library(ZamanianLabThemes)
library(here)

# frame correlation -------------------------------------------------------

# snake order of the IX recording
record_order <- tibble(
  row = LETTERS[1:8], `01` = seq(1, 8, 1), `02` = seq(16, 9, -1),
  `03` = `01` + 16, `04` = `02` + 16, `05` = `03` + 16, `06` = `04` + 16,
  `07` = `05` + 16, `08` = `06` + 16, `09` = `07` + 16, `10` = `08` + 16,
  `11` = `09` + 16, `12` = `10` + 16
) %>%
  pivot_longer(`01`:`12`, names_to = "col", values_to = "well_number")

frame_data <- read_csv(here("Fig2/data/frames.csv")) %>%
  pivot_longer(cols = 4:11, names_to = "plate_frame", values_to = "optical_flow") %>%
  separate(plate_frame, c("plate", "n_frames")) %>%
  mutate(col = as.character(col))

frame_cor <- frame_data %>%
  pivot_wider(well:plate, names_from = n_frames, values_from = optical_flow) %>%
  group_nest(plate) %>%
  mutate(
    cor_matrix = map(data, ~ correlate(select(.x, 4:7))),
    cor_matrix = map(cor_matrix, shave),
    cor_matrix = map(cor_matrix, stretch)
  ) %>%
  unnest(cor_matrix) %>%
  select(-data) %>%
  filter(!is.na(r), x == 40) %>%
  unite(comparison, x:y, remove = FALSE) %>%
  rename(n_frames = y) %>%
  select(-x) %>%
  filter(plate == 172) %>%
  mutate(
    comparison = str_replace(comparison, "_", " vs "),
    r = str_c("*", round(r, 3), "*"),
    y = c(7.75, 7.05, 6.35),
    n_frames = factor(n_frames, levels = c("5", "10", "20", "40"))
  )

label <- tibble(
  n_frames = c("40", "20", "10", "5"),
  x = 0,
  y = c(104514356, 49801758, 23594987, 10426864)
)

# plot the data and correlation coefficients
(cor_plot <- frame_data %>%
    filter(plate == 172) %>%
    mutate(
      col = str_pad(col, 2, pad = "0"),
      n_frames = factor(n_frames, levels = c("5", "10", "20", "40"))
    ) %>%
    left_join(record_order) %>%
    ggplot() +
    geom_line(aes(x = well_number, y = log10(optical_flow), color = n_frames)) +
    geom_richtext(
      data = label,
      aes(label = n_frames, x = x, y = log10(y), color = n_frames),
      fill = NA, label.color = NA, size = 2.5, label.padding = unit(c(0, 0, 0, 0), "lines"),
      hjust = 0.7, show.legend = FALSE
    ) +
    scale_color_manual(values = wes_palette("Moonrise2"), guide = "none") +
    scale_x_continuous(limits = c(0, 96), breaks = seq(1, 96, 8)) +
    labs(x = "Well number", y = "log<sub>10</sub>Motility", color = "Number of frames") +
    theme_nw2() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "empty",
      axis.text.x = element_markdown(angle = 0, hjust = 0.5)
    ) +
    NULL)

save_plot(here("Fig2/subplots/Fig2b.pdf"), cor_plot, base_width = 4, base_height = 1)
write_rds(cor_plot, here("Fig2/subplots/frame_optimization.rds"))

# calculate autocorrelation  ----------------------------------------------

m_data <- readRDS("data/motility_data.rds") %>%
  mutate(norm_flow = optical_flow / log2(worm_area), .after = worm_area) %>%
  mutate(assay_date = str_extract(plate, "20[0-9]{6}"), .before = plate) %>%
  filter(assay_date %in% c("20201118", "20210819")) %>%
  mutate(
    treatment = case_when(
      treatment %in% c("HeatKilled", "SodiumAzide") ~ "Positive Control",
      treatment == "DMSO" ~ "Negative Control",
      TRUE ~ treatment
    )
  ) %>%
  mutate(
    type = case_when(
      !treatment %in% c("Negative Control", "Positive Control") ~ "Treated",
      TRUE ~ treatment
    )
  ) %>%
  left_join(record_order)

# show the drift in flow
(drift <- m_data %>%
  select(assay_date, plate, row, col, well, treatment, type, norm_flow) %>%
  mutate(smooth_row = as.numeric(factor(row))) %>%
  filter(
    assay_date == "20201118",
    !col %in% c("01", "12")
  ) %>%
  ggplot(aes(x = smooth_row, y = norm_flow)) +
  geom_point(size = 1, alpha = 0.75, color = "#273046") +
  geom_smooth(data = . %>% filter(norm_flow > 4e6), method = "lm", color = "indianred") +
  geom_hline(yintercept = 4e6, linetype = "dashed") +
  scale_x_continuous(labels = c("A", "C", "E", "G"), breaks = seq(1, 8, 2)) +
  facet_grid(cols = vars(col)) +
  labs(x = "Row", y = "Normalized flow") +
  theme_nw2() +
  theme(legend.position = "empty") +
  NULL)

(uncorrected_plot <- m_data %>%
  filter(
    assay_date == "20201118",
    (type == "Positive Control" & well_number > 80) |
      (type == "Negative Control" & well_number < 10) |
      type == "Treated"
  ) %>%
  ggplot(aes(x = well_number, y = optical_flow)) +
  geom_point(aes(color = type), alpha = 0.75, size = 1.35) +
  # geom_smooth(data = . %>% filter(norm_flow > 4e6),
  #             color = 'black', method = "lm", formula = y ~ log(x)) +
  labs(x = "Well number", y = "Raw motility") +
  scale_x_continuous(limits = c(1, 96), breaks = seq(1, 96, 8)) +
  scale_color_manual(values = c("#FAD510", "#CB2314", "#273046")) +
  theme_nw2() +
  theme(legend.position = "empty",
        axis.text.x = element_markdown(angle = 0, hjust = 0.5)) +
  NULL)

save_plot(here("Fig2/subplots/Fig2a.pdf"), uncorrected_plot, base_width = 5, base_height = 3)
write_rds(uncorrected_plot, here("Fig2/subplots/uncorrected_drift_hundred.rds"))
