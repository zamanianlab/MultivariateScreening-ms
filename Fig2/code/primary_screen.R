# data wrangling/plotting
library(tidyverse)
library(janitor)

# stats
library(broom)

# other plotting
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(ZamanianLabThemes)
library(ggrepel)

# misc
library(conflicted)
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

# import ------------------------------------------------------------------

# snake order of the IX recording
record_order <- tibble(
  row = LETTERS[1:8], `01` = seq(1, 8, 1), `02` = seq(16, 9, -1),
  `03` = `01` + 16, `04` = `02` + 16, `05` = `03` + 16, `06` = `04` + 16,
  `07` = `05` + 16, `08` = `06` + 16, `09` = `07` + 16, `10` = `08` + 16,
  `11` = `09` + 16, `12` = `10` + 16
) %>%
  pivot_longer(`01`:`12`, names_to = "col", values_to = "well_number")

# norm_flow is the raw flow divided by the log of the worm area
m_data <- read_rds(here("Fig2/data/motility_data.rds")) %>%
  mutate(assay_date = str_extract(plate, "20[0-9]{6}"), .before = plate) %>%
  mutate(
    norm_flow = optical_flow / log2(worm_area), .after = worm_area,
    plate = str_extract(plate, "p[0-9]{2}")
  ) %>%
  # keep only the 12hr data using the ARPs from the SMSF
  filter(assay_date == "20210819") %>%
  # remove time points for joining
  select(assay_date:conc, strains, optical_flow, worm_area, norm_flow, -other)

v_data <- read_rds(here("Fig2/data/celltox_data.rds")) %>%
  # change the date so the motility and viability dates alitn
  mutate(assay_date = as.character(as.numeric(str_extract(plate, "20[0-9]{6}")) - 1), .before = plate) %>%
  mutate(plate = str_extract(plate, "p[0-9]{2}")) %>%
  # keep only the data using the ARPs from the SMSF
  filter(assay_date == "20210819") %>%
  # recode some incorrect metadata
  mutate(
    treatment = case_when(
      well %in% c("A01", "B01", "C01", "D01", "E12", "F12", "G12", "H12") ~ "DMSO",
      well %in% c("E01", "F01", "G01", "H01", "A12", "B12", "C12", "D12") ~ "HeatKilled",
      TRUE ~ treatment
    ),
    conc = case_when(
      well %in% c("A01", "B01", "C01", "D01", "E12", "F12", "G12", "H12") ~ "1p",
      well %in% c("E01", "F01", "G01", "H01", "A12", "B12", "C12", "D12") ~ "HK",
      TRUE ~ conc
    )
  ) %>%
  # remove time points for joining
  select(assay_date, plate:conc, strains, area_occupied_area_occupied_green_worms, -other) %>%
  group_by(across(!where(is.numeric))) %>%
  summarise(celltox = mean(area_occupied_area_occupied_green_worms)) %>%
  ungroup()

data <- left_join(m_data, v_data) %>%
  left_join(., record_order) %>%
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
  select(assay_date:strains, type, well_number, everything())

# correlation of phenotypes
pheno_corr <- cor(data$norm_flow, data$celltox)

# trim outliers -----------------------------------------------------------

pc_stats_v <- as.numeric(quantile(
  filter(data, treatment == "Positive Control")$celltox,
  c(0, 0.25, 0.5, 0.75, 1)
))
nc_stats_v <- as.numeric(quantile(
  filter(data, treatment == "Negative Control")$celltox,
  c(0, 0.25, 0.5, 0.75, 1)
))

pc_iqr_v <- diff(pc_stats_v[c(2, 4)])
nc_iqr_v <- diff(nc_stats_v[c(2, 4)])

coef <- 3

pc_outliers_v <- pc_stats_v[2] - coef * pc_iqr_v
nc_outliers_v <- nc_stats_v[4] + coef * nc_iqr_v

(v_distribution <- data %>%
    ggplot(aes(x = 1, y = celltox)) +
    geom_quasirandom(aes(color = type), alpha = 0.25, size = 3) +
    geom_boxplot(
      data = . %>% filter(str_detect(type, "Control")),
      aes(color = type),
      x = 0, position = position_dodge(0)
    ) +
    geom_hline(yintercept = nc_outliers_v, color = "#FAD510", linetype = "dashed") +
    geom_hline(yintercept = pc_outliers_v, color = "#CB2314", linetype = "dashed") +
    scale_color_manual(values = c("#FAD510", "#CB2314", "#273046")) +
    theme_nw2() +
    NULL)

pc_stats_m <- as.numeric(quantile(
  filter(data, treatment == "Positive Control")$norm_flow,
  c(0, 0.25, 0.5, 0.75, 1)
))
nc_stats_m <- as.numeric(quantile(
  filter(data, treatment == "Negative Control")$norm_flow,
  c(0, 0.25, 0.5, 0.75, 1)
))

pc_iqr_m <- diff(pc_stats_m[c(2, 4)])
nc_iqr_m <- diff(nc_stats_m[c(2, 4)])

pc_outliers_m <- pc_stats_m[4] + coef * pc_iqr_m
nc_outliers_m <- nc_stats_m[2] - coef * nc_iqr_m

(m_distribution <- data %>%
    ggplot(aes(x = 1, y = norm_flow)) +
    geom_quasirandom(aes(color = type), alpha = 0.25, size = 3) +
    geom_boxplot(
      data = . %>% filter(str_detect(type, "Control")),
      aes(color = type),
      x = 0, position = position_dodge(0), coef = 2.5
    ) +
    geom_hline(yintercept = nc_outliers_m, color = "#FAD510", linetype = "dashed") +
    geom_hline(yintercept = pc_outliers_m, color = "#CB2314", linetype = "dashed") +
    scale_color_manual(values = c("#FAD510", "#CB2314", "#273046")) +
    theme_nw2() +
    NULL)

# remove 7 outliers in the controls
trimmed <- data %>%
  # viability trim
  filter(case_when(
    type == "Negative Control" & celltox < nc_outliers_v ~ TRUE,
    type == "Positive Control" & celltox > pc_outliers_v ~ TRUE,
    type == "Treated" & is.numeric(celltox) ~ TRUE
  )) %>%
  # motility trim
  filter(case_when(
    type == "Negative Control" & norm_flow > nc_outliers_m ~ TRUE,
    type == "Positive Control" & norm_flow < pc_outliers_m ~ TRUE,
    type == "Treated" & is.numeric(norm_flow) ~ TRUE
  ))

# correct for drift and normalize -----------------------------------------

# show the drift in flow
drift <- trimmed %>%
  select(assay_date, plate, row, col, well, treatment, type, well_number, optical_flow, norm_flow, celltox) %>%
  arrange(well_number) %>%
  # filter(str_detect(treatment, 'Control')) %>%
  mutate(smooth_row = as.numeric(factor(row))) %>%
  filter(!col %in% c("01", "12")) %>%
  ggplot(aes(x = smooth_row, y = optical_flow)) +
  geom_point(size = 1, alpha = 0.75) +
  geom_smooth(method = "lm", color = "indianred") +
  scale_x_continuous(labels = c("A", "C", "E", "G"), breaks = seq(1, 8, 2)) +
  facet_grid(cols = vars(col)) +
  labs(x = "Row", y = "Raw motility") +
  theme_nw2() +
  theme(legend.position = "empty") +
  NULL

c_model <- trimmed %>%
  mutate(smooth_row = as.numeric(factor(row))) %>%
  mutate(smooth_col = as.numeric(col)) %>%
  filter(str_detect(treatment, "Control")) %>%
  lm(norm_flow ~ smooth_row + smooth_col, data = .) %>%
  tidy()

row_coef <- as.numeric(c_model[2, 2])
col_coef <- as.numeric(c_model[3, 2])

corrected <- trimmed %>%
  mutate(smooth_row = as.numeric(factor(row))) %>%
  mutate(smooth_col = as.numeric(col)) %>%
  mutate(corrected_flow = norm_flow - smooth_row * row_coef - smooth_col * col_coef)

(corrected_drift <- corrected %>%
    select(assay_date, plate, row, col, well, treatment, type, well_number, corrected_flow) %>%
    arrange(well_number) %>%
    mutate(smooth_row = as.numeric(factor(row))) %>%
    filter(!col %in% c("01", "12")) %>%
    ggplot(aes(x = smooth_row, y = corrected_flow)) +
    geom_point(size = 1, alpha = 0.75, color = "#273046") +
    geom_smooth(method = "lm", color = "indianred") +
    scale_x_continuous(labels = c("A", "C", "E", "G"), breaks = seq(1, 8, 2)) +
    facet_grid(cols = vars(col)) +
    labs(x = "Row", y = "Corrected motility") +
    theme_nw2() +
    theme(legend.position = "empty") +
    NULL)

(corrected_drift_total <- corrected %>%
    ggplot(aes(x = well_number, y = corrected_flow, color = type)) +
    geom_point(aes(color = type), alpha = 0.75, size = 1.35) +
    scale_color_manual(values = c("#FAD510", "#CB2314", "#273046")) +
    scale_x_continuous(limits = c(1, 96), breaks = seq(1, 96, 8)) +
    labs(x = "Well number", y = "Regressed motility") +
    theme_nw2() +
    theme(legend.position = "empty",
          axis.text.x = element_markdown(angle = 0, hjust = 0.5)) +
    NULL)

write_rds(corrected_drift_total, here("Fig2/subplots/corrected_drift.rds"))
save_plot(here('Fig2/subplots/Fig2d.pdf'), corrected_drift_total,
          base_height = 3, base_width = 5)

(correct_dist <- corrected %>%
    ggplot() +
    geom_histogram(aes(x = corrected_flow, fill = type)) +
    scale_fill_manual(values = c("#FAD510", "#CB2314", "#273046")) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Normalized, regressed motility", y = "", color = "", fill = "") +
    coord_flip() +
    theme_nw2() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "empty"
    ) +
    NULL)

write_rds(correct_dist, here("Fig2/subplots/corrected_dist.rds"))
save_plot(here('Fig2/subplots/Fig2d_inset.pdf'), correct_dist,
          base_height = 2, base_width = 2.5)

# rescale -----------------------------------------------------------------

control_summary <- corrected %>%
  filter(str_detect(treatment, "Control") == TRUE) %>%
  pivot_longer(c(corrected_flow, celltox), names_to = "assay_type", values_to = "value") %>%
  group_by(assay_date, plate, assay_type, treatment) %>%
  summarise(control_mean = mean(value)) %>%
  pivot_wider(names_from = treatment, values_from = control_mean) %>%
  rename(n_mean = 4, p_mean = 5)

rescaled <- corrected %>%
  pivot_longer(c(corrected_flow, celltox), names_to = "assay_type", values_to = "value") %>%
  left_join(., control_summary) %>%
  mutate(rescaled = (value - n_mean) / (p_mean - n_mean)) %>%
  pivot_wider(id_cols = c(assay_date:conc, strains, type), names_from = assay_type, values_from = rescaled) %>%
  rename(motility = corrected_flow, viability = celltox)

write_rds(rescaled, here('Fig2/data/rescaled_data.rds'))

# z plot ------------------------------------------------------------------

plate_summary <- rescaled %>%
  group_by(assay_date, plate) %>%
  summarise(
    mean_motility = mean(motility), sd_motility = sd(motility),
    mean_viability = mean(viability), sd_viability = sd(viability)
  )

z <- rescaled %>%
  select(assay_date:conc, type, motility, viability) %>%
  left_join(., plate_summary) %>%
  mutate(
    z_motility = (motility - mean_motility) / sd_motility,
    z_viability = (viability - mean_viability) / sd_viability
  )

motility_cutoff <- 1
viability_cutoff <- 1

hits <- z %>%
  mutate(hit_type = case_when(
    z_motility > motility_cutoff & z_viability > viability_cutoff & !treatment == "Positive Control" ~ "Both",
    z_motility > motility_cutoff & z_viability < viability_cutoff & !treatment == "Positive Control" ~ "Motility",
    z_motility < motility_cutoff & z_viability > viability_cutoff & !treatment == "Positive Control" ~ "Viability"
  )) %>%
  filter(!is.na(hit_type)) %>%
  mutate(z_total = z_motility + z_viability) %>%
  # filter(z_total > 1) %>%
  arrange(desc(z_total)) %>%
  mutate(number = seq(1, n(), 1))

# correlation of phenotypes in hits
hits_corr <- data %>% 
  filter(treatment %in% hits$treatment)

hits_corr <- cor(hits_corr$norm_flow, hits_corr$celltox)

hit_type <- hits %>%
  group_by(assay_date, hit_type) %>%
  summarize(n = n()) %>%
  mutate(
    x = c(3.25, 3.25, -1),
    y = c(5, 0.8, 5)
  )

(z_plot <- z %>%
    left_join(., select(hits, treatment, hit_type)) %>%
    mutate(
      hit_type = case_when(
        is.na(hit_type) == TRUE ~ "Non-Hit",
        TRUE ~ hit_type
      ),
      type = case_when(
        hit_type != "Non-Hit" ~ "Hit",
        TRUE ~ type
      )
    ) %>%
    ggplot() +
    annotate("rect",
             xmin = motility_cutoff, xmax = Inf, ymin = -Inf, ymax = Inf,
             alpha = .2, color = NA, fill = "#354823"
    ) +
    annotate("rect",
             xmin = -Inf, xmax = Inf, ymin = viability_cutoff, ymax = Inf,
             alpha = .2, color = NA, fill = "#354823"
    ) +
    geom_point(aes(x = z_motility, y = z_viability, color = type, shape = hit_type),
               size = 1.35
    ) +
    geom_rug(aes(x = z_motility, y = z_viability, color = type)) +
    geom_text(data = hit_type, aes(x = x, y = y, label = n), fontface = "italic") +
    geom_text_repel(
      data = hits, aes(x = z_motility, y = z_viability, label = number),
      color = "#00AFB5", size = 2
    ) +
    scale_color_manual(values = c("#00AFB5", "#FAD510", "#CB2314", "#273046")) +
    scale_shape_manual(values = c(18, 15, 16, 17)) +
    labs(x = "Motility Z-score", y = "Viability Z-score", shape = "Hit Type", color = "Compound", fill = "Compound") +
    theme_nw2() +
    theme(
      legend.position = "right",
      axis.text.x = element_markdown(angle = 0, hjust = 0.5)
    ) +
    NULL)

write_rds(z_plot, here('Fig2/subplots/primary_screen.rds'))
save_plot(here('Fig2/subplots/Fig2g.pdf'), z_plot,
          base_height = 4, base_width = 5)


# 100 uM screen -----------------------------------------------------------

hundred_data <- readRDS("data/motility_data.rds") %>%
  mutate(assay_date = str_extract(plate, "20[0-9]{6}"), .before = plate) %>%
  mutate(
    norm_flow = optical_flow / log2(worm_area), .after = worm_area,
    plate = str_extract(plate, "p[0-9]{2}")
  ) %>%
  filter(assay_date == "20201118") %>%
  select(-other)

hundred_summary <- hundred_data %>%
  filter(!treatment %in% c("DMSO", "SodiumAzide")) %>%
  group_by(assay_date, plate) %>%
  summarise(mean_flow = mean(norm_flow), sd_flow = sd(norm_flow))

z_hundred <- hundred_data %>%
  select(assay_date:conc, optical_flow, norm_flow) %>%
  left_join(., hundred_summary) %>%
  mutate(z_flow = (norm_flow - mean_flow) / sd_flow)

(hundred_plot <- z_hundred %>%
    mutate(
      treatment = case_when(
        treatment == "SodiumAzide" ~ "Positive Control",
        treatment == "DMSO" ~ "Negative Control",
        TRUE ~ treatment
      ),
      type = case_when(
        str_detect(treatment, "Control") == TRUE ~ treatment,
        TRUE ~ "Treated"
      )
    ) %>%
    ggplot() +
    geom_histogram(aes(x = optical_flow, fill = type)) +
    scale_fill_manual(values = c("#FAD510", "#CB2314", "#273046")) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Raw motility", y = "", color = "", fill = "") +
    coord_flip() +
    theme_nw2() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "empty"
    ) +
    NULL)

write_rds(hundred_plot, "subplots/hundred_screen.rds")
save_plot(here('Fig2/subplots/Fig2a_inset.pdf'), hundred_plot,
          base_height = 2, base_width = 2.5)

# Hits table --------------------------------------------------------------

compound_info <- read_csv("metadata/compound_information_mini.csv") %>%
  clean_names() %>%
  mutate(plate_number = str_c("p", str_pad(plate_number, 2, "left", "0"))) %>%
  drop_na(cat_no)

compound_summary <- compound_info %>%
  group_by(target_class) %>%
  summarise(compound_sum = n())

hits_table <- hits %>%
  left_join(., compound_info, by = c("treatment" = "compound_name", "plate" = "plate_number", "well" = "plate_location")) %>%
  select(number, treatment, target_class, primary_target, primary_action, contains("z"), hit_type)

write_csv(hits_table, "tables/Table1.csv")



# Hits plot ---------------------------------------------------------------

(hit_bar <- hits_table %>%
   select(-hit_type) %>%
   ungroup() %>%
   group_by(target_class) %>%
   summarize(target_sum = n()) %>%
   left_join(compound_summary) %>%
   mutate(
     target_class = case_when(
       target_class == "Enzyme-Linked Receptors" ~ "Enzyme-Linked<br>Receptors",
       TRUE ~ target_class
     ),
     prop = target_sum / compound_sum
   ) %>%
   ggplot(aes(x = reorder(target_class, target_sum), y = prop)) +
   geom_label(aes(label = target_sum),
              fill = "grey11", color = "white", size = 2
   ) +
   labs(x = "Target class", y = "Percent of class") +
   scale_y_continuous(labels = scales::percent, limits = c(0.009, 0.045)) +
   coord_flip() +
   theme_nw2() +
   theme(
     axis.text.x = element_markdown(angle = 0, hjust = 0.5)
   ) +
   NULL)

write_rds(hit_bar, here("Fig2/subplots/hit_bar.rds"))
save_plot(here('Fig2/subplots/Fig2h.pdf'), hit_bar,
          base_height = 2, base_width = 3)

# export data -------------------------------------------------------------

# hits %>%
#   select(number, treatment, z_motility, z_viability) %>%
#   mutate(treatment = str_remove(treatment, " hydrochloride| maleate| citrate")) %>%
#   write_rds("data/z_scores.rds")
