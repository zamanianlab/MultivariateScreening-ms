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


# analysis ----------------------------------------------------------------

all_data <- read_rds(here('Fig3/data/metabolism_data.rds'))

(dilution <- all_data %>%
  filter(status != "HK") %>%
  ggplot(aes(x = dilution, y = rfu, color = treatment)) +
  geom_point() +
  facet_wrap(facets = vars(stages, treatment), scales = "free_y", nrow = 2) +
  theme_nw() +
  theme(
    legend.position = "none"
  ) +
  NULL
)

(hk <- all_data %>%
  filter(
    treatment == "Water",
    (stages == "AM" & dilution == 50) |
      (stages == "AF" & dilution == 150)
  ) %>%
  mutate(
    stages = case_when(
      stages == "AF" ~ "Female",
      stages == "AM" ~ "Male"
    ),
    status = case_when(
      status == "A" ~ "Alive",
      status == "HK" ~ "Heat-<br>killed"
    )
  ) %>%
  ggplot(aes(x = status, y = rfu)) +
  geom_errorbar(
    position = position_dodge(width = 1),
    stat = "summary", width = 0.5, alpha = 0.75,
    color = "black"
  ) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean", fill = "black") +
  geom_quasirandom(
    shape = 21, fill = "white", size = 2,
    alpha = 0.75, dodge.width = 0.9, width = 0.1,
    show.legend = FALSE, color = "black"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "RFU") +
  facet_wrap(facets = vars(stages)) +
  theme_nw2() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_markdown(angle = 0, hjust = 0.5)
  ) +
  NULL)

save_plot(here("Fig3/subplots/Fig3b.pdf"), hk, base_width = 3, base_height = 2)
write_rds(hk, here("Fig3/subplots/heatkilled_metabolism.rds"))
