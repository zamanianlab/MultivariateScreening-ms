library(tidyverse)
library(ggtext)
library(cowplot)
library(magick)
library(patchwork)
library(here)

conflict_prefer("filter", "dplyr")


# Figure 6 ----------------------------------------------------------------

cutoff <- -1

a_screen <- read_rds(here('Fig6/subplots/cel_screen_dev.rds')) +
  labs(color = '') +
  theme(legend.position = 'bottom',
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

b_followup <- read_rds(here('Fig6/subplots/cel_followup.rds'))  +
  scale_fill_manual(values = c("#FAD510", "#CB2314", '#00AFB5'), guide = 'none') 

drug_label <- tibble(x = c(3.5, 1.5),
                     y = 1.7,
                     label = c('ML 228', 'NAV 2729'),
                     panel = c('ML 228', 'NAV 2729'))

c_deconvolution <- read_rds(here('Fig6/subplots/deconvolution.rds')) +
  geom_text(data = drug_label, aes(x = x, y = y, label = label)) +
  labs(fill = '') +
  theme(
    strip.text.x = element_blank(),
    # axis.text.x = element_markdown(angle = 0, hjust = 0.5),
    legend.position = 'right'
  )

hypoxia <- image_read_pdf(here('Fig6/subplots/hypoxia_reporter.pdf'))

d_hypoxia <- ggdraw() +
  draw_image(hypoxia)

e_iron <- read_rds(here('Fig6/subplots/k_s_media.rds'))


# cowplot -----------------------------------------------------------------

legend <- get_legend(a_screen)

(top <- plot_grid(a_screen, b_followup, axis = 't', align = 'h', rel_widths = c(0.85, 1),
                 labels = c('a', 'b'), label_size = 12))

bottom <- plot_grid(d_hypoxia, e_iron, axis = 't', align = 'h', rel_widths = c(0.65, 1),
                    labels = c('d', 'e'), label_x = c(0, -0.03), label_size = 12)

fig6 <- plot_grid(top, c_deconvolution, bottom, nrow = 3, labels = c('', 'c', ''),
                  rel_heights = c(1.15, 0.35, 0.5), label_y = 1.15, label_size = 12)

save_plot(here('Fig6/Fig6.pdf'),
          fig6, base_width = 8, base_height = 8)

save_plot(here('Fig6/Fig6.png'),
          fig6, base_width = 8, base_height = 8)



