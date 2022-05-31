library(tidyverse)
library(cowplot)
library(patchwork)
library(magick)
library(here)


# Figure 3 ----------------------------------------------------------------

a_motility <- read_rds(here('Fig3/subplots/optimized_motility.rds')) +
  theme(
    # axis.title.x = element_blank(),
    legend.position = 'empty')

b_metabolism <- read_rds(here('Fig3/subplots/heatkilled_metabolism.rds'))

c_fecundity <- read_rds(here('Fig3/subplots/optimized_fecundity.rds')) +
  theme(legend.position = 'empty')

d_fecundity <- read_rds(here('Fig3/subplots/baseline_fecundity.rds')) +
  theme(legend.position = 'right')

e_fecundity <- read_rds(here('Fig3/subplots/heatkilled_fecundity.rds'))

f_scheme <- image_read_pdf(here('Fig3/subplots/Fig3f.pdf'))
f <- ggdraw() + draw_image(f_scheme)

b_align <- cowplot::align_plots(c_fecundity, e_fecundity, 
                                align = 'h', axis = 'b')

left <- plot_grid(a_motility, b_align[[1]], 
                  align = 'v', axis = 'lr', rel_heights = c(1, 1.1),
                  nrow = 2, labels = c('a', 'c'), label_size = 12)

right <- plot_grid(b_metabolism, d_fecundity, b_align[[2]],
                   align = 'v', axis = 'l',
                   nrow = 3, labels = c('b', 'd', 'e'), label_size = 12, hjust = 0)

top <- plot_grid(left, right,
                 rel_widths = c(1.25, 1), ncol = 2)

fig3 <- plot_grid(
  top, f, rel_heights = c(2.7, 1), nrow = 2, labels = c('', 'f'), label_size = 12)

save_plot(here('Fig3/Fig3.pdf'), 
          fig3, base_height = 7, base_width = 7)
save_plot(here('Fig3/Fig3.png'), 
          fig3, base_height = 7, base_width = 7)
