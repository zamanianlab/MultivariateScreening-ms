library(tidyverse)
library(ggtext)
library(cowplot)
library(magick)
library(patchwork)
library(here)


# Figure 4 ----------------------------------------------------------------

a_motility <- read_rds(here('Fig4/subplots/motility.rds')) +
  theme(
    axis.title.x = element_blank(),
    legend.position = 'empty'
    )

b_distinct_lines <- read_rds(here('Fig4/subplots/distinct_lines.rds'))

c_fecundity <- read_rds(here('Fig4/subplots/proportions.rds')) +
  theme(
    legend.position = 'empty'
  )

c_inset <- read_rds(here('Fig4/subplots/baseline.rds')) +
  theme(
    axis.text.y =  element_markdown(size = 7),
    axis.title.y = element_markdown(size = 9),
    # axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = 'white', colour = 'black', size = 0.25),
    plot.background = element_rect(fill = NA, colour = NA)
    )

c <- c_fecundity + 
  inset_element(c_inset,
                left = 0.75, bottom = 0.5, right = 1, top = 1, align_to = 'full', ignore_tag = TRUE)

d_lactate <- read_rds(here('Fig4/subplots/lactate.rds')) +
  theme(
    legend.position = 'bottom'
  )

e_viability <- image_read_pdf(here('Fig4/subplots/Fig4e.pdf'))
e <- ggdraw() + draw_image(e_viability)

legend <- get_legend(d_lactate)


# cowplot -----------------------------------------------------------------

left <- plot_grid(a_motility, c, d_lactate + theme(legend.position = 'empty'), legend,
                 nrow = 4, align = 'v', axis = 'l', rel_heights = c(1, 0.5, 1, 0.1),
                 labels = c('a', 'c', 'd'), label_size = 12)

right <- plot_grid(b_distinct_lines, e,
                   # scale = c(1, 1.05),
                   nrow = 2, rel_heights = c(1, 0.5),
                   # align = 'v', axis = 'lr',
                   labels = c('b', 'e'), label_size = 12)

fig4 <- plot_grid(left, right,
                  ncol = 2, rel_widths = c(1, 0.5))

save_plot(here('Fig4/Fig4.pdf'),
          fig4, base_width = 8.66, base_height = 11)
save_plot(here('Fig4/Fig4.png'),
          fig4, base_width = 8.66, base_height = 11)