library(tidyverse)
library(cowplot)
library(magick)
library(conflicted)
library(patchwork)
library(here)

conflict_prefer("area", "patchwork")

# Figure 2 ----------------------------------------------------------------

### A

a_hundred_drift <- read_rds(here('Fig2/subplots/uncorrected_drift_hundred.rds')) +
  geom_vline(data = tibble(xintercept = c(1, 17, 33, 49, 65, 81)),
             aes(xintercept = xintercept), 
             linetype = 'dashed',
             color = 'grey') +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())

a_inset <- read_rds(here('Fig2/subplots/hundred_screen.rds')) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        plot.background = element_rect(fill = NA, colour = NA))

a <- a_hundred_drift + 
  inset_element(a_inset,
                left = 0.7, bottom = 0.5, right = 1, top = 1, align_to = 'full',
                ignore_tag = TRUE)

### B
b_frames <- read_rds(here('Fig2/subplots/frame_optimization.rds'))

### C
c_celltox <- image_read_pdf((here('Fig2/subplots/Fig2c.pdf')))
c_celltox <- ggdraw() + draw_image(c_celltox)  #+ theme(plot.margin = margin(-1, -1, -1, -1, "cm"))

### D
d_corrected <- read_rds(here('Fig2/subplots/corrected_drift.rds')) +
  geom_vline(data = tibble(xintercept = c(1, 17, 33, 49, 65, 81)),
             aes(xintercept = xintercept), 
             linetype = 'dashed',
             color = 'grey') +
  theme(
    panel.grid.major.x = element_blank())

d_inset <- read_rds(here('Fig2/subplots/corrected_dist.rds')) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        plot.background = element_rect(fill = NA, colour = NA))

d <- d_corrected + 
  inset_element(d_inset,
               left = 0.6, bottom = 0.15, right = 0.9, top = 0.65, align_to = 'full',
               ignore_tag = TRUE)
### E
e_z_factor <- read_rds(here('Fig2/subplots/primary_screen_z_factor.rds')) +
  scale_color_manual(values = c("#FAD510", "#CB2314"), guide = 'none') +
  scale_fill_manual(values = c("#FAD510", "#CB2314"), guide = 'none') +
  theme(legend.position = 'empty')


### F
f_design <- image_read_pdf((here('Fig2/subplots/Fig2f.pdf')))
f <- ggdraw() + draw_image(f_design)

### G
g_screen <- read_rds(here('Fig2/subplots/primary_screen.rds'))

### H
h_hit_bar <- read_rds(here('Fig2/subplots/hit_bar.rds'))

### I
i_ec <- read_rds(here('Fig2/subplots/ec50_bar.rds'))

# set layout

abc_layout <- c(
  area(t = 1, l = 1, b = 8, r = 12),   # a
  area(t = 1, l = 13, b = 3, r = 20),  # b
  area(t = 4, l = 13, b = 8, r = 20)   # c
  
)

abc_patch <- 
    a + 
    b_frames + 
    c_celltox +
    plot_layout(design = abc_layout)

de_layout <- c(
  area(t = 1, l = 1, b = 4, r = 12),   # d
  area(t = 1, l = 13, b = 4, r = 20),  # e
  area(t = 5, l = 1, b = 5.1, r = 20)
)

de_patch <-
    d + e_z_factor + get_legend(a_hundred_drift + theme(legend.position = 'bottom',
                                                        legend.title = element_blank())) +
    plot_layout(design = de_layout)

# abc_patch / de_patch / f

ghi_layout <- c(
  area(t = 1, l = 1, b = 4, r = 14),   # a
  area(t = 1, l = 15, b = 2, r = 20),  # b
  area(t = 3, l = 15, b = 4, r = 20)   # c
)

ghi_patch <-
  g_screen +
  h_hit_bar + 
  i_ec +
  plot_layout(design = ghi_layout)

fig2 <- abc_patch / de_patch / f / ghi_patch +
  plot_layout(heights = c(0.6, 0.7, 1, 1)) +
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(size = 12, face = 'bold')
  )

save_plot(here('Fig2/Fig2.pdf'),
          fig2, base_width = 8.66, base_height = 12)
save_plot(here('Fig2/Fig2.png'),
          fig2, base_width = 8.66, base_height = 12)


# Supplementary Figure 2 --------------------------------------------------

sf_a <- read_rds(here('Fig2/supplementary/SupplementaryFig3a.rds'))
sf_b <- read_rds(here('Fig2/supplementary/SupplementaryFig3b.rds'))

sf2 <- plot_grid(sf_a, sf_b, nrow = 2, labels = 'auto')

save_plot(here('Fig2', 'supplementary', 'SupplementaryFig3.pdf'),
          sf2, base_width = 6.5, base_height = 8)
