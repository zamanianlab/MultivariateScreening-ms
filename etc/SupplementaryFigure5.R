library(tidyverse)
library(cowplot)
library(here)
library(ggtree)
library(treeio)
library(tidytree)



# arf6 --------------------------------------------------------------------

arf6_t <- read.iqtree(here('etc', 'data', 'trees', 'ARF6_output', 'alignments', 'ARF6_HUMAN.combined_final.aln.contree'))
arf6_t %>% as_tibble() %>% View()

arf6_gg <- arf6_t %>% 
  groupClade(785) %>% 
  ggtree(
    aes(color = group),
    size = 0.2
    # layout = 'equal_angle'
    ) +
  # geom_nodelab(,
  #   size = 1) +
  geom_nodepoint(
    aes(fill = UFboot),
    shape = 21,
    size = 0.5,
    color = 'white'
  ) +
  geom_tiplab(size = 0.5) +
  labs(title = 'ARF6') +
  scale_color_manual(values = c('black', 'indianred'),
                     guide = 'none') +
  scale_fill_viridis_c() +
  NULL

save_plot(here('etc', 'plots', 'arf6_tree.pdf'),
          arf6_gg, base_width = 15, base_height = 30)


# 5HT2A -------------------------------------------------------------------

ht2a_t <- read.iqtree(here('etc', 'data', 'trees', '5HT2A_output', 'alignments', '5HT2A_HUMAN.combined_final.aln.contree'))
ht2a_t %>% as_tibble() %>% View()

ht2a_gg <- ht2a_t %>% 
  groupClade(169) %>% 
  ggtree(
    aes(color = group),
    size = 0.2,
    # layout = 'daylight', branch.length = 'none'
  ) +
  # geom_nodelab(
  #   aes(label = node),
  #   size = 1) +
  geom_nodepoint(
    aes(fill = UFboot),
    shape = 21,
    size = 2,
    color = 'white'
  ) +
  geom_tiplab(size = 1) +
  labs(title = '5HT2A') +
  scale_color_manual(values = c('black', 'indianred'),
                     guide = 'none') +
  scale_fill_viridis_c() +
  NULL

save_plot(here('etc', 'plots', '5ht2a_tree.pdf'),
          ht2a_gg, base_width = 5, base_height = 10)


# AA3R --------------------------------------------------------------------

aa3r_t <- read.iqtree(here('etc', 'data', 'trees', 'AA3R_output', 'alignments', 'AA3R_HUMAN.combined_final.aln.contree'))
aa3r_t %>% as_tibble() %>% View()

aa3r_gg <- aa3r_t %>% 
  groupClade(152) %>% 
  ggtree(
    aes(color = group),
    size = 0.2,
    # branch.length = 'none'
    # layout = 'equal_angle'
  ) +
  # geom_nodelab(
  #   aes(label = node),
  #   size = 1) +
  geom_nodepoint(
    aes(fill = UFboot),
    shape = 21,
    size = 2,
    color = 'white'
  ) +
  geom_tiplab(size = 1) +
  labs(title = 'AA3R') +
  scale_color_manual(values = c('black', 'indianred'),
                     guide = 'none') +
  scale_fill_viridis_c() +
  coord_cartesian(clip = "off") +
  NULL

save_plot(here('etc', 'plots', 'aa3r_tree.pdf'),
          aa3r_gg, base_width = 6, base_height = 8)

# KCNH2 --------------------------------------------------------------------

kcnh2_t <- read.iqtree(here('etc', 'data', 'trees', 'KCNH2_output', 'alignments', 'KCNH2_HUMAN.combined_final.aln.contree')) %>% 
  root(node = 59)
kcnh2_t %>% as_tibble() %>% View()

kcnh2_gg <- kcnh2_t %>% 
  groupClade(42) %>% 
  ggtree(
    aes(color = group),
    size = 0.2,
    # branch.length = 'none'
    # layout = 'circular'
  ) +
  # geom_nodelab(
  #   aes(label = node),
  #   size = 1) +
  geom_nodepoint(
    aes(fill = UFboot),
    shape = 21,
    size = 2,
    color = 'white'
  ) +
  geom_tiplab(size = 1) +
  labs(title = 'KCNH2') +
  scale_color_manual(values = c('black', 'indianred'),
                     guide = 'none') +
  scale_fill_viridis_c() +
  coord_cartesian(clip = "off") +
  NULL

save_plot(here('etc', 'plots', 'kcnh2_tree.pdf'),
          kcnh2_gg, base_width = 6, base_height = 6)


# OPRM --------------------------------------------------------------------

oprm_t <- read.iqtree(here('etc', 'data', 'trees', 'OPRM_output', 'alignments', 'OPRM_HUMAN.combined_final.aln.contree')) %>% 
  root(node = 166)
oprm_t %>% as_tibble() %>% View()

oprm_gg <- oprm_t %>% 
  groupClade(153) %>% 
  ggtree(
    aes(color = group),
    size = 0.2,
    # branch.length = 'none'
    # layout = 'circular'
  ) +
  # geom_nodelab(
  #   aes(label = node),
  #   size = 1) + 
  geom_nodepoint(
    aes(fill = UFboot),
    shape = 21,
    size = 2,
    color = 'white'
  ) +
  geom_tiplab(size = 1) +
  labs(title = 'OPRM') +
  scale_color_manual(values = c('black', 'indianred'),
                     guide = 'none') +
  scale_fill_viridis_c() +
  coord_cartesian(clip = "off") +
  NULL

save_plot(here('etc', 'plots', 'oprm_tree.pdf'),
          oprm_gg, base_width = 6, base_height = 10)

# SC6A3 --------------------------------------------------------------------

sc6a3_t <- read.iqtree(here('etc', 'data', 'trees', 'SC6A3_output', 'alignments', 'SC6A3_HUMAN.combined_final.aln.contree')) 
sc6a3_t %>% as_tibble() %>% View()

sc6a3_gg <- sc6a3_t %>% 
  groupClade(104) %>% 
  ggtree(
    aes(color = group),
    size = 0.2,
    # branch.length = 'none'
    # layout = 'circular'
  ) +
  # geom_nodelab(
  #   aes(label = node),
  #   size = 1) +
  geom_nodepoint(
    aes(fill = UFboot),
    shape = 21,
    size = 2,
    color = 'white'
  ) +
  geom_tiplab(size = 1) +
  labs(title = 'SC6A3') +
  scale_color_manual(values = c('black', 'indianred'),
                     guide = 'none') +
  scale_fill_viridis_c() +
  coord_cartesian(clip = "off") +
  NULL

save_plot(here('etc', 'plots', 'sc6a3_tree.pdf'),
          sc6a3_gg, base_width = 5, base_height = 8)


pdf(width = 12, height = 8, here('etc', 'supplementary', 'SupplementaryFigure5.pdf'))
aa3r_gg
arf6_gg
ht2a_gg
kcnh2_gg
oprm_gg
sc6a3_gg
dev.off()


