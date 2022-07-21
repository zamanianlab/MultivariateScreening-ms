library(tidyverse)
library(ggtext)
library(ZamanianLabThemes)
library(here)

hits <- read_csv(here('Fig2', 'tables', 'Table1.csv'))

id_map <- read_tsv(here('Fig3', 'metadata', 'id_map.tab')) %>% 
  janitor::clean_names() %>% 
  select(-protein_names) %>% 
  mutate(gene_name = str_remove(gene_names, '\\s.*$'))

clusters <- read_csv(here('Fig3', 'metadata', 'human_clusters.csv')) %>% 
  janitor::clean_names() %>% 
  rename(target = uni_prot_kb)

truth <- read_csv(here('Fig3', 'metadata', 'compound_target_info.csv')) %>% 
  filter(!is.na(cat_no)) %>% 
  rename(molid = cat_no) %>% 
  left_join(id_map) %>% 
  left_join(clusters) %>% 
  rename(tocris_target_class = target_class,
         tocris_primary_target = primary_target,
         nw_target_class = cluster,
         target_name = entry_name) %>% 
  mutate(is_hit = compound_name %in% hits$treatment)

blastout <- read_tsv(here('Fig2', 'data', 'human_targets_vs_brugia.out'),
                     col_names = c('query', 'subject', 'length', 'percent_positive', 'percent_identity', 'evalue', 'bitscore')) %>% 
  group_by(query) %>% 
  arrange(-bitscore, .by_group = TRUE) %>% 
  slice_head(n = 1) %>% 
  mutate(query = str_remove(query, 'sp\\|.*\\|')) %>% 
  left_join(truth, by = c('query' = 'target_name')) %>% 
  filter(!is.na(gene_id))

blastout %>% 
  group_by(is_hit) %>% 
  summarise(mean_bitscore = mean(bitscore))

t.test(bitscore ~ is_hit, data = blastout)

summary <- blastout %>%
  group_by(is_hit) %>% 
  summarise(Mean = mean(bitscore),
            Median = median(bitscore)) %>% 
  pivot_longer(cols = c(Mean, Median), names_to = 'measurement', values_to = 'value') %>% 
  mutate(
    is_hit = case_when(
      is_hit == TRUE ~ 'Target of primary hit',
      is_hit == FALSE ~ 'Target of non-hit'
    ))

(hist <- blastout %>% 
    mutate(
      is_hit = case_when(
        is_hit == TRUE ~ 'Target of primary hit',
        is_hit == FALSE ~ 'Target of non-hit'
      )
    ) %>% 
    ggplot() +
    # geom_point(aes(x = n, y = 1 / evalue, color = is_hit)) +
    geom_density(aes(x = bitscore),
                 alpha = 0.6,
                 fill = 'black') +
    geom_vline(data = summary,
               aes(xintercept = value, linetype = measurement)) +
    geom_rug(aes(x = bitscore),
             sides = 't',
             color = 'indianred', size = 0.5, length = unit(0.05, "npc")) +
    labs(x = 'Bitscore', y = 'Density', linetype = 'Measurement') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.00275)) +
    facet_grid(rows = vars(is_hit)) +
    theme_nw2() +
    theme(
      legend.position = c(0.8, 0.8),
      axis.text.x = element_markdown(angle = 0, hjust = 0.5)
    ) +
    NULL
)

cowplot::save_plot(here('Fig2', 'supplementary', 'SupplementaryFig2.pdf'),
                   hist, base_width = 7)

tibble(query = blastout$query %in% truth$target_name) %>% group_by(query) %>% tally()



"Does the target have an ortholog in Brugia?
Is the ortholog highly similar?
Do hits have higher similarity than non-hits?"