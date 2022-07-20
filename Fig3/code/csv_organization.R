library(tidyverse)
library(readxl)

compound_target_combos <- read_excel(here('Fig3', 'metadata', '7151 - Gene IDs - Aug 21.xlsx'), 
                                     skip = 1, na = 'N/A') %>% 
  janitor::clean_names() %>% 
  dplyr::select(-fda_approved) %>% 
  pivot_longer(cols = contains('target'), names_to = 'target_number', values_to = 'gene_name') %>% 
  mutate(target_number = str_remove(target_number, 'target_[0-9]*_')) %>% 
  drop_na() %>% 
  pivot_wider(names_from = target_number, values_from = gene_name) %>% 
  unnest(cols = c(gene_symbol, uni_prot)) %>% 
  rename(gene_name = gene_symbol, uniprot = uni_prot)

compound_info <- read_excel(here('Fig3', 'metadata', '7151 - Tocriscreen 2.0 Mini - CUSTOMER FILE - May 2021.xlsx'),
                            skip = 5, sheet = 'Compound Information') %>% 
  janitor::clean_names()

human_targets <- read_csv(here('Fig3', 'metadata', 'human_target_dataset.csv')) %>% 
  rename(gene_id = wbps_gene_id, gene_name = external_gene_id)

all <- left_join(human_targets, compound_target_combos) %>%
  left_join(compound_info) %>% 
  dplyr::select(gene_id:primary_action) %>% 
  write_csv(here('Fig3', 'metadata', 'compound_target_info.csv'))

