library(tidyverse)
library(here)


primary_hits <- read_csv(here("Fig2/tables/Table1.csv")) %>%
  mutate(treatment = str_remove(treatment, " hydrochloride| maleate| citrate"))

dose_response <- read_rds(here('Fig2/data/ec50.rds')) %>%
  mutate(treatment = str_remove(treatment, " hydrochloride| maleate| citrate"))

cel_hits <- read_rds(here('Fig7/data/cel_hits.rds')) %>%
  mutate(treatment = str_remove(treatment, " hydrochloride| maleate| citrate")) %>% 
  rename(cel = conc)

motility_hits <- read_rds(here("Fig4/data/motility_hits.rds")) %>%
  pivot_wider(id_cols = treatment, names_from = stages, values_from = hit) %>%
  mutate(adult_hit_type = case_when(
    `Adult female` == "Non-Hit" & `Adult male` == "Non-Hit" ~ "Non-Hit",
    `Adult female` == "Hit" & `Adult male` == "Hit" ~ "Both",
    `Adult female` == "Hit" ~ "Female",
    `Adult male` == "Hit" ~ "Male",
    TRUE ~ treatment
  )) %>%
  select(-contains("male"), adult_motility_hit_type = adult_hit_type)

fecundity_hits <- read_rds(here('Fig4/data/fecundity_hits.rds')) %>% 
  select(treatment, fecundity_mean = treatment_mean)
  
metabolism_hits <- read_rds(here('Fig4/data/metabolism_hits.rds')) %>% 
  select(stages, treatment, metabolism_mean = treatment_mean) %>% 
  pivot_wider(names_from = stages, values_from = metabolism_mean) %>% 
  rename('male_metabolism' = 'Adult male', 'female_metabolism' = 'Adult female')

final_table <- left_join(primary_hits, dose_response) %>%
  left_join(motility_hits) %>%
  left_join(fecundity_hits) %>%
  left_join(metabolism_hits) %>% 
  left_join(cel_hits) %>%
  select(-adult_motility_hit_type) %>%
  filter(treatment != "GNF 5837") %>% 
  rename(Number = number, Compound = treatment, 'Target class' = target_class,
         'Primary target' = primary_target, 'Primary action' = primary_action,
         'Z-motility' = z_motility, 'Z-viability' = z_viability, 'Z-total' = z_total,
         'Primary screen hit type' = hit_type, 'EC50 (12 hr. motility)' = motility_12hr, 'EC50 (36 hr. motility)' = motility_36hr, 
         'EC50 (36 hr. viability)' = viability_36hr, 'Fecundity mean (% control)' = fecundity_mean,
         'Female metabolism (% control)' = female_metabolism, 'Male metabolism (% control)' = male_metabolism,
         'C. elegans hit (concentration)' = cel)

write_csv(final_table, here('Fig2/tables/Table1_final.csv'))

cor.test(final_table$`Z-total`, final_table$`EC50 (36 hr. motility)`)


